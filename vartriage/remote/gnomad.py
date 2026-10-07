"""Remote tabix gnomAD allele frequency backend.

Queries gnomAD allele frequencies from remote tabix-indexed VCF files
using HTTP byte-range requests via pysam.TabixFile. gnomAD distributes
one VCF per chromosome, so the URL template uses a {chrom} placeholder.

Satisfies the FrequencyDatabase protocol from vartriage.protocols.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import logging
import time
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import pysam

from vartriage.remote.cache import RemoteScoreCache
from vartriage.remote.circuit_breaker import CircuitBreaker
from vartriage.remote.config import RemoteTabixConfig
from vartriage.remote.presets import resolve_cache_source_id, resolve_preset

logger = logging.getLogger(__name__)

_SOURCE_ID = "gnomad-remote"

# gnomAD v4 ancestry groups carried as AF_<pop> subfields in the VCF INFO column,
# alongside the global "AF". These are the seven continental groups gnomAD reports
# per-ancestry allele frequencies for.
GNOMAD_POPULATIONS: tuple[str, ...] = (
    "afr",
    "amr",
    "asj",
    "eas",
    "fin",
    "nfe",
    "sas",
)

# INFO keys parsed for per-population lookups: global AF plus AF_<pop> for each group.
_POP_AF_KEYS: tuple[str, ...] = (
    "AF",
    *tuple(f"AF_{pop}" for pop in GNOMAD_POPULATIONS),
)

# Homozygote-count key carried alongside the AF subfields for BS2 evaluation.
# gnomAD v4 names it nhomalt; it is split per-ALT like the AF_* fields.
_HOM_COUNT_KEY = "nhomalt"

_VariantKey = tuple[str, int, str, str]


class RemoteTabixGnomAD:
    """gnomAD frequency lookups via remote tabix-indexed VCF.

    Satisfies the FrequencyDatabase protocol. Handles per-chromosome
    URL templates and multi-allelic VCF records. Uses retry with
    exponential backoff on transient network failures, consistent
    with the CADD backend.

    Parameters
    ----------
    config : RemoteTabixConfig
        Remote backend configuration. Must have gnomad_remote_url set.

    Raises
    ------
    ValueError
        If gnomad_remote_url is None.
    """

    def __init__(self, config: RemoteTabixConfig) -> None:
        if config.gnomad_remote_url is None:
            raise ValueError("gnomad_remote_url is required for RemoteTabixGnomAD")

        self._config = config
        self._url_template = resolve_preset(config.gnomad_remote_url)
        self._source_id = resolve_cache_source_id(_SOURCE_ID, config.gnomad_remote_url)
        self._cache = RemoteScoreCache(
            db_path=config.cache_path,
            ttl_days=config.cache_ttl_days,
        )
        self._breaker = CircuitBreaker()
        self._tabix_handles: dict[str, pysam.TabixFile] = {}

        # Stats for audit trail
        self._network_fetches = 0
        self._cache_hits = 0

    def load(self, reference_path: Path) -> None:
        """No-op for protocol compatibility.

        Remote backends don't load from a local file. The URL is
        configured via RemoteTabixConfig at construction time.
        """

    def lookup_batch(self, variants: list[_VariantKey]) -> list[float | None]:
        """Query allele frequencies for a batch of variants.

        Checks the local cache first, then queries uncached variants
        via remote tabix. Returns frequencies in the same positional
        order as the input list.

        Parameters
        ----------
        variants : list[tuple[str, int, str, str]]
            (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[float | None]
            Allele frequencies positionally matched. None for
            variants not found or when the circuit breaker is open.
        """
        if not variants:
            return []

        results: list[float | None] = [None] * len(variants)
        uncached_indices: list[int] = []

        # Phase 1: check cache. A score hit returns the frequency; a confirmed
        # absence returns None without a remote query. Only variants in neither
        # table are queried remotely.
        variant_list = list(variants)
        cached_scores = self._cache.get_batch(self._source_id, variant_list)
        absent_flags: list[bool] | None = None

        for i, score in enumerate(cached_scores):
            if score is not None:
                results[i] = score
                self._cache_hits += 1
            else:
                if absent_flags is None:
                    absent_flags = self._cache.get_absent_batch(
                        self._source_id, variant_list
                    )
                if absent_flags[i]:
                    self._cache_hits += 1
                else:
                    uncached_indices.append(i)

        if not uncached_indices:
            return results

        # Phase 2: query remote for cache misses
        if self._breaker.is_open:
            logger.debug(
                "Circuit breaker open — skipping %d remote gnomAD queries",
                len(uncached_indices),
            )
            return results

        uncached_variants = [variants[i] for i in uncached_indices]
        remote_results, confirmed_absent = self._query_remote_batched(uncached_variants)

        # Phase 3: merge results and cache. Present variants cache their score;
        # variants a successful fetch did not find cache as confirmed absent so
        # the next run is a hit, not a re-query. A failed fetch yields neither.
        cache_entries: list[tuple[str, int, str, str, float]] = []

        for idx, variant in zip(uncached_indices, uncached_variants, strict=True):
            af = remote_results.get(variant)
            if af is not None:
                results[idx] = af
                chrom, pos, ref, alt = variant
                cache_entries.append((chrom, pos, ref, alt, af))

        if cache_entries:
            self._cache.put_batch(self._source_id, cache_entries)
        if confirmed_absent:
            self._cache.put_absent_batch(self._source_id, confirmed_absent)

        return results

    def lookup_batch_populations(
        self, variants: list[_VariantKey]
    ) -> list[dict[str, float] | None]:
        """Query per-population allele frequencies for a batch of variants.

        Returns, positionally matched to the input, a map per variant holding the
        global "AF" plus each "AF_<pop>" subfield present (keys: "AF" and
        "AF_afr".."AF_sas"). None marks a variant not found or a query skipped
        because the circuit breaker is open.

        This is the ancestry-aware companion to lookup_batch. It is backed by a
        dedicated population-map cache (JSON per variant) with the same TTL and
        pinning semantics as the single-float cache, so repeat runs are cache
        hits rather than fresh remote queries. The global FrequencyDatabase
        protocol path (lookup_batch) is unchanged.
        """
        if not variants:
            return []

        results: list[dict[str, float] | None] = [None] * len(variants)
        uncached_indices: list[int] = []

        # Phase 1: check the population cache. A stored empty map is a hit
        # (confirmed queried, nothing found) and leaves the result at None.
        cached_maps = self._cache.get_population_batch(self._source_id, list(variants))
        for i, af_map in enumerate(cached_maps):
            if af_map is None:
                uncached_indices.append(i)
            elif af_map:
                results[i] = af_map

        if not uncached_indices:
            return results

        if self._breaker.is_open:
            logger.debug(
                "Circuit breaker open — skipping %d remote gnomAD population queries",
                len(uncached_indices),
            )
            return results

        # Phase 2+3: query remote for cache misses, persisting each group's
        # results as soon as they are fetched. Per-group persistence means an
        # interrupted run resumes from the last completed group rather than
        # losing a whole batch of work.
        uncached_variants = [variants[i] for i in uncached_indices]
        indices_by_variant: dict[_VariantKey, list[int]] = defaultdict(list)
        for idx in uncached_indices:
            indices_by_variant[variants[idx]].append(idx)

        for chrom, group in self._iter_groups(uncached_variants):
            fetch_ok, group_maps = self._query_range_populations(chrom, group)
            for variant, af_map in group_maps.items():
                for idx in indices_by_variant.get(variant, ()):
                    results[idx] = af_map

            # Only cache when the range query actually completed. A failed or
            # timed-out fetch returns no maps, which is indistinguishable from a
            # genuine absence by result alone — caching it as an empty map under
            # clinical pinning (cache_ttl_days=-1) would permanently record a
            # transiently-unreachable variant as confirmed-absent and never
            # re-query it, silently suppressing frequency evidence. Persist a
            # confirmed absence (empty map) only on a successful query.
            if not fetch_ok:
                continue
            group_entries: list[tuple[str, int, str, str, dict[str, float]]] = []
            for variant in group:
                chrom_g, pos_g, ref_g, alt_g = variant
                group_entries.append(
                    (chrom_g, pos_g, ref_g, alt_g, group_maps.get(variant, {}))
                )
            self._cache.put_population_batch(self._source_id, group_entries)

        return results

    def _query_range_populations(
        self, chrom: str, group: list[_VariantKey]
    ) -> tuple[bool, dict[_VariantKey, dict[str, float]]]:
        """Range-query the remote gnomAD VCF, returning per-population AF maps.

        Returns a (fetch_ok, maps) pair. fetch_ok is True when the remote range
        query completed (including a successful query that found nothing), and
        False when the fetch exhausted its retries or degraded. The caller must
        not cache results from a failed fetch, since an empty map then means
        "query failed", not "confirmed absent".
        """
        results: dict[_VariantKey, dict[str, float]] = {}

        wanted: dict[tuple[int, str, str], _VariantKey] = {}
        for variant in group:
            _, pos, ref, alt = variant
            wanted[(pos, ref, alt)] = variant

        start_pos = min(v[1] for v in group)
        end_pos = max(v[1] for v in group)
        query_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"

        records = self._fetch_records(chrom, query_chrom, start_pos - 1, end_pos)
        if records is None:
            return False, results

        for record_line in records:
            parsed = self._parse_gnomad_record_populations(record_line)
            if parsed is None:
                continue
            for pos, ref, alt, af_map in parsed:
                lookup_key = (pos, ref, alt)
                if lookup_key in wanted:
                    results[wanted[lookup_key]] = af_map
                    self._network_fetches += 1

        return True, results

    def close(self) -> None:
        """Release resources."""
        for handle in self._tabix_handles.values():
            with contextlib.suppress(Exception):
                handle.close()
        self._tabix_handles.clear()
        self._cache.close()

    @property
    def network_fetches(self) -> int:
        """Count of variants fetched from the remote server."""
        return self._network_fetches

    @property
    def cache_hits(self) -> int:
        """Count of cache hits."""
        return self._cache_hits

    # ------------------------------------------------------------------
    # Batching and grouping
    # ------------------------------------------------------------------

    def _query_remote_batched(
        self, variants: list[_VariantKey]
    ) -> tuple[dict[_VariantKey, float], list[_VariantKey]]:
        """Group variants by chromosome, batch by window, query remote.

        Returns a (results, confirmed_absent) pair. confirmed_absent lists the
        variants a successful range query did not find, so the caller can record
        them as cache hits rather than re-querying every run. A group whose fetch
        failed contributes nothing to either collection, so a transient failure
        is never mistaken for an absence.
        """
        results: dict[_VariantKey, float] = {}
        confirmed_absent: list[_VariantKey] = []
        for chrom, group in self._iter_groups(variants):
            fetch_ok, group_results = self._query_range(chrom, group)
            results.update(group_results)
            if fetch_ok:
                confirmed_absent.extend(v for v in group if v not in group_results)
        return results, confirmed_absent

    def _iter_groups(
        self, variants: list[_VariantKey]
    ) -> Iterator[tuple[str, list[_VariantKey]]]:
        """Yield (chrom, group) for variants grouped by chromosome and window."""
        by_chrom: dict[str, list[_VariantKey]] = defaultdict(list)
        for variant in variants:
            by_chrom[variant[0]].append(variant)

        window = self._config.batch_window_bp
        for chrom, chrom_variants in by_chrom.items():
            sorted_vars = sorted(chrom_variants, key=lambda v: v[1])
            for group in self._group_by_window(sorted_vars, window):
                yield chrom, group

    @staticmethod
    def _group_by_window(
        sorted_variants: list[_VariantKey],
        window: int,
    ) -> list[list[_VariantKey]]:
        """Split sorted variants into groups within `window` bp."""
        if not sorted_variants:
            return []

        groups: list[list[_VariantKey]] = []
        current_group: list[_VariantKey] = [sorted_variants[0]]

        for variant in sorted_variants[1:]:
            if variant[1] - current_group[-1][1] <= window:
                current_group.append(variant)
            else:
                groups.append(current_group)
                current_group = [variant]

        groups.append(current_group)
        return groups

    # ------------------------------------------------------------------
    # Range query with retry
    # ------------------------------------------------------------------

    def _query_range(
        self, chrom: str, group: list[_VariantKey]
    ) -> tuple[bool, dict[_VariantKey, float]]:
        """Query a range from the remote gnomAD VCF.

        Returns a (fetch_ok, results) pair. fetch_ok is True when the remote
        range query completed (including a successful query that found nothing),
        and False when the fetch exhausted its retries or degraded. The caller
        must not record absences from a failed fetch, since an empty result then
        means "query failed", not "confirmed absent".
        """
        results: dict[_VariantKey, float] = {}

        wanted: dict[tuple[int, str, str], _VariantKey] = {}
        for variant in group:
            _, pos, ref, alt = variant
            wanted[(pos, ref, alt)] = variant

        start_pos = min(v[1] for v in group)
        end_pos = max(v[1] for v in group)

        # gnomAD VCFs use "chr" prefix for GRCh38
        query_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"

        records = self._fetch_records(chrom, query_chrom, start_pos - 1, end_pos)
        if records is None:
            return False, results

        for record_line in records:
            parsed = self._parse_gnomad_record(record_line)
            if parsed is None:
                continue

            for pos, ref, alt, af in parsed:
                lookup_key = (pos, ref, alt)
                if lookup_key in wanted:
                    original_variant = wanted[lookup_key]
                    results[original_variant] = af
                    self._network_fetches += 1

        return True, results

    def _fetch_records(
        self, chrom: str, query_chrom: str, start: int, end: int
    ) -> list[str] | None:
        """Fetch records with retries, exponential backoff, and circuit breaker.

        Returns list of record lines on success, None on exhausted retries.
        Records breaker success on first successful fetch.
        Uses a thread-based timeout to prevent indefinite hangs on stalled
        S3 connections (pysam/htslib has no built-in timeout).
        """
        max_retries = self._config.max_retries
        timeout = self._config.read_timeout
        backoff = 1.0

        for attempt in range(max_retries + 1):
            try:
                tabix = self._get_tabix_for_chrom(chrom)

                pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                try:

                    def _do_fetch(
                        _t: pysam.TabixFile = tabix,
                        _c: str = query_chrom,
                        _s: int = start,
                        _e: int = end,
                    ) -> list[str]:
                        return list(_t.fetch(_c, _s, _e))

                    future = pool.submit(_do_fetch)
                    records = future.result(timeout=timeout)
                finally:
                    pool.shutdown(wait=False, cancel_futures=True)

                self._breaker.record_success()
                return records
            except concurrent.futures.TimeoutError:
                logger.warning(
                    "Remote gnomAD query timed out after %.0fs for %s:%d-%d "
                    "(attempt %d/%d)",
                    timeout,
                    chrom,
                    start,
                    end,
                    attempt + 1,
                    max_retries + 1,
                )
                self._reset_chrom_connection(chrom)
                if attempt == max_retries:
                    self._breaker.record_failure()
                    return None
                time.sleep(backoff)
                backoff *= 2.0
            except (OSError, ValueError, pysam.utils.SamtoolsError) as exc:
                if attempt == max_retries:
                    self._breaker.record_failure()
                    logger.warning(
                        "Remote gnomAD query exhausted retries for %s:%d-%d: %s",
                        chrom,
                        start,
                        end,
                        exc,
                    )
                    return None

                logger.debug(
                    "Remote gnomAD query attempt %d failed, retrying in %.1fs: %s",
                    attempt + 1,
                    backoff,
                    exc,
                )
                time.sleep(backoff)
                backoff *= 2.0
                self._reset_chrom_connection(chrom)
            except Exception as exc:
                # htslib surfaces transient S3/BGZF read failures (libcurl socket
                # errors, truncated BGZF blocks) as exceptions outside the OSError
                # hierarchy. A single bad range read must never abort a
                # genome-wide run, so treat any remote read failure the same as a
                # retryable fetch error: back off and retry, then degrade to no
                # frequency for this range.
                if attempt == max_retries:
                    self._breaker.record_failure()
                    logger.warning(
                        "Remote gnomAD query failed (non-standard error) for "
                        "%s:%d-%d: %s",
                        chrom,
                        start,
                        end,
                        exc,
                    )
                    return None

                logger.debug(
                    "Remote gnomAD query attempt %d hit a non-standard error, "
                    "retrying in %.1fs: %s",
                    attempt + 1,
                    backoff,
                    exc,
                )
                time.sleep(backoff)
                backoff *= 2.0
                self._reset_chrom_connection(chrom)

        return None  # unreachable, satisfies type checker

    # ------------------------------------------------------------------
    # Parsing and connection management
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_gnomad_record(
        record_line: str,
    ) -> list[tuple[int, str, str, float]] | None:
        """Parse a gnomAD VCF record into list of (pos, ref, alt, af).

        Handles multi-allelic records by splitting ALT and AF fields
        and returning one entry per alternate allele.
        """
        fields = record_line.split("\t")
        if len(fields) < 8:
            return None

        try:
            pos = int(fields[1])
        except ValueError:
            return None

        ref = fields[3]
        alts = fields[4].split(",")
        info_field = fields[7]

        af_str = _extract_info_field(info_field, "AF")
        if af_str is None:
            return None

        af_values = af_str.split(",")
        entries: list[tuple[int, str, str, float]] = []

        for i, alt in enumerate(alts):
            if i >= len(af_values):
                break
            try:
                af = float(af_values[i])
                entries.append((pos, ref, alt, af))
            except ValueError:
                continue

        return entries or None

    @staticmethod
    def _parse_gnomad_record_populations(
        record_line: str,
    ) -> list[tuple[int, str, str, dict[str, float]]] | None:
        """Parse a gnomAD VCF record into (pos, ref, alt, per-population AF map).

        Returns one entry per alternate allele. The map holds the global "AF"
        plus each "AF_<pop>" subfield present for that allele and, when present,
        the "nhomalt" homozygote count (carried as a float, read back as an int
        for BS2). Missing or malformed subfields are omitted rather than
        defaulted, so a caller can tell absent from zero. Multi-allelic records
        split every field on comma and index by ALT position, matching the
        global-AF parser.
        """
        fields = record_line.split("\t")
        if len(fields) < 8:
            return None

        try:
            pos = int(fields[1])
        except ValueError:
            return None

        ref = fields[3]
        alts = fields[4].split(",")
        info_field = fields[7]

        per_key_values: dict[str, list[str]] = {}
        for key in (*_POP_AF_KEYS, _HOM_COUNT_KEY):
            raw = _extract_info_field(info_field, key)
            if raw is not None:
                per_key_values[key] = raw.split(",")

        if "AF" not in per_key_values:
            return None

        entries: list[tuple[int, str, str, dict[str, float]]] = []
        for i, alt in enumerate(alts):
            af_map: dict[str, float] = {}
            for key, values in per_key_values.items():
                if i >= len(values):
                    continue
                token = values[i]
                if token in (".", ""):
                    continue
                try:
                    af_map[key] = float(token)
                except ValueError:
                    continue
            if "AF" in af_map:
                entries.append((pos, ref, alt, af_map))

        return entries or None

    def _get_tabix_for_chrom(self, chrom: str) -> pysam.TabixFile:
        """Get or open the tabix handle for a chromosome.

        gnomAD distributes per-chromosome VCFs. The URL template uses
        a {chrom} placeholder that gets formatted per chromosome.
        """
        if chrom in self._tabix_handles:
            return self._tabix_handles[chrom]

        url = self._url_template.format(chrom=chrom)
        logger.info("Opening remote gnomAD tabix for %s: %s", chrom, url)

        # pysam.TabixFile opens the remote index over HTTP with no connect
        # timeout (htslib has none), so a stalled S3 connection would block the
        # whole run forever. Bound the open with the configured connect timeout;
        # a timeout raises TimeoutError into the caller's retry/degrade path.
        #
        # The open runs in a C call that cannot be cancelled, so on timeout the
        # worker thread stays alive until htslib finally returns or errors. A
        # done-callback closes any handle that arrives after we have given up,
        # so the connection is reclaimed when the stalled open eventually
        # completes rather than leaking for the life of the process. The pool is
        # shut down without waiting so the caller is not blocked by the stall.
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        future = pool.submit(pysam.TabixFile, url)
        try:
            handle: pysam.TabixFile = future.result(
                timeout=self._config.connect_timeout
            )
        except BaseException:
            future.add_done_callback(_close_late_tabix_handle)
            raise
        finally:
            pool.shutdown(wait=False, cancel_futures=True)

        self._tabix_handles[chrom] = handle
        return handle

    def _reset_chrom_connection(self, chrom: str) -> None:
        """Close and discard the tabix handle for a chromosome."""
        handle = self._tabix_handles.pop(chrom, None)
        if handle is not None:
            with contextlib.suppress(Exception):
                handle.close()


def _close_late_tabix_handle(
    future: concurrent.futures.Future[pysam.TabixFile],
) -> None:
    """Close a TabixFile handle that opened after its wait timed out.

    When a remote index open exceeds the connect timeout, the caller abandons
    the wait but the htslib open cannot be cancelled and keeps running. This
    callback runs when that open finally resolves: it closes a successfully
    opened handle so the connection is released, and swallows the stored error
    of a failed open. Without it, every stalled open leaks an open connection
    for the life of the process.
    """
    if future.cancelled():
        return
    exc = future.exception()
    if exc is not None:
        return
    handle = future.result()
    with contextlib.suppress(Exception):
        handle.close()


def _extract_info_field(info: str, key: str) -> str | None:
    """Extract a key's value from the VCF INFO column."""
    prefix = f"{key}="
    return next(
        (e[len(prefix) :] for e in info.split(";") if e.startswith(prefix)),
        None,
    )
