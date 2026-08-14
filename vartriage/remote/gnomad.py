"""Remote tabix gnomAD allele frequency backend.

Queries gnomAD allele frequencies from remote tabix-indexed VCF files
using HTTP byte-range requests via pysam.TabixFile. gnomAD distributes
one VCF per chromosome, so the URL template uses a {chrom} placeholder.

Satisfies the FrequencyDatabase protocol from vartriage.protocols.
"""

from __future__ import annotations

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
from vartriage.remote.presets import resolve_preset

logger = logging.getLogger(__name__)

_SOURCE_ID = "gnomad-remote"

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

        # Phase 1: check cache
        cached_scores = self._cache.get_batch(_SOURCE_ID, list(variants))

        for i, score in enumerate(cached_scores):
            if score is not None:
                results[i] = score
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
        remote_results = self._query_remote_batched(uncached_variants)

        # Phase 3: merge results and cache
        cache_entries: list[tuple[str, int, str, str, float]] = []

        for idx, variant in zip(uncached_indices, uncached_variants, strict=True):
            af = remote_results.get(variant)
            if af is not None:
                results[idx] = af
                chrom, pos, ref, alt = variant
                cache_entries.append((chrom, pos, ref, alt, af))

        if cache_entries:
            self._cache.put_batch(_SOURCE_ID, cache_entries)

        return results

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
    ) -> dict[_VariantKey, float]:
        """Group variants by chromosome, batch by window, query remote."""
        results: dict[_VariantKey, float] = {}
        for chrom, group in self._iter_groups(variants):
            group_results = self._query_range(chrom, group)
            results.update(group_results)
        return results

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
    ) -> dict[_VariantKey, float]:
        """Query a range from the remote gnomAD VCF."""
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
            return results

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

        return results

    def _fetch_records(
        self, chrom: str, query_chrom: str, start: int, end: int
    ) -> list[str] | None:
        """Fetch records with retries, exponential backoff, and circuit breaker.

        Returns list of record lines on success, None on exhausted retries.
        Records breaker success on first successful fetch.
        """
        max_retries = self._config.max_retries
        backoff = 1.0

        for attempt in range(max_retries + 1):
            try:
                tabix = self._get_tabix_for_chrom(chrom)
                records = list(tabix.fetch(query_chrom, start, end))
                self._breaker.record_success()
                return records
            except (OSError, ValueError) as exc:
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

    def _get_tabix_for_chrom(self, chrom: str) -> pysam.TabixFile:
        """Get or open the tabix handle for a chromosome.

        gnomAD distributes per-chromosome VCFs. The URL template uses
        a {chrom} placeholder that gets formatted per chromosome.
        """
        if chrom in self._tabix_handles:
            return self._tabix_handles[chrom]

        url = self._url_template.format(chrom=chrom)
        logger.info("Opening remote gnomAD tabix for %s: %s", chrom, url)
        handle = pysam.TabixFile(url)
        self._tabix_handles[chrom] = handle
        return handle

    def _reset_chrom_connection(self, chrom: str) -> None:
        """Close and discard the tabix handle for a chromosome."""
        handle = self._tabix_handles.pop(chrom, None)
        if handle is not None:
            with contextlib.suppress(Exception):
                handle.close()


def _extract_info_field(info: str, key: str) -> str | None:
    """Extract a key's value from the VCF INFO column."""
    prefix = f"{key}="
    return next(
        (e[len(prefix) :] for e in info.split(";") if e.startswith(prefix)),
        None,
    )
