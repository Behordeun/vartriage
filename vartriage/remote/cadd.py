"""Remote tabix CADD Phred score backend.

Queries CADD scores from a remote bgzipped/tabix-indexed TSV using HTTP
byte-range requests via pysam.TabixFile. Only the compressed blocks covering
each query region are transferred — the full 80+ GB file stays on the server.

CADD TSV format (whole_genome_SNVs.tsv.gz):
    Column 1: Chrom (bare number, no "chr" prefix)
    Column 2: Pos (1-based)
    Column 3: Ref
    Column 4: Alt
    Column 5: RawScore
    Column 6: PHRED

Multi-allelic handling: CADD pre-computes all 3 possible substitutions per
position. We match on exact ref+alt to return the correct score.
"""

from __future__ import annotations

import contextlib
import logging
import time
from collections import defaultdict

import pysam

from vartriage.prioritization.score_loader import CoordinateKey
from vartriage.remote.cache import RemoteScoreCache
from vartriage.remote.circuit_breaker import CircuitBreaker
from vartriage.remote.config import RemoteTabixConfig
from vartriage.remote.presets import resolve_preset

logger = logging.getLogger(__name__)

_SOURCE_ID = "cadd-remote"


class RemoteTabixCADD:
    """CADD score lookups via remote tabix-indexed TSV.

    Opens a pysam.TabixFile pointed at a remote URL (or named preset).
    Batches nearby variants into range queries to reduce HTTP round-trips.
    Results are cached in SQLite to skip the network on repeated runs.

    Parameters
    ----------
    config : RemoteTabixConfig
        Remote backend configuration. Must have cadd_remote_url set.

    Raises
    ------
    ValueError
        If cadd_remote_url is None.
    OSError
        If the remote URL cannot be opened (connection failure).
    """

    def __init__(self, config: RemoteTabixConfig) -> None:
        if config.cadd_remote_url is None:
            raise ValueError("cadd_remote_url is required for RemoteTabixCADD")

        self._config = config
        self._url = resolve_preset(config.cadd_remote_url)
        self._cache = RemoteScoreCache(
            db_path=config.cache_path,
            ttl_days=config.cache_ttl_days,
        )
        self._breaker = CircuitBreaker()
        self._tabix: pysam.TabixFile | None = None

        # Stats for audit trail
        self._network_fetches = 0
        self._cache_hits = 0

    def lookup_batch(self, variants: list[CoordinateKey]) -> dict[CoordinateKey, float]:
        """Look up CADD Phred scores for a batch of variants.

        Checks the local cache first, then queries uncached variants
        via remote tabix. Nearby variants within batch_window_bp are
        grouped into single range queries.

        Parameters
        ----------
        variants : list[CoordinateKey]
            (chrom, pos, ref, alt) tuples.

        Returns
        -------
        dict[CoordinateKey, float]
            Scores keyed by coordinate. Variants without scores
            (not found or network failure) are absent.
        """
        if not variants:
            return {}

        results: dict[CoordinateKey, float] = {}
        uncached: list[CoordinateKey] = []

        # Phase 1: check cache
        cache_keys = list(variants)
        cached_scores = self._cache.get_batch(_SOURCE_ID, cache_keys)

        for i, score in enumerate(cached_scores):
            if score is not None:
                results[variants[i]] = score
                self._cache_hits += 1
            else:
                uncached.append(variants[i])

        if not uncached:
            return results

        # Phase 2: query remote for cache misses
        if self._breaker.is_open:
            logger.debug(
                "Circuit breaker open — skipping %d remote CADD queries",
                len(uncached),
            )
            return results

        remote_scores = self._query_remote_batched(uncached)
        results.update(remote_scores)

        # Phase 3: cache newly fetched scores
        if cache_entries := [
            (chrom, pos, ref, alt, score)
            for (chrom, pos, ref, alt), score in remote_scores.items()
        ]:
            self._cache.put_batch(_SOURCE_ID, cache_entries)

        return results

    def close(self) -> None:
        """Release resources."""
        if self._tabix is not None:
            self._tabix.close()
            self._tabix = None
        self._cache.close()

    @property
    def network_fetches(self) -> int:
        """Count of variants fetched from the remote server."""
        return self._network_fetches

    @property
    def cache_hits(self) -> int:
        """Count of cache hits."""
        return self._cache_hits

    def _query_remote_batched(
        self, variants: list[CoordinateKey]
    ) -> dict[CoordinateKey, float]:
        """Group nearby variants and query remote tabix.

        Variants within batch_window_bp on the same chromosome are
        merged into a single range query to reduce HTTP round-trips.
        """
        results: dict[CoordinateKey, float] = {}

        # Group by chromosome
        by_chrom: dict[str, list[CoordinateKey]] = defaultdict(list)
        for key in variants:
            by_chrom[key[0]].append(key)

        window = self._config.batch_window_bp

        for chrom, chrom_variants in by_chrom.items():
            # Sort by position for range grouping
            sorted_vars = sorted(chrom_variants, key=lambda k: k[1])

            # Group into windows
            groups = self._group_by_window(sorted_vars, window)

            for group in groups:
                group_results = self._query_range(chrom, group)
                results.update(group_results)

        return results

    def _group_by_window(
        self, sorted_variants: list[CoordinateKey], window: int
    ) -> list[list[CoordinateKey]]:
        """Split sorted variants into groups where consecutive positions
        are within `window` bp of each other."""
        if not sorted_variants:
            return []

        groups: list[list[CoordinateKey]] = []
        current_group: list[CoordinateKey] = [sorted_variants[0]]

        for variant in sorted_variants[1:]:
            if variant[1] - current_group[0][1] <= window:
                current_group.append(variant)
            else:
                groups.append(current_group)
                current_group = [variant]

        groups.append(current_group)
        return groups

    def _query_range(
        self, chrom: str, group: list[CoordinateKey]
    ) -> dict[CoordinateKey, float]:
        """Query a range of positions from the remote tabix file.

        CADD files use bare chromosome numbers (no "chr" prefix).
        We strip "chr" before querying and match results back using
        the original key format.
        """
        results: dict[CoordinateKey, float] = {}

        # Build a lookup set for fast matching
        wanted: dict[tuple[int, str, str], CoordinateKey] = {}
        for key in group:
            _, pos, ref, alt = key
            wanted[(pos, ref, alt)] = key

        # CADD uses bare chromosome numbers
        query_chrom = chrom.removeprefix("chr")

        start_pos = min(k[1] for k in group)
        end_pos = max(k[1] for k in group)

        try:
            tabix = self._get_tabix()
            records = tabix.fetch(query_chrom, start_pos - 1, end_pos)
        except (OSError, ValueError) as exc:
            self._breaker.record_failure()
            logger.warning(
                "Remote CADD query failed for %s:%d-%d: %s",
                query_chrom,
                start_pos,
                end_pos,
                exc,
            )
            return results

        self._breaker.record_success()

        for record_line in records:
            parsed = self._parse_cadd_record(record_line)
            if parsed is None:
                continue

            pos, ref, alt, phred = parsed
            lookup_key = (pos, ref, alt)

            if lookup_key in wanted:
                original_key = wanted[lookup_key]
                results[original_key] = phred
                self._network_fetches += 1

        return results

    def _parse_cadd_record(
        self, record_line: str
    ) -> tuple[int, str, str, float] | None:
        """Parse a CADD TSV record into (pos, ref, alt, phred).

        Expected columns: Chrom, Pos, Ref, Alt, RawScore, PHRED
        """
        fields = record_line.split("\t")
        if len(fields) < 6:
            return None

        try:
            pos = int(fields[1])
            ref = fields[2]
            alt = fields[3]
            phred = float(fields[5])
        except (ValueError, IndexError):
            return None

        return pos, ref, alt, phred

    def _get_tabix(self) -> pysam.TabixFile:
        """Lazy-open the remote tabix file."""
        if self._tabix is None:
            logger.info("Opening remote CADD tabix: %s", self._url)
            self._tabix = pysam.TabixFile(self._url)
        return self._tabix

    def _retry_query(self, query_chrom: str, start: int, end: int) -> list[str] | None:
        """Query with retries and exponential backoff.

        Returns list of record lines on success, None on exhausted retries.
        """
        max_retries = self._config.max_retries
        backoff = 1.0

        for attempt in range(max_retries + 1):
            try:
                tabix = self._get_tabix()
                return list(tabix.fetch(query_chrom, start, end))
            except (OSError, ValueError) as exc:
                if attempt == max_retries:
                    self._breaker.record_failure()
                    logger.warning(
                        "Remote CADD query exhausted retries for %s:%d-%d: %s",
                        query_chrom,
                        start,
                        end,
                        exc,
                    )
                    return None

                logger.debug(
                    "Remote CADD query attempt %d failed, retrying in %.1fs: %s",
                    attempt + 1,
                    backoff,
                    exc,
                )
                time.sleep(backoff)
                backoff *= 2.0

                # Reset connection on failure
                if self._tabix is not None:
                    with contextlib.suppress(Exception):
                        self._tabix.close()
                    self._tabix = None

        return None
