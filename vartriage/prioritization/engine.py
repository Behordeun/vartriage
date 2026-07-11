"""Prioritization engine orchestrator.

Composes allele frequency filtering with pathogenicity score normalization
and composite ranking to produce a prioritized stream of scored variants.
Processes variants in configurable batches and implements chunked fallback
on MemoryError to handle large datasets gracefully.
"""

from __future__ import annotations

import logging
from itertools import islice
from typing import Iterator

from vartriage.models.config import PrioritizationConfig
from vartriage.models.variant import AnnotatedVariant, ScoredVariant
from vartriage.prioritization.frequency_filter import FrequencyFilter
from vartriage.prioritization.score_loader import CoordinateKey, ScoreLoader
from vartriage.prioritization.scoring import score_variants

logger = logging.getLogger(__name__)

_MAX_CHUNK_SIZE: int = 500_000


class PrioritizationEngine:
    """Filter and rank variants by frequency and pathogenicity scores.

    Processes an iterator of annotated variants through two stages:

    1. **Frequency filtering**: excludes variants with allele frequency
       above the configured maximum threshold (retaining frequency-unknown
       variants).
    2. **Pathogenicity scoring**: normalizes CADD/REVEL scores, computes
       a composite rank, and sorts each batch in descending order by
       composite rank (nulls last).

    Variants are processed in configurable batches to bound memory usage.
    On MemoryError, the engine retries with a reduced chunk size (capped at
    500,000 variants per chunk).

    Parameters
    ----------
    config : PrioritizationConfig, optional
        Configuration containing ``max_allele_frequency``, ``batch_size``,
        and optional score file paths. When None, defaults are used
        (max_af=0.01, batch_size=10,000).

    Raises
    ------
    ValueError
        If ``config.max_allele_frequency`` is outside [0.0, 1.0] or
        ``config.batch_size`` is outside [1,000, 100,000]. Enforced at
        config construction time via ``PrioritizationConfig.__post_init__``.
    """

    def __init__(self, config: PrioritizationConfig | None = None) -> None:
        if config is None:
            config = PrioritizationConfig()
        self._config = config
        self._batch_size = config.batch_size
        self._frequency_filter = FrequencyFilter(config)
        self._score_loader = ScoreLoader()
        self._cadd_scores: dict[CoordinateKey, float] = {}
        self._revel_scores: dict[CoordinateKey, float] = {}
        self._spliceai_scores: dict[CoordinateKey, float] = {}

        if config.cadd_scores_path is not None:
            self._cadd_scores = self._score_loader.load_cadd(
                config.cadd_scores_path
            )
        if config.revel_scores_path is not None:
            self._revel_scores = self._score_loader.load_revel(
                config.revel_scores_path
            )
        if config.spliceai_scores_path is not None:
            self._spliceai_scores = self._score_loader.load_spliceai(
                config.spliceai_scores_path
            )

    def prioritize(
        self, variants: Iterator[AnnotatedVariant]
    ) -> Iterator[ScoredVariant]:
        """Filter by allele frequency and score remaining variants.

        Parameters
        ----------
        variants : Iterator[AnnotatedVariant]
            Input stream of annotated variant records. May be empty.

        Yields
        ------
        ScoredVariant
            Variants that pass frequency filtering, scored and sorted in
            descending order by composite pathogenicity rank within each
            batch. Variants with null composite rank appear last.
        """
        filtered = self._frequency_filter.apply(variants)
        yield from self._process_in_batches(filtered)

    def _process_in_batches(
        self, variants: Iterator[AnnotatedVariant]
    ) -> Iterator[ScoredVariant]:
        """Score variants in configurable batches.

        Parameters
        ----------
        variants : Iterator[AnnotatedVariant]
            Filtered variant stream.

        Yields
        ------
        ScoredVariant
            Scored variants from each batch, sorted within batch.
        """
        batch_size = self._batch_size

        while True:
            batch = list(islice(variants, batch_size))
            if not batch:
                break

            try:
                scored = self._score_batch(batch)
            except MemoryError:
                logger.warning(
                    "MemoryError during scoring of batch (size=%d). "
                    "Falling back to chunked processing.",
                    len(batch),
                )
                scored = self._chunked_fallback(batch)

            yield from scored

    def _score_batch(
        self, batch: list[AnnotatedVariant]
    ) -> list[ScoredVariant]:
        """Score a single batch of variants.

        Extracts coordinate keys from each variant and performs lookups
        against pre-loaded CADD and REVEL score dictionaries. Falls back
        to None for variants without a matching score entry.

        Parameters
        ----------
        batch : list[AnnotatedVariant]
            Batch of annotated variants to score.

        Returns
        -------
        list[ScoredVariant]
            Scored variants sorted descending by composite rank, nulls last.
        """
        keys: list[CoordinateKey] = [
            (v.variant.chrom, v.variant.pos, v.variant.ref, v.variant.alt)
            for v in batch
        ]
        cadd_scores = self._score_loader.lookup_batch(
            keys, self._cadd_scores
        )
        revel_scores = self._score_loader.lookup_batch(
            keys, self._revel_scores
        )

        spliceai_scores = None
        if self._spliceai_scores:
            spliceai_scores = self._score_loader.lookup_batch(
                keys, self._spliceai_scores
            )

        return score_variants(
            batch, cadd_scores, revel_scores, spliceai_scores
        )

    def _chunked_fallback(
        self, batch: list[AnnotatedVariant]
    ) -> list[ScoredVariant]:
        """Process a batch in smaller chunks after a MemoryError.

        Splits the batch into chunks of at most ``_MAX_CHUNK_SIZE``
        (500,000) variants and scores each chunk independently. Results
        from all chunks are merged and re-sorted by composite rank.

        Parameters
        ----------
        batch : list[AnnotatedVariant]
            The full batch that triggered MemoryError.

        Returns
        -------
        list[ScoredVariant]
            All scored variants merged and sorted.
        """
        chunk_size = min(_MAX_CHUNK_SIZE, max(1, len(batch) // 2))
        all_scored: list[ScoredVariant] = []

        for start in range(0, len(batch), chunk_size):
            chunk = batch[start : start + chunk_size]
            try:
                scored_chunk = self._score_batch(chunk)
                all_scored.extend(scored_chunk)
            except MemoryError:
                logger.error(
                    "MemoryError persists even with chunk_size=%d. "
                    "Reducing further.",
                    chunk_size,
                )
                smaller_size = max(1, chunk_size // 2)
                for sub_start in range(0, len(chunk), smaller_size):
                    sub_chunk = chunk[sub_start : sub_start + smaller_size]
                    scored_sub = self._score_batch(sub_chunk)
                    all_scored.extend(scored_sub)

        return _merge_sort_scored(all_scored)


def _merge_sort_scored(variants: list[ScoredVariant]) -> list[ScoredVariant]:
    """Sort scored variants descending by composite_rank, nulls last.

    Parameters
    ----------
    variants : list[ScoredVariant]
        Variants to sort.

    Returns
    -------
    list[ScoredVariant]
        Sorted variants.
    """

    def sort_key(v: ScoredVariant) -> tuple[int, float]:
        if v.composite_rank is None:
            return (1, 0.0)
        return (0, -v.composite_rank)

    return sorted(variants, key=sort_key)
