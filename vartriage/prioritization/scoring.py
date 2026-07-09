"""Pathogenicity score normalization and composite ranking.

This module normalizes CADD Phred and REVEL scores to a common 0.0-1.0
scale, computes a weighted composite pathogenicity rank, and sorts
variants in descending order by that rank with nulls placed last.

Numpy vectorized operations are used for batch normalization when
processing 1,000+ variants to avoid iterative Python loops.
"""

from __future__ import annotations

import warnings
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from vartriage.exceptions import VarTriageWarning
from vartriage.models.variant import AnnotatedVariant, ScoredVariant
from vartriage.models.warnings import MissingDataWarning


REVEL_WEIGHT: float = 0.6
CADD_WEIGHT: float = 0.4
CADD_MAX_PHRED: float = 99.0


class ScoreValidationWarning(VarTriageWarning):
    """Emitted when a pathogenicity score is out of its valid range."""


def normalize_cadd_scores(
    scores: Sequence[Optional[float]],
) -> list[Optional[float]]:
    """Normalize CADD Phred scores to the 0.0-1.0 scale using vectorized ops.

    Parameters
    ----------
    scores : Sequence[Optional[float]]
        Raw CADD Phred scores. None entries represent missing data.

    Returns
    -------
    list[Optional[float]]
        Normalized scores where each valid score is divided by 99.0 and
        capped at 1.0. Negative scores are rejected (returned as None)
        and trigger a ScoreValidationWarning. None inputs remain None.
    """
    n = len(scores)
    if n == 0:
        return []

    result: list[Optional[float]] = [None] * n
    valid_indices: list[int] = []
    valid_values: list[float] = []

    for i, score in enumerate(scores):
        if score is None:
            continue
        if score < 0.0:
            warnings.warn(
                f"Negative CADD Phred score ({score}) at index {i}; "
                f"excluding from composite calculation.",
                ScoreValidationWarning,
                stacklevel=2,
            )
            continue
        valid_indices.append(i)
        valid_values.append(score)

    if valid_values:
        arr: NDArray[np.float64] = np.array(valid_values, dtype=np.float64)
        normalized = np.minimum(arr / CADD_MAX_PHRED, 1.0)
        for idx, norm_val in zip(valid_indices, normalized):
            result[idx] = float(norm_val)

    return result


def validate_revel_scores(
    scores: Sequence[Optional[float]],
) -> list[Optional[float]]:
    """Validate REVEL scores are within the 0.0-1.0 range.

    Parameters
    ----------
    scores : Sequence[Optional[float]]
        REVEL scores. None entries represent missing data.

    Returns
    -------
    list[Optional[float]]
        Validated scores. Out-of-range scores are rejected (returned as
        None) and trigger a ScoreValidationWarning. None inputs remain None.
    """
    n = len(scores)
    if n == 0:
        return []

    result: list[Optional[float]] = [None] * n

    for i, score in enumerate(scores):
        if score is None:
            continue
        if score < 0.0 or score > 1.0:
            warnings.warn(
                f"REVEL score ({score}) at index {i} outside valid range "
                f"[0.0, 1.0]; excluding from composite calculation.",
                ScoreValidationWarning,
                stacklevel=2,
            )
            continue
        result[i] = score

    return result


def compute_composite_ranks(
    cadd_normalized: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
) -> list[Optional[float]]:
    """Compute composite pathogenicity ranks from normalized scores.

    Uses the formula: (REVEL * 0.6) + (CADD_normalized * 0.4) when both
    are available. Falls back to the single available score when only one
    is present. Returns None when neither score is available.

    Parameters
    ----------
    cadd_normalized : Sequence[Optional[float]]
        Normalized CADD scores (0.0-1.0 scale).
    revel_scores : Sequence[Optional[float]]
        Validated REVEL scores (0.0-1.0 scale).

    Returns
    -------
    list[Optional[float]]
        Composite ranks for each variant.

    Raises
    ------
    ValueError
        If input sequences differ in length.
    """
    n = len(cadd_normalized)
    if n != len(revel_scores):
        raise ValueError(
            f"Score sequences must have equal length: "
            f"CADD has {n}, REVEL has {len(revel_scores)}"
        )

    if n == 0:
        return []

    cadd_arr = np.array(
        [c if c is not None else np.nan for c in cadd_normalized],
        dtype=np.float64,
    )
    revel_arr = np.array(
        [r if r is not None else np.nan for r in revel_scores],
        dtype=np.float64,
    )

    both_mask = ~np.isnan(cadd_arr) & ~np.isnan(revel_arr)
    cadd_only_mask = ~np.isnan(cadd_arr) & np.isnan(revel_arr)
    revel_only_mask = np.isnan(cadd_arr) & ~np.isnan(revel_arr)

    composite = np.full(n, np.nan, dtype=np.float64)
    composite[both_mask] = (
        revel_arr[both_mask] * REVEL_WEIGHT + cadd_arr[both_mask] * CADD_WEIGHT
    )
    composite[cadd_only_mask] = cadd_arr[cadd_only_mask]
    composite[revel_only_mask] = revel_arr[revel_only_mask]

    result: list[Optional[float]] = []
    for val in composite:
        if np.isnan(val):
            result.append(None)
        else:
            result.append(float(val))

    return result


def score_variants(
    variants: Sequence[AnnotatedVariant],
    cadd_scores: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
) -> list[ScoredVariant]:
    """Score a batch of annotated variants and sort by composite rank.

    Normalizes CADD Phred scores, validates REVEL scores, computes composite
    pathogenicity ranks, and sorts the result in descending order by
    composite_rank with null-ranked variants placed last.

    Emits MissingDataWarning for variants with no pathogenicity scores.

    Parameters
    ----------
    variants : Sequence[AnnotatedVariant]
        Annotated variants to score.
    cadd_scores : Sequence[Optional[float]]
        Raw CADD Phred scores aligned with the variants sequence.
    revel_scores : Sequence[Optional[float]]
        Raw REVEL scores aligned with the variants sequence.

    Returns
    -------
    list[ScoredVariant]
        Scored variants sorted descending by composite_rank, nulls last.

    Raises
    ------
    ValueError
        If input sequences differ in length.
    """
    n = len(variants)
    if n != len(cadd_scores) or n != len(revel_scores):
        raise ValueError(
            f"All input sequences must have equal length. "
            f"variants={n}, cadd_scores={len(cadd_scores)}, "
            f"revel_scores={len(revel_scores)}"
        )

    if n == 0:
        return []

    cadd_normalized = normalize_cadd_scores(cadd_scores)
    revel_validated = validate_revel_scores(revel_scores)
    composites = compute_composite_ranks(cadd_normalized, revel_validated)

    missing_data_warnings: list[MissingDataWarning] = []
    scored: list[ScoredVariant] = []

    for i, variant in enumerate(variants):
        raw_cadd = cadd_scores[i] if cadd_normalized[i] is not None else None
        composite = composites[i]

        if composite is None:
            warning = MissingDataWarning(
                chrom=variant.variant.chrom,
                pos=variant.variant.pos,
                ref=variant.variant.ref,
                alt=variant.variant.alt,
                source="pathogenicity_scores",
                reason="no CADD or REVEL scores available",
            )
            missing_data_warnings.append(warning)

        scored_variant = ScoredVariant(
            annotated=variant,
            cadd_phred=raw_cadd,
            cadd_normalized=cadd_normalized[i],
            revel_score=revel_validated[i],
            composite_rank=composite,
        )
        scored.append(scored_variant)

    for w in missing_data_warnings:
        warnings.warn(
            f"MissingDataWarning: {w.chrom}:{w.pos} {w.ref}>{w.alt} - "
            f"{w.reason}",
            UserWarning,
            stacklevel=2,
        )

    return sort_by_composite_rank(scored)


def sort_by_composite_rank(variants: list[ScoredVariant]) -> list[ScoredVariant]:
    """Sort scored variants descending by composite_rank, nulls last.

    Parameters
    ----------
    variants : list[ScoredVariant]
        Scored variants to sort.

    Returns
    -------
    list[ScoredVariant]
        Variants sorted with highest composite_rank first and null-ranked
        variants at the end.
    """
    def sort_key(v: ScoredVariant) -> tuple[int, float]:
        if v.composite_rank is None:
            return (1, 0.0)
        return (0, -v.composite_rank)

    return sorted(variants, key=sort_key)
