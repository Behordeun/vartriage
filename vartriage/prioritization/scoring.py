"""Pathogenicity score normalization and composite ranking.

This module normalizes CADD Phred, REVEL, and SpliceAI scores to a
common 0.0-1.0 scale, computes a weighted composite pathogenicity rank,
and sorts variants in descending order by that rank with nulls placed
last.

Numpy vectorized operations are used for batch normalization when
processing 1,000+ variants to avoid iterative Python loops.
"""

from __future__ import annotations

import warnings
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from vartriage.exceptions import VarTriageWarning
from vartriage.models.variant import (
    AnnotatedVariant,
    FunctionalConsequence,
    ScoredVariant,
)
from vartriage.models.warnings import MissingDataWarning

REVEL_WEIGHT: float = 0.6
CADD_WEIGHT: float = 0.4
CADD_MAX_PHRED: float = 99.0

REVEL_WEIGHT_3: float = 0.5
CADD_WEIGHT_3: float = 0.3
SPLICEAI_WEIGHT_3: float = 0.2


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


def validate_spliceai_scores(
    scores: Sequence[Optional[float]],
) -> list[Optional[float]]:
    """Validate SpliceAI scores are within the 0.0-1.0 range.

    Parameters
    ----------
    scores : Sequence[Optional[float]]
        SpliceAI scores. None entries represent missing data.

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
                f"SpliceAI score ({score}) at index {i} outside valid range "
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
    spliceai_scores: Optional[Sequence[Optional[float]]] = None,
) -> list[Optional[float]]:
    """Compute composite pathogenicity ranks from normalized scores.

    When ``spliceai_scores`` is None (SpliceAI not configured), the legacy
    two-score formula is used: REVEL * 0.6 + CADD * 0.4 when both are
    present, with single-score fallback otherwise.

    When ``spliceai_scores`` is provided as a sequence, dynamic proportional
    weight redistribution is applied based on the base 3-score weights
    (REVEL 0.5, CADD 0.3, SpliceAI 0.2).

    Parameters
    ----------
    cadd_normalized : Sequence[Optional[float]]
        Normalized CADD scores (0.0-1.0 scale).
    revel_scores : Sequence[Optional[float]]
        Validated REVEL scores (0.0-1.0 scale).
    spliceai_scores : Optional[Sequence[Optional[float]]]
        Validated SpliceAI scores (0.0-1.0 scale). None means SpliceAI
        is not configured and the legacy formula should be used.

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

    if spliceai_scores is not None and n != len(spliceai_scores):
        raise ValueError(
            f"Score sequences must have equal length: "
            f"CADD has {n}, SpliceAI has {len(spliceai_scores)}"
        )

    if n == 0:
        return []

    # Legacy path: SpliceAI not configured, use original 0.6/0.4 formula
    if spliceai_scores is None:
        return _compute_legacy_ranks(cadd_normalized, revel_scores)

    # 3-score path: dynamic proportional weight redistribution
    return _compute_three_score_ranks(cadd_normalized, revel_scores, spliceai_scores)


def _compute_legacy_ranks(
    cadd_normalized: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
) -> list[Optional[float]]:
    """Legacy two-score composite rank using REVEL*0.6 + CADD*0.4."""
    n = len(cadd_normalized)

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


# Base weights for three-score composite: (score_weight, score_index)
_THREE_SCORE_WEIGHTS: tuple[tuple[float, int], ...] = (
    (REVEL_WEIGHT_3, 0),
    (CADD_WEIGHT_3, 1),
    (SPLICEAI_WEIGHT_3, 2),
)


def _rank_three_scores(
    revel: Optional[float], cadd: Optional[float], splice: Optional[float]
) -> Optional[float]:
    """Compute composite rank for one variant with proportional weight redistribution."""
    scores = (revel, cadd, splice)
    present = [(w, s) for (w, _), s in zip(_THREE_SCORE_WEIGHTS, scores) if s is not None]
    if not present:
        return None
    if len(present) == 1:
        return present[0][1]
    weight_sum = sum(w for w, _ in present)
    return float(sum(s * (w / weight_sum) for w, s in present))


def _compute_three_score_ranks(
    cadd_normalized: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
    spliceai_scores: Sequence[Optional[float]],
) -> list[Optional[float]]:
    """Three-score composite rank with dynamic proportional weight redistribution.

    Base weights: REVEL=0.5, CADD=0.3, SpliceAI=0.2.
    When one or more scores are missing for a variant, the available scores'
    base weights are rescaled to sum to 1.0.
    """
    return [
        _rank_three_scores(revel_scores[i], cadd_normalized[i], spliceai_scores[i])
        for i in range(len(cadd_normalized))
    ]


def score_variants(
    variants: Sequence[AnnotatedVariant],
    cadd_scores: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
    spliceai_scores: Optional[Sequence[Optional[float]]] = None,
) -> list[ScoredVariant]:
    """Score a batch of annotated variants and sort by composite rank.

    Normalizes CADD Phred scores, validates REVEL and SpliceAI scores,
    computes composite pathogenicity ranks, and sorts the result in
    descending order by composite_rank with null-ranked variants placed last.

    Emits MissingDataWarning for variants with no pathogenicity scores.

    Parameters
    ----------
    variants : Sequence[AnnotatedVariant]
        Annotated variants to score.
    cadd_scores : Sequence[Optional[float]]
        Raw CADD Phred scores aligned with the variants sequence.
    revel_scores : Sequence[Optional[float]]
        Raw REVEL scores aligned with the variants sequence.
    spliceai_scores : Optional[Sequence[Optional[float]]]
        Raw SpliceAI scores aligned with the variants sequence. None
        means SpliceAI is not configured (legacy mode).

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

    if spliceai_scores is not None and n != len(spliceai_scores):
        raise ValueError(
            f"All input sequences must have equal length. "
            f"variants={n}, spliceai_scores={len(spliceai_scores)}"
        )

    if n == 0:
        return []

    cadd_normalized = normalize_cadd_scores(cadd_scores)
    revel_validated = validate_revel_scores(revel_scores)

    spliceai_validated: Optional[list[Optional[float]]] = None
    if spliceai_scores is not None:
        spliceai_validated = validate_spliceai_scores(spliceai_scores)

    composites = compute_composite_ranks(
        cadd_normalized, revel_validated, spliceai_validated
    )

    missing_data_warnings: list[MissingDataWarning] = []
    scored: list[ScoredVariant] = []

    for i, variant in enumerate(variants):
        raw_cadd = cadd_scores[i] if cadd_normalized[i] is not None else None
        composite = composites[i]

        splice_score: Optional[float] = None
        if spliceai_validated is not None:
            splice_score = spliceai_validated[i]

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
            spliceai_score=splice_score,
            composite_rank=composite,
            prioritization_score=compute_prioritization_score(
                consequence=variant.consequence,
                revel_score=revel_validated[i],
                spliceai_score=splice_score,
                cadd_phred=raw_cadd,
            ),
        )
        scored.append(scored_variant)

    for w in missing_data_warnings:
        warnings.warn(
            f"MissingDataWarning: {w.chrom}:{w.pos} {w.ref}>{w.alt} - " f"{w.reason}",
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


def compute_prioritization_score(
    consequence: "FunctionalConsequence",
    revel_score: Optional[float],
    spliceai_score: Optional[float],
    cadd_phred: Optional[float],
) -> Optional[float]:
    """Compute a prioritization score using validated, literature-backed logic.

    Strategy:
    - Missense: use REVEL directly (validated 0.7 threshold for PP3)
    - Splice-adjacent: use SpliceAI delta directly (validated 0.5 threshold)
    - Other coding: use CADD Phred normalized to 0-1 (Phred / 60, capped)
    - When multiple apply: take the maximum (most concerning wins)
    - No scores available: return None

    This replaces the unvalidated 0.4/0.6 weighted average. The score is
    for sorting and triage, not a clinical classification threshold.

    Parameters
    ----------
    consequence
        Functional consequence of the variant.
    revel_score
        Validated REVEL score (0.0-1.0) or None.
    spliceai_score
        Validated SpliceAI delta score (0.0-1.0) or None.
    cadd_phred
        Raw CADD Phred score or None.

    Returns
    -------
    Optional[float]
        Prioritization score (0.0-1.0), or None if no scores available.
    """
    scores: list[float] = []

    # REVEL: primary for missense (literature-validated)
    if consequence == FunctionalConsequence.MISSENSE and revel_score is not None:
        scores.append(revel_score)

    # SpliceAI: primary for splice-adjacent
    if consequence == FunctionalConsequence.SPLICE_SITE and spliceai_score is not None:
        scores.append(spliceai_score)

    # CADD: general deleteriousness (normalized to 0-1 scale)
    if cadd_phred is not None:
        scores.append(min(cadd_phred / 60.0, 1.0))

    # REVEL as supplementary even for non-missense if available
    if revel_score is not None and consequence != FunctionalConsequence.MISSENSE:
        scores.append(revel_score)

    # SpliceAI as supplementary for missense near splice sites
    if spliceai_score is not None and consequence != FunctionalConsequence.SPLICE_SITE:
        scores.append(spliceai_score)

    return max(scores) if scores else None
