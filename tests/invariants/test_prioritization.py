"""Hypothesis tests for prioritization engine.

Covers allele frequency filtering, composite pathogenicity score computation,
and rank ordering of scored variants.
"""

from __future__ import annotations

import warnings
from typing import Optional

from hypothesis import assume, given, settings
from hypothesis import strategies as st

from tests.generators.variants import chromosome, genomic_position, snv_allele
from vartriage.models.config import PrioritizationConfig
from vartriage.models.variant import (AnnotatedVariant, FunctionalConsequence,
                                      ScoredVariant, Variant)
from vartriage.prioritization.frequency_filter import FrequencyFilter
from vartriage.prioritization.scoring import (compute_composite_ranks,
                                              normalize_cadd_scores,
                                              score_variants,
                                              sort_by_composite_rank,
                                              validate_revel_scores)

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

VALID_MAX_AF = st.floats(
    min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False
)

VALID_ALLELE_FREQUENCY = st.floats(
    min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False
)

VALID_CADD_PHRED = st.floats(
    min_value=0.0, max_value=150.0, allow_nan=False, allow_infinity=False
)

VALID_REVEL = st.floats(
    min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False
)


@st.composite
def annotated_variant_for_frequency_filter(
    draw: st.DrawFn,
    af: Optional[float] = None,
    frequency_unknown: Optional[bool] = None,
) -> AnnotatedVariant:
    """Generate an AnnotatedVariant with configurable AF and frequency_unknown."""
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())
    consequence = draw(st.sampled_from(list(FunctionalConsequence)))

    if af is None:
        allele_frequency: Optional[float] = draw(
            st.one_of(st.none(), VALID_ALLELE_FREQUENCY)
        )
    else:
        allele_frequency = af

    if frequency_unknown is None:
        freq_unknown = draw(st.booleans())
    else:
        freq_unknown = frequency_unknown

    variant = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
        info={},
    )

    return AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=allele_frequency,
        frequency_unknown=freq_unknown,
    )


@st.composite
def annotated_variant_with_known_af(
    draw: st.DrawFn,
) -> AnnotatedVariant:
    """Generate an AnnotatedVariant with a known (non-None) AF and frequency_unknown=False."""
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())
    consequence = draw(st.sampled_from(list(FunctionalConsequence)))
    allele_frequency = draw(VALID_ALLELE_FREQUENCY)

    variant = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
        info={},
    )

    return AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=allele_frequency,
        frequency_unknown=False,
    )


# ---------------------------------------------------------------------------


@given(
    variant=annotated_variant_with_known_af(),
    max_af=VALID_MAX_AF,
)
@settings(max_examples=200)
def test_af_above_threshold_excluded(variant: AnnotatedVariant, max_af: float) -> None:
    """Variants with allele_frequency > max_af are excluded."""
    assume(variant.allele_frequency is not None)
    assume(variant.allele_frequency > max_af)

    config = PrioritizationConfig(max_allele_frequency=max_af)
    ff = FrequencyFilter(config)
    result = list(ff.apply(iter([variant])))

    assert result == [], (
        f"Variant with AF={variant.allele_frequency} should be excluded "
        f"when max_af={max_af}"
    )


@given(
    variant=annotated_variant_with_known_af(),
    max_af=VALID_MAX_AF,
)
@settings(max_examples=200)
def test_af_at_or_below_threshold_retained(
    variant: AnnotatedVariant, max_af: float
) -> None:
    """Variants with allele_frequency <= max_af are retained."""
    assume(variant.allele_frequency is not None)
    assume(variant.allele_frequency <= max_af)

    config = PrioritizationConfig(max_allele_frequency=max_af)
    ff = FrequencyFilter(config)
    result = list(ff.apply(iter([variant])))

    assert result == [variant], (
        f"Variant with AF={variant.allele_frequency} should be retained "
        f"when max_af={max_af}"
    )


@given(
    max_af=VALID_MAX_AF,
    data=st.data(),
)
@settings(max_examples=200)
def test_frequency_unknown_always_retained(max_af: float, data: st.DataObject) -> None:
    """Variants with frequency_unknown=True are always retained regardless of threshold."""
    variant = data.draw(annotated_variant_for_frequency_filter(frequency_unknown=True))

    config = PrioritizationConfig(max_allele_frequency=max_af)
    ff = FrequencyFilter(config)
    result = list(ff.apply(iter([variant])))

    assert result == [variant], (
        f"Variant with frequency_unknown=True should always be retained "
        f"regardless of AF={variant.allele_frequency} and max_af={max_af}"
    )


# ---------------------------------------------------------------------------


@given(
    cadd_phred=VALID_CADD_PHRED,
    revel_score=VALID_REVEL,
)
@settings(max_examples=200)
def test_composite_both_scores_available(cadd_phred: float, revel_score: float) -> None:
    """When both CADD and REVEL are available, composite = (revel * 0.6) + (cadd_norm * 0.4)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cadd_normalized = normalize_cadd_scores([cadd_phred])
        revel_validated = validate_revel_scores([revel_score])

    assert cadd_normalized[0] is not None
    assert revel_validated[0] is not None

    composites = compute_composite_ranks(cadd_normalized, revel_validated)
    result = composites[0]

    expected_cadd_norm = min(cadd_phred / 99.0, 1.0)
    expected = (revel_score * 0.6) + (expected_cadd_norm * 0.4)

    assert result is not None
    assert abs(result - expected) < 1e-10, (
        f"Composite {result} != expected {expected} for "
        f"CADD={cadd_phred}, REVEL={revel_score}"
    )


@given(revel_score=VALID_REVEL)
@settings(max_examples=200)
def test_composite_revel_only(revel_score: float) -> None:
    """When only REVEL is available, composite equals the REVEL score."""
    cadd_normalized = [None]
    revel_validated = validate_revel_scores([revel_score])

    assert revel_validated[0] is not None

    composites = compute_composite_ranks(cadd_normalized, revel_validated)
    result = composites[0]

    assert result is not None
    assert (
        abs(result - revel_score) < 1e-10
    ), f"Composite {result} != REVEL {revel_score} when CADD is missing"


@given(cadd_phred=VALID_CADD_PHRED)
@settings(max_examples=200)
def test_composite_cadd_only(cadd_phred: float) -> None:
    """When only CADD is available, composite equals min(cadd_phred / 99, 1.0)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cadd_normalized = normalize_cadd_scores([cadd_phred])

    assert cadd_normalized[0] is not None

    revel_validated: list[Optional[float]] = [None]
    composites = compute_composite_ranks(cadd_normalized, revel_validated)
    result = composites[0]

    expected = min(cadd_phred / 99.0, 1.0)

    assert result is not None
    assert (
        abs(result - expected) < 1e-10
    ), f"Composite {result} != expected CADD_norm {expected} when REVEL is missing"


@settings(max_examples=100)
@given(data=st.data())
def test_composite_neither_score_is_null(data: st.DataObject) -> None:
    """When neither CADD nor REVEL is available, composite is null."""
    cadd_normalized: list[Optional[float]] = [None]
    revel_validated: list[Optional[float]] = [None]

    composites = compute_composite_ranks(cadd_normalized, revel_validated)
    result = composites[0]

    assert (
        result is None
    ), f"Composite should be None when neither score is available, got {result}"


@given(
    cadd_phred=VALID_CADD_PHRED,
)
@settings(max_examples=200)
def test_cadd_normalization_formula(cadd_phred: float) -> None:
    """CADD normalization applies min(score/99.0, 1.0)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = normalize_cadd_scores([cadd_phred])

    expected = min(cadd_phred / 99.0, 1.0)
    assert result[0] is not None
    assert (
        abs(result[0] - expected) < 1e-10
    ), f"Normalized CADD {result[0]} != expected {expected} for input {cadd_phred}"


# ---------------------------------------------------------------------------


@st.composite
def scored_variant_with_rank(
    draw: st.DrawFn,
    rank: Optional[float] = None,
) -> ScoredVariant:
    """Generate a ScoredVariant with a specified or random composite_rank."""
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())
    consequence = draw(st.sampled_from(list(FunctionalConsequence)))

    variant = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
        info={},
    )

    annotated = AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=0.005,
        frequency_unknown=False,
    )

    if rank is None:
        composite_rank: Optional[float] = draw(
            st.one_of(
                st.none(),
                st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
            )
        )
    else:
        composite_rank = rank

    return ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
    )


@given(
    variants=st.lists(scored_variant_with_rank(), min_size=0, max_size=50),
)
@settings(max_examples=200)
def test_rank_ordering_descending_nulls_last(
    variants: list[ScoredVariant],
) -> None:
    """Output is strictly descending by composite_rank with nulls at the end."""
    result = sort_by_composite_rank(variants)

    # Split into ranked and null-ranked
    ranked = [v for v in result if v.composite_rank is not None]
    null_ranked = [v for v in result if v.composite_rank is None]

    # All ranked variants come before null-ranked
    expected_length = len(ranked) + len(null_ranked)
    assert len(result) == expected_length == len(variants)

    if ranked and null_ranked:
        ranked_end_idx = len(ranked) - 1
        assert result[ranked_end_idx].composite_rank is not None
        assert result[ranked_end_idx + 1].composite_rank is None

    # Ranked variants are in descending order
    for i in range(len(ranked) - 1):
        assert ranked[i].composite_rank is not None
        assert ranked[i + 1].composite_rank is not None
        assert ranked[i].composite_rank >= ranked[i + 1].composite_rank, (
            f"Rank ordering violated: {ranked[i].composite_rank} < "
            f"{ranked[i + 1].composite_rank} at positions {i}, {i+1}"
        )


@given(
    data=st.data(),
)
@settings(max_examples=200)
def test_score_variants_output_ordering(data: st.DataObject) -> None:
    """score_variants produces output sorted descending by composite_rank, nulls last."""
    n = data.draw(st.integers(min_value=1, max_value=20))

    variants: list[AnnotatedVariant] = []
    cadd_scores: list[Optional[float]] = []
    revel_scores: list[Optional[float]] = []

    for _ in range(n):
        v = data.draw(annotated_variant_for_frequency_filter(frequency_unknown=False))
        variants.append(v)
        cadd_scores.append(data.draw(st.one_of(st.none(), VALID_CADD_PHRED)))
        revel_scores.append(data.draw(st.one_of(st.none(), VALID_REVEL)))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = score_variants(variants, cadd_scores, revel_scores)

    # Verify descending order with nulls last
    ranked = [sv for sv in result if sv.composite_rank is not None]
    null_ranked = [sv for sv in result if sv.composite_rank is None]

    # All ranked come before nulls
    for i, sv in enumerate(result):
        if sv.composite_rank is None:
            # Everything after this should also be None
            for j in range(i, len(result)):
                assert (
                    result[j].composite_rank is None
                ), f"Found non-null rank after null at position {j}"
            break

    # Ranked variants are in descending order
    for i in range(len(ranked) - 1):
        assert ranked[i].composite_rank is not None
        assert ranked[i + 1].composite_rank is not None
        assert ranked[i].composite_rank >= ranked[i + 1].composite_rank, (
            f"score_variants output ordering violated: "
            f"{ranked[i].composite_rank} < {ranked[i + 1].composite_rank}"
        )
