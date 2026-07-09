"""Unit tests for pathogenicity score normalization and composite ranking."""

from __future__ import annotations

import warnings

import pytest

from vartriage.models.variant import (
    AnnotatedVariant,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.prioritization.scoring import (
    CADD_WEIGHT,
    REVEL_WEIGHT,
    ScoreValidationWarning,
    compute_composite_ranks,
    normalize_cadd_scores,
    score_variants,
    sort_by_composite_rank,
    validate_revel_scores,
)


def _make_annotated_variant(
    chrom: str = "chr1", pos: int = 100, ref: str = "A", alt: str = "T"
) -> AnnotatedVariant:
    """Helper to create a minimal AnnotatedVariant for testing."""
    v = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
    )
    return AnnotatedVariant(
        variant=v,
        consequence=FunctionalConsequence.MISSENSE,
    )


class TestNormalizeCaddScores:
    """Tests for CADD Phred score normalization."""

    def test_normalizes_score_to_fraction_of_99(self) -> None:
        result = normalize_cadd_scores([33.0])
        assert result == [pytest.approx(33.0 / 99.0)]

    def test_caps_at_one_for_scores_above_99(self) -> None:
        result = normalize_cadd_scores([100.0, 150.0, 99.0])
        assert result == [1.0, 1.0, 1.0]

    def test_zero_normalizes_to_zero(self) -> None:
        result = normalize_cadd_scores([0.0])
        assert result == [0.0]

    def test_none_values_remain_none(self) -> None:
        result = normalize_cadd_scores([None, 50.0, None])
        assert result[0] is None
        assert result[1] == pytest.approx(50.0 / 99.0)
        assert result[2] is None

    def test_negative_scores_rejected_with_warning(self) -> None:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = normalize_cadd_scores([-5.0, 30.0])
            assert result[0] is None
            assert result[1] == pytest.approx(30.0 / 99.0)
            assert len(w) == 1
            assert issubclass(w[0].category, ScoreValidationWarning)
            assert "-5.0" in str(w[0].message)

    def test_empty_list_returns_empty(self) -> None:
        assert normalize_cadd_scores([]) == []

    def test_all_none_returns_all_none(self) -> None:
        result = normalize_cadd_scores([None, None, None])
        assert result == [None, None, None]


class TestValidateRevelScores:
    """Tests for REVEL score validation."""

    def test_valid_scores_pass_through(self) -> None:
        result = validate_revel_scores([0.0, 0.5, 1.0])
        assert result == [0.0, 0.5, 1.0]

    def test_out_of_range_high_rejected(self) -> None:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_revel_scores([1.5])
            assert result == [None]
            assert len(w) == 1
            assert issubclass(w[0].category, ScoreValidationWarning)

    def test_negative_rejected(self) -> None:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_revel_scores([-0.1])
            assert result == [None]
            assert len(w) == 1
            assert issubclass(w[0].category, ScoreValidationWarning)

    def test_none_remains_none(self) -> None:
        result = validate_revel_scores([None, 0.8, None])
        assert result[0] is None
        assert result[1] == 0.8
        assert result[2] is None

    def test_empty_list_returns_empty(self) -> None:
        assert validate_revel_scores([]) == []


class TestComputeCompositeRanks:
    """Tests for composite rank computation."""

    def test_both_available_uses_weighted_formula(self) -> None:
        cadd = [0.5]
        revel = [0.8]
        result = compute_composite_ranks(cadd, revel)
        expected = 0.8 * REVEL_WEIGHT + 0.5 * CADD_WEIGHT
        assert result == [pytest.approx(expected)]

    def test_cadd_only_uses_cadd_directly(self) -> None:
        result = compute_composite_ranks([0.5], [None])
        assert result == [0.5]

    def test_revel_only_uses_revel_directly(self) -> None:
        result = compute_composite_ranks([None], [0.7])
        assert result == [0.7]

    def test_neither_available_returns_none(self) -> None:
        result = compute_composite_ranks([None], [None])
        assert result == [None]

    def test_mixed_batch(self) -> None:
        cadd = [0.5, None, 0.3, None]
        revel = [0.8, 0.6, None, None]
        result = compute_composite_ranks(cadd, revel)
        assert result[0] == pytest.approx(0.8 * 0.6 + 0.5 * 0.4)
        assert result[1] == pytest.approx(0.6)
        assert result[2] == pytest.approx(0.3)
        assert result[3] is None

    def test_mismatched_lengths_raises_valueerror(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            compute_composite_ranks([0.5, 0.3], [0.8])

    def test_empty_returns_empty(self) -> None:
        assert compute_composite_ranks([], []) == []


class TestSortByCompositeRank:
    """Tests for descending sort with nulls last."""

    def test_descending_order(self) -> None:
        v1 = _make_annotated_variant(pos=1)
        v2 = _make_annotated_variant(pos=2)
        v3 = _make_annotated_variant(pos=3)

        scored = [
            ScoredVariant(annotated=v1, composite_rank=0.3),
            ScoredVariant(annotated=v2, composite_rank=0.9),
            ScoredVariant(annotated=v3, composite_rank=0.6),
        ]
        result = sort_by_composite_rank(scored)
        ranks = [sv.composite_rank for sv in result]
        assert ranks == [0.9, 0.6, 0.3]

    def test_nulls_placed_last(self) -> None:
        v1 = _make_annotated_variant(pos=1)
        v2 = _make_annotated_variant(pos=2)
        v3 = _make_annotated_variant(pos=3)

        scored = [
            ScoredVariant(annotated=v1, composite_rank=None),
            ScoredVariant(annotated=v2, composite_rank=0.5),
            ScoredVariant(annotated=v3, composite_rank=None),
        ]
        result = sort_by_composite_rank(scored)
        assert result[0].composite_rank == 0.5
        assert result[1].composite_rank is None
        assert result[2].composite_rank is None

    def test_empty_list(self) -> None:
        assert sort_by_composite_rank([]) == []

    def test_all_nulls(self) -> None:
        v1 = _make_annotated_variant(pos=1)
        v2 = _make_annotated_variant(pos=2)
        scored = [
            ScoredVariant(annotated=v1, composite_rank=None),
            ScoredVariant(annotated=v2, composite_rank=None),
        ]
        result = sort_by_composite_rank(scored)
        assert len(result) == 2
        assert all(sv.composite_rank is None for sv in result)


class TestScoreVariants:
    """Integration tests for the full scoring pipeline."""

    def test_scores_and_sorts_batch(self) -> None:
        variants = [
            _make_annotated_variant(pos=1),
            _make_annotated_variant(pos=2),
            _make_annotated_variant(pos=3),
        ]
        cadd = [10.0, 50.0, 30.0]
        revel = [0.9, 0.3, 0.7]

        result = score_variants(variants, cadd, revel)

        assert len(result) == 3
        ranks = [sv.composite_rank for sv in result]
        for i in range(len(ranks) - 1):
            assert ranks[i] is not None
            assert ranks[i + 1] is not None
            assert ranks[i] >= ranks[i + 1]

    def test_emits_warning_when_no_scores(self) -> None:
        variants = [_make_annotated_variant(pos=1)]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = score_variants(variants, [None], [None])
            assert result[0].composite_rank is None
            warning_messages = [str(x.message) for x in w]
            assert any("MissingDataWarning" in msg for msg in warning_messages)

    def test_negative_cadd_excluded_from_composite(self) -> None:
        variants = [_make_annotated_variant(pos=1)]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = score_variants(variants, [-5.0], [0.8])
            assert result[0].cadd_normalized is None
            assert result[0].cadd_phred is None
            assert result[0].revel_score == 0.8
            assert result[0].composite_rank == pytest.approx(0.8)

    def test_out_of_range_revel_excluded(self) -> None:
        variants = [_make_annotated_variant(pos=1)]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = score_variants(variants, [30.0], [1.5])
            assert result[0].revel_score is None
            assert result[0].cadd_normalized == pytest.approx(30.0 / 99.0)
            assert result[0].composite_rank == pytest.approx(30.0 / 99.0)

    def test_mismatched_lengths_raises(self) -> None:
        variants = [_make_annotated_variant(pos=1)]
        with pytest.raises(ValueError, match="equal length"):
            score_variants(variants, [30.0, 40.0], [0.5])

    def test_empty_input(self) -> None:
        assert score_variants([], [], []) == []

    def test_preserves_annotated_variant_data(self) -> None:
        variant = _make_annotated_variant(chrom="chrX", pos=42)
        result = score_variants([variant], [20.0], [0.5])
        assert result[0].annotated.variant.chrom == "chrX"
        assert result[0].annotated.variant.pos == 42
