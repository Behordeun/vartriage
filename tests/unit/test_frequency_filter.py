"""Unit tests for FrequencyFilter."""

import pytest

from vartriage.models.config import PrioritizationConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    FunctionalConsequence,
    Variant,
)
from vartriage.prioritization.frequency_filter import FrequencyFilter


def _make_annotated(
    chrom: str = "chr1",
    pos: int = 100,
    allele_frequency: float | None = 0.005,
    frequency_unknown: bool = False,
) -> AnnotatedVariant:
    variant = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
    )
    return AnnotatedVariant(
        variant=variant,
        consequence=FunctionalConsequence.MISSENSE,
        allele_frequency=allele_frequency,
        frequency_unknown=frequency_unknown,
    )


class TestFrequencyFilterBasic:
    """Core allele frequency filtering logic."""

    def test_retains_variant_below_threshold(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(allele_frequency=0.005)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_retains_variant_at_exact_threshold(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(allele_frequency=0.01)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_excludes_variant_above_threshold(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(allele_frequency=0.05)
        result = list(ff.apply(iter([variant])))
        assert result == []

    def test_uses_default_threshold_of_0_01(self) -> None:
        ff = FrequencyFilter()
        below = _make_annotated(allele_frequency=0.009)
        above = _make_annotated(allele_frequency=0.02)
        result = list(ff.apply(iter([below, above])))
        assert result == [below]

    def test_retains_variant_with_zero_frequency(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(allele_frequency=0.0)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]


class TestFrequencyFilterUnknown:
    """Handling of frequency_unknown flag."""

    def test_retains_frequency_unknown_variant_regardless_of_af(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(
            allele_frequency=None, frequency_unknown=True
        )
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_retains_frequency_unknown_even_with_high_af(self) -> None:
        """Edge case: frequency_unknown=True should bypass even if AF is set high."""
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(
            allele_frequency=0.5, frequency_unknown=True
        )
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_retains_none_frequency_without_unknown_flag(self) -> None:
        """Variants with None AF but frequency_unknown=False still pass through."""
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variant = _make_annotated(allele_frequency=None, frequency_unknown=False)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]


class TestFrequencyFilterOrdering:
    """Ordering preservation."""

    def test_preserves_input_order_of_passing_variants(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variants = [
            _make_annotated(chrom="chr1", pos=100, allele_frequency=0.001),
            _make_annotated(chrom="chr2", pos=200, allele_frequency=0.05),
            _make_annotated(chrom="chr3", pos=300, allele_frequency=0.008),
            _make_annotated(chrom="chr4", pos=400, allele_frequency=0.9),
            _make_annotated(chrom="chr5", pos=500, allele_frequency=0.0001),
        ]
        result = list(ff.apply(iter(variants)))
        assert [v.variant.chrom for v in result] == ["chr1", "chr3", "chr5"]

    def test_mixed_known_and_unknown_preserves_order(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        variants = [
            _make_annotated(chrom="chr1", pos=1, allele_frequency=0.001),
            _make_annotated(chrom="chr2", pos=2, allele_frequency=0.5),
            _make_annotated(chrom="chr3", pos=3, frequency_unknown=True, allele_frequency=None),
            _make_annotated(chrom="chr4", pos=4, allele_frequency=0.002),
        ]
        result = list(ff.apply(iter(variants)))
        assert [v.variant.chrom for v in result] == ["chr1", "chr3", "chr4"]


class TestFrequencyFilterEmptyInput:
    """Empty input handling."""

    def test_empty_input_produces_empty_output(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.01))
        result = list(ff.apply(iter([])))
        assert result == []


class TestFrequencyFilterConfigValidation:
    """Configuration validation at construction time."""

    def test_rejects_negative_threshold(self) -> None:
        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            PrioritizationConfig(max_allele_frequency=-0.1)

    def test_rejects_threshold_above_one(self) -> None:
        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            PrioritizationConfig(max_allele_frequency=1.5)

    def test_accepts_threshold_at_zero(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.0))
        variant = _make_annotated(allele_frequency=0.0)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_accepts_threshold_at_one(self) -> None:
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=1.0))
        variant = _make_annotated(allele_frequency=1.0)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]
