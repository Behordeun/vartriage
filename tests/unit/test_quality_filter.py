"""Unit tests for QualityFilter."""

import warnings

import pytest

from vartriage.filter.quality_filter import QualityFilter
from vartriage.models.config import QualityFilterConfig
from vartriage.models.variant import Variant
from vartriage.models.warnings import MissingDataWarning


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    qual: float | None = 30.0,
    filter_status: str = "PASS",
) -> Variant:
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=qual,
        filter_status=filter_status,
    )


class TestQualityFilterBasic:
    """Core filtering logic."""

    def test_passes_variant_with_pass_filter_and_sufficient_qual(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variant = _make_variant(qual=25.0, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == [variant]

    def test_passes_variant_with_dot_filter(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variant = _make_variant(qual=25.0, filter_status=".")
        result = list(qf.apply(iter([variant])))
        assert result == [variant]

    def test_excludes_variant_with_failing_filter(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variant = _make_variant(qual=50.0, filter_status="LowQual")
        result = list(qf.apply(iter([variant])))
        assert result == []

    def test_excludes_variant_below_qual_threshold(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
        variant = _make_variant(qual=29.9, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == []

    def test_passes_variant_at_exact_threshold(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
        variant = _make_variant(qual=30.0, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == [variant]

    def test_uses_default_min_qual_of_20(self) -> None:
        qf = QualityFilter()
        below = _make_variant(qual=19.9, filter_status="PASS")
        at_threshold = _make_variant(qual=20.0, filter_status="PASS")
        result = list(qf.apply(iter([below, at_threshold])))
        assert result == [at_threshold]


class TestQualityFilterMissingQual:
    """Handling of missing QUAL scores."""

    def test_excludes_variant_with_none_qual(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variant = _make_variant(qual=None, filter_status="PASS")
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = list(qf.apply(iter([variant])))
        assert result == []

    def test_emits_warning_for_missing_qual(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variant = _make_variant(chrom="chr7", pos=54321, qual=None, filter_status="PASS")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            list(qf.apply(iter([variant])))

        assert len(caught) == 1
        warning_msg = caught[0].message
        assert isinstance(warning_msg, UserWarning)
        # The warning args contain our MissingDataWarning dataclass
        args = warning_msg.args
        assert len(args) == 1
        missing_data_warning = args[0]
        assert isinstance(missing_data_warning, MissingDataWarning)
        assert missing_data_warning.chrom == "chr7"
        assert missing_data_warning.pos == 54321

    def test_no_warning_when_filter_excludes_before_qual_check(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variant = _make_variant(qual=None, filter_status="LowQual")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            list(qf.apply(iter([variant])))
        assert len(caught) == 0


class TestQualityFilterOrdering:
    """Ordering preservation."""

    def test_preserves_input_order_of_passing_variants(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        variants = [
            _make_variant(chrom="chr1", pos=100, qual=50.0),
            _make_variant(chrom="chr2", pos=200, qual=10.0),  # excluded
            _make_variant(chrom="chr3", pos=300, qual=40.0),
            _make_variant(chrom="chr4", pos=400, qual=5.0),  # excluded
            _make_variant(chrom="chr5", pos=500, qual=60.0),
        ]
        result = list(qf.apply(iter(variants)))
        assert [v.chrom for v in result] == ["chr1", "chr3", "chr5"]


class TestQualityFilterEmptyInput:
    """Empty input handling."""

    def test_empty_input_produces_empty_output(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        result = list(qf.apply(iter([])))
        assert result == []

    def test_empty_input_no_warnings(self) -> None:
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            list(qf.apply(iter([])))
        assert len(caught) == 0
