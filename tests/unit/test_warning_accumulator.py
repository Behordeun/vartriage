"""Unit tests for the WarningAccumulator class."""

from __future__ import annotations

import warnings

import pytest

from vartriage._internal.warning_accumulator import (
    MissingDataSummaryWarning,
    WarningAccumulator,
    is_connection_failure,
)
from vartriage.models.config import MissingDataConfig
from vartriage.models.warnings import MissingDataWarning


class TestIsConnectionFailure:
    """Tests for the is_connection_failure helper function."""

    def test_none_reason_is_not_failure(self) -> None:
        assert is_connection_failure(None) is False

    def test_not_found_is_not_failure(self) -> None:
        assert is_connection_failure("not_found") is False

    def test_connection_error_is_failure(self) -> None:
        assert is_connection_failure("connection_error") is True

    def test_timeout_is_failure(self) -> None:
        assert is_connection_failure("timeout") is True

    def test_connection_timeout_is_failure(self) -> None:
        assert is_connection_failure("connection_timeout") is True

    def test_read_timeout_is_failure(self) -> None:
        assert is_connection_failure("read_timeout") is True

    def test_dns_error_is_failure(self) -> None:
        assert is_connection_failure("dns_error") is True

    def test_unknown_reason_is_not_failure(self) -> None:
        assert is_connection_failure("some_other_reason") is False


class TestWarningAccumulatorBasics:
    """Tests for basic accumulation behavior."""

    def test_initial_state(self) -> None:
        acc = WarningAccumulator()
        assert acc.total_count == 0
        assert acc.threshold == 1000
        assert acc.threshold_exceeded is False
        assert acc.sources == frozenset()
        assert acc.not_found_count == 0
        assert acc.connection_failure_count == 0
        assert acc.summary_emitted is False

    def test_custom_threshold(self) -> None:
        config = MissingDataConfig(warning_threshold=50)
        acc = WarningAccumulator(config)
        assert acc.threshold == 50

    def test_add_single_warning(self) -> None:
        acc = WarningAccumulator()
        warning = MissingDataWarning(
            chrom="chr1", pos=100, ref="A", alt="T",
            source="gnomAD", reason="not_found",
        )
        acc.add(warning)
        assert acc.total_count == 1
        assert acc.sources == frozenset({"gnomAD"})
        assert acc.not_found_count == 1
        assert acc.connection_failure_count == 0

    def test_add_multiple_sources(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found"))
        acc.add(MissingDataWarning("chr1", 200, "G", "C", "ClinVar", "not_found"))
        acc.add(MissingDataWarning("chr2", 300, "T", "A", "REVEL", "timeout"))
        assert acc.total_count == 3
        assert acc.sources == frozenset({"gnomAD", "ClinVar", "REVEL"})
        assert acc.not_found_count == 2
        assert acc.connection_failure_count == 1

    def test_count_by_source(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found"))
        acc.add(MissingDataWarning("chr1", 200, "G", "C", "gnomAD", "timeout"))
        acc.add(MissingDataWarning("chr2", 300, "T", "A", "ClinVar", "not_found"))
        assert acc.count_by_source == {"gnomAD": 2, "ClinVar": 1}

    def test_add_batch(self) -> None:
        acc = WarningAccumulator()
        batch = [
            MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found"),
            MissingDataWarning("chr1", 200, "G", "C", "ClinVar", "not_found"),
            MissingDataWarning("chr2", 300, "T", "A", "REVEL", "connection_error"),
        ]
        acc.add_batch(batch)
        assert acc.total_count == 3
        assert acc.connection_failure_count == 1

    def test_warnings_list_returns_copy(self) -> None:
        acc = WarningAccumulator()
        w = MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found")
        acc.add(w)
        result = acc.warnings_list
        assert result == [w]
        result.append(w)
        assert acc.total_count == 1


class TestWarningAccumulatorThreshold:
    """Tests for threshold behavior and summary emission."""

    def test_threshold_not_exceeded_at_boundary(self) -> None:
        config = MissingDataConfig(warning_threshold=3)
        acc = WarningAccumulator(config)
        for i in range(3):
            acc.add(MissingDataWarning(
                "chr1", i + 1, "A", "T", "gnomAD", "not_found"
            ))
        assert acc.total_count == 3
        assert acc.threshold_exceeded is False

    def test_threshold_exceeded_above_boundary(self) -> None:
        config = MissingDataConfig(warning_threshold=3)
        acc = WarningAccumulator(config)
        for i in range(4):
            acc.add(MissingDataWarning(
                "chr1", i + 1, "A", "T", "gnomAD", "not_found"
            ))
        assert acc.total_count == 4
        assert acc.threshold_exceeded is True

    def test_summary_emitted_on_threshold_exceeded(self) -> None:
        config = MissingDataConfig(warning_threshold=2)
        acc = WarningAccumulator(config)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            acc.add(MissingDataWarning("chr1", 1, "A", "T", "gnomAD", "not_found"))
            acc.add(MissingDataWarning("chr1", 2, "G", "C", "ClinVar", "not_found"))
            assert len(caught) == 0

            acc.add(MissingDataWarning("chr2", 3, "T", "A", "REVEL", "timeout"))
            assert len(caught) == 1
            assert issubclass(caught[0].category, MissingDataSummaryWarning)
            summary = caught[0].message
            assert summary.total_count == 3
            assert summary.sources == frozenset({"gnomAD", "ClinVar", "REVEL"})
            assert summary.not_found_count == 2
            assert summary.connection_failure_count == 1

    def test_summary_emitted_only_once(self) -> None:
        config = MissingDataConfig(warning_threshold=1)
        acc = WarningAccumulator(config)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            acc.add(MissingDataWarning("chr1", 1, "A", "T", "gnomAD", "not_found"))
            acc.add(MissingDataWarning("chr1", 2, "G", "C", "gnomAD", "not_found"))
            acc.add(MissingDataWarning("chr1", 3, "T", "A", "gnomAD", "not_found"))
            acc.add(MissingDataWarning("chr1", 4, "C", "G", "gnomAD", "not_found"))

            summary_warnings = [
                w for w in caught
                if issubclass(w.category, MissingDataSummaryWarning)
            ]
            assert len(summary_warnings) == 1

    def test_summary_emitted_flag(self) -> None:
        config = MissingDataConfig(warning_threshold=1)
        acc = WarningAccumulator(config)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            acc.add(MissingDataWarning("chr1", 1, "A", "T", "gnomAD", "not_found"))
            assert acc.summary_emitted is False

            acc.add(MissingDataWarning("chr1", 2, "G", "C", "gnomAD", "not_found"))
            assert acc.summary_emitted is True


class TestWarningAccumulatorReset:
    """Tests for the reset method."""

    def test_reset_clears_all_state(self) -> None:
        config = MissingDataConfig(warning_threshold=1)
        acc = WarningAccumulator(config)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            acc.add(MissingDataWarning("chr1", 1, "A", "T", "gnomAD", "not_found"))
            acc.add(MissingDataWarning("chr1", 2, "G", "C", "ClinVar", "timeout"))

        acc.reset()
        assert acc.total_count == 0
        assert acc.sources == frozenset()
        assert acc.not_found_count == 0
        assert acc.connection_failure_count == 0
        assert acc.summary_emitted is False
        assert acc.warnings_list == []


class TestMissingDataSummaryWarning:
    """Tests for the MissingDataSummaryWarning dataclass."""

    def test_str_representation(self) -> None:
        summary = MissingDataSummaryWarning(
            total_count=1500,
            sources=frozenset({"gnomAD", "ClinVar"}),
            not_found_count=1200,
            connection_failure_count=300,
        )
        text = str(summary)
        assert "1500" in text
        assert "gnomAD" in text
        assert "ClinVar" in text
        assert "1200" in text
        assert "300" in text

    def test_is_user_warning(self) -> None:
        summary = MissingDataSummaryWarning(total_count=5, sources=frozenset({"gnomAD"}))
        assert isinstance(summary, UserWarning)


class TestReasonDistinction:
    """Tests that the accumulator properly distinguishes not_found from connection failures."""

    def test_none_reason_counted_as_not_found(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", None))
        assert acc.not_found_count == 1
        assert acc.connection_failure_count == 0

    def test_not_found_reason_counted_correctly(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found"))
        assert acc.not_found_count == 1
        assert acc.connection_failure_count == 0

    def test_timeout_counted_as_connection_failure(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "timeout"))
        assert acc.not_found_count == 0
        assert acc.connection_failure_count == 1

    def test_connection_error_counted_as_connection_failure(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "connection_error"))
        assert acc.not_found_count == 0
        assert acc.connection_failure_count == 1

    def test_mixed_reasons(self) -> None:
        acc = WarningAccumulator()
        acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found"))
        acc.add(MissingDataWarning("chr1", 200, "G", "C", "gnomAD", "timeout"))
        acc.add(MissingDataWarning("chr1", 300, "T", "A", "ClinVar", None))
        acc.add(MissingDataWarning("chr1", 400, "C", "G", "REVEL", "connection_error"))
        assert acc.not_found_count == 2
        assert acc.connection_failure_count == 2
