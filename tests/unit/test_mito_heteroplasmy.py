"""Unit tests for heteroplasmy extraction and classification."""

from __future__ import annotations

from vartriage.mito.heteroplasmy import (
    classify_level,
    extract_heteroplasmy,
)


class TestClassifyLevel:
    """Verify heteroplasmy threshold boundaries."""

    def test_homoplasmic_at_95_percent(self):
        assert classify_level(95.0) == "homoplasmic"

    def test_homoplasmic_at_100_percent(self):
        assert classify_level(100.0) == "homoplasmic"

    def test_high_at_60_percent(self):
        assert classify_level(60.0) == "high"

    def test_high_at_94_percent(self):
        assert classify_level(94.9) == "high"

    def test_moderate_at_20_percent(self):
        assert classify_level(20.0) == "moderate"

    def test_moderate_at_59_percent(self):
        assert classify_level(59.9) == "moderate"

    def test_low_at_1_percent(self):
        assert classify_level(1.0) == "low"

    def test_low_at_19_percent(self):
        assert classify_level(19.9) == "low"

    def test_sub_threshold_below_1_percent(self):
        assert classify_level(0.5) == "sub_threshold"

    def test_sub_threshold_at_zero(self):
        assert classify_level(0.0) == "sub_threshold"


class TestExtractHeteroplasmyFromAD:
    """Verify AD-field-based heteroplasmy extraction."""

    def test_high_heteroplasmy_from_ad(self):
        info = {"AD": [100, 900]}
        result = extract_heteroplasmy(info)
        assert result is not None
        assert abs(result.percentage - 90.0) < 0.1
        assert result.category == "high"
        assert result.depth == 1000

    def test_low_heteroplasmy_from_ad(self):
        info = {"AD": [950, 50]}
        result = extract_heteroplasmy(info)
        assert result is not None
        assert abs(result.percentage - 5.0) < 0.1
        assert result.category == "low"
        assert result.depth == 1000

    def test_moderate_heteroplasmy_from_ad(self):
        info = {"AD": [600, 400]}
        result = extract_heteroplasmy(info)
        assert result is not None
        assert abs(result.percentage - 40.0) < 0.1
        assert result.category == "moderate"

    def test_zero_depth_returns_none(self):
        info = {"AD": [0, 0]}
        result = extract_heteroplasmy(info)
        assert result is None

    def test_single_element_ad_returns_none(self):
        info = {"AD": [100]}
        result = extract_heteroplasmy(info)
        assert result is None


class TestExtractHeteroplasmyFromAF:
    """Verify AF-field fallback extraction."""

    def test_af_as_float(self):
        info = {"AF": 0.85}
        result = extract_heteroplasmy(info)
        assert result is not None
        assert abs(result.percentage - 85.0) < 0.1
        assert result.category == "high"
        assert result.depth == 0

    def test_af_as_list(self):
        info = {"AF": [0.45]}
        result = extract_heteroplasmy(info)
        assert result is not None
        assert abs(result.percentage - 45.0) < 0.1
        assert result.category == "moderate"

    def test_af_out_of_range_returns_none(self):
        info = {"AF": 1.5}
        result = extract_heteroplasmy(info)
        assert result is None

    def test_missing_ad_and_af_returns_none(self):
        info = {"DP": 100}
        result = extract_heteroplasmy(info)
        assert result is None

    def test_ad_takes_priority_over_af(self):
        info = {"AD": [200, 800], "AF": 0.5}
        result = extract_heteroplasmy(info)
        assert result is not None
        # AD gives 80%, AF would give 50% -- AD wins
        assert abs(result.percentage - 80.0) < 0.1

    def test_negative_ad_values_returns_none(self):
        info = {"AD": [-1, 2]}
        result = extract_heteroplasmy(info)
        assert result is None
