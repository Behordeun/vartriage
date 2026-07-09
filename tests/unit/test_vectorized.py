"""Unit tests for the vectorized helper utilities module.

Tests numpy-based score normalization helpers and polars-based batch join
helpers for coordinate overlaps and frequency lookups.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import pytest
from numpy.typing import NDArray

from vartriage._internal.vectorized import (
    _FALLBACK_CHUNK_SIZE,
    _MAX_MEMORY_BYTES,
    batch_coordinate_overlap_join,
    batch_frequency_join,
    batch_normalize_scores,
    compute_composite_vectorized,
    normalize_cadd_phred_vectorized,
    polars_available,
    validate_revel_vectorized,
)


class TestNormalizeCaddPhredVectorized:
    """Tests for normalize_cadd_phred_vectorized."""

    def test_zero_score_normalizes_to_zero(self) -> None:
        scores = np.array([0.0], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert result[0] == pytest.approx(0.0)

    def test_score_at_99_normalizes_to_one(self) -> None:
        scores = np.array([99.0], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert result[0] == pytest.approx(1.0)

    def test_score_above_99_caps_at_one(self) -> None:
        scores = np.array([150.0, 200.0], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert result[0] == pytest.approx(1.0)
        assert result[1] == pytest.approx(1.0)

    def test_intermediate_score_divides_by_99(self) -> None:
        scores = np.array([30.0, 60.0], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert result[0] == pytest.approx(30.0 / 99.0)
        assert result[1] == pytest.approx(60.0 / 99.0)

    def test_negative_scores_become_nan(self) -> None:
        scores = np.array([-1.0, -50.0], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert np.isnan(result[0])
        assert np.isnan(result[1])

    def test_nan_remains_nan(self) -> None:
        scores = np.array([np.nan, 30.0, np.nan], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert np.isnan(result[0])
        assert result[1] == pytest.approx(30.0 / 99.0)
        assert np.isnan(result[2])

    def test_empty_array(self) -> None:
        scores = np.array([], dtype=np.float64)
        result = normalize_cadd_phred_vectorized(scores)
        assert len(result) == 0

    def test_does_not_modify_input(self) -> None:
        scores = np.array([30.0, -5.0], dtype=np.float64)
        original = scores.copy()
        normalize_cadd_phred_vectorized(scores)
        np.testing.assert_array_equal(scores, original)


class TestValidateRevelVectorized:
    """Tests for validate_revel_vectorized."""

    def test_valid_range_preserved(self) -> None:
        scores = np.array([0.0, 0.5, 1.0], dtype=np.float64)
        result = validate_revel_vectorized(scores)
        assert result[0] == pytest.approx(0.0)
        assert result[1] == pytest.approx(0.5)
        assert result[2] == pytest.approx(1.0)

    def test_below_zero_becomes_nan(self) -> None:
        scores = np.array([-0.1, -1.0], dtype=np.float64)
        result = validate_revel_vectorized(scores)
        assert np.isnan(result[0])
        assert np.isnan(result[1])

    def test_above_one_becomes_nan(self) -> None:
        scores = np.array([1.01, 2.0], dtype=np.float64)
        result = validate_revel_vectorized(scores)
        assert np.isnan(result[0])
        assert np.isnan(result[1])

    def test_nan_remains_nan(self) -> None:
        scores = np.array([np.nan, 0.5], dtype=np.float64)
        result = validate_revel_vectorized(scores)
        assert np.isnan(result[0])
        assert result[1] == pytest.approx(0.5)

    def test_empty_array(self) -> None:
        scores = np.array([], dtype=np.float64)
        result = validate_revel_vectorized(scores)
        assert len(result) == 0


class TestComputeCompositeVectorized:
    """Tests for compute_composite_vectorized."""

    def test_both_available_uses_weighted_formula(self) -> None:
        cadd = np.array([0.5], dtype=np.float64)
        revel = np.array([0.8], dtype=np.float64)
        result = compute_composite_vectorized(cadd, revel)
        expected = 0.8 * 0.6 + 0.5 * 0.4
        assert result[0] == pytest.approx(expected)

    def test_cadd_only_uses_cadd(self) -> None:
        cadd = np.array([0.5], dtype=np.float64)
        revel = np.array([np.nan], dtype=np.float64)
        result = compute_composite_vectorized(cadd, revel)
        assert result[0] == pytest.approx(0.5)

    def test_revel_only_uses_revel(self) -> None:
        cadd = np.array([np.nan], dtype=np.float64)
        revel = np.array([0.7], dtype=np.float64)
        result = compute_composite_vectorized(cadd, revel)
        assert result[0] == pytest.approx(0.7)

    def test_neither_available_is_nan(self) -> None:
        cadd = np.array([np.nan], dtype=np.float64)
        revel = np.array([np.nan], dtype=np.float64)
        result = compute_composite_vectorized(cadd, revel)
        assert np.isnan(result[0])

    def test_custom_weights(self) -> None:
        cadd = np.array([0.5], dtype=np.float64)
        revel = np.array([0.8], dtype=np.float64)
        result = compute_composite_vectorized(cadd, revel, revel_weight=0.5, cadd_weight=0.5)
        expected = 0.8 * 0.5 + 0.5 * 0.5
        assert result[0] == pytest.approx(expected)

    def test_shape_mismatch_raises_value_error(self) -> None:
        cadd = np.array([0.5, 0.6], dtype=np.float64)
        revel = np.array([0.8], dtype=np.float64)
        with pytest.raises(ValueError, match="Array shapes must match"):
            compute_composite_vectorized(cadd, revel)

    def test_multiple_variants_mixed(self) -> None:
        cadd = np.array([0.5, 0.5, np.nan, np.nan], dtype=np.float64)
        revel = np.array([0.8, np.nan, 0.6, np.nan], dtype=np.float64)
        result = compute_composite_vectorized(cadd, revel)
        assert result[0] == pytest.approx(0.8 * 0.6 + 0.5 * 0.4)
        assert result[1] == pytest.approx(0.5)
        assert result[2] == pytest.approx(0.6)
        assert np.isnan(result[3])


class TestBatchNormalizeScores:
    """Tests for batch_normalize_scores."""

    def test_full_pipeline_both_scores(self) -> None:
        cadd_out, revel_out, comp_out = batch_normalize_scores(
            [30.0], [0.8]
        )
        assert cadd_out[0] == pytest.approx(30.0 / 99.0)
        assert revel_out[0] == pytest.approx(0.8)
        assert comp_out[0] == pytest.approx(0.8 * 0.6 + (30.0 / 99.0) * 0.4)

    def test_none_values_remain_none(self) -> None:
        cadd_out, revel_out, comp_out = batch_normalize_scores(
            [None], [None]
        )
        assert cadd_out[0] is None
        assert revel_out[0] is None
        assert comp_out[0] is None

    def test_negative_cadd_rejected(self) -> None:
        cadd_out, revel_out, comp_out = batch_normalize_scores(
            [-5.0], [0.5]
        )
        assert cadd_out[0] is None
        assert revel_out[0] == pytest.approx(0.5)
        assert comp_out[0] == pytest.approx(0.5)

    def test_invalid_revel_rejected(self) -> None:
        cadd_out, revel_out, comp_out = batch_normalize_scores(
            [30.0], [1.5]
        )
        assert cadd_out[0] == pytest.approx(30.0 / 99.0)
        assert revel_out[0] is None
        assert comp_out[0] == pytest.approx(30.0 / 99.0)

    def test_empty_inputs(self) -> None:
        cadd_out, revel_out, comp_out = batch_normalize_scores([], [])
        assert cadd_out == []
        assert revel_out == []
        assert comp_out == []

    def test_mismatched_lengths_raises(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            batch_normalize_scores([1.0], [1.0, 2.0])

    def test_large_batch_processes_correctly(self) -> None:
        n = 5000
        cadd_in: list[Optional[float]] = [float(i % 100) for i in range(n)]
        revel_in: list[Optional[float]] = [float(i % 100) / 100.0 for i in range(n)]
        cadd_out, revel_out, comp_out = batch_normalize_scores(cadd_in, revel_in)
        assert len(cadd_out) == n
        assert len(revel_out) == n
        assert len(comp_out) == n
        assert all(c is not None for c in comp_out)


class TestMemoryConstants:
    """Tests for memory budget constants."""

    def test_max_memory_is_2gb(self) -> None:
        assert _MAX_MEMORY_BYTES == 2 * 1024 * 1024 * 1024

    def test_fallback_chunk_size_is_500k(self) -> None:
        assert _FALLBACK_CHUNK_SIZE == 500_000


class TestPolarsAvailability:
    """Tests for polars availability check."""

    def test_polars_available_returns_bool(self) -> None:
        result = polars_available()
        assert isinstance(result, bool)


@pytest.mark.skipif(not polars_available(), reason="polars not installed")
class TestBatchFrequencyJoin:
    """Tests for polars-based batch_frequency_join."""

    def test_matching_variants_get_frequency(self) -> None:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("chrom\tpos\tref\talt\taf\n")
            f.write("chr1\t100\tA\tT\t0.05\n")
            f.write("chr1\t200\tG\tC\t0.001\n")
            ref_path = f.name

        variants = [("chr1", 100, "A", "T"), ("chr1", 200, "G", "C")]
        result = batch_frequency_join(variants, ref_path)
        assert result[0] == pytest.approx(0.05)
        assert result[1] == pytest.approx(0.001)
        Path(ref_path).unlink()

    def test_missing_variants_get_none(self) -> None:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("chrom\tpos\tref\talt\taf\n")
            f.write("chr1\t100\tA\tT\t0.05\n")
            ref_path = f.name

        variants = [("chr2", 500, "G", "A")]
        result = batch_frequency_join(variants, ref_path)
        assert result[0] is None
        Path(ref_path).unlink()

    def test_empty_variants_returns_empty(self) -> None:
        result = batch_frequency_join([], "/dev/null")
        assert result == []

    def test_unreadable_reference_raises_runtime_error(self) -> None:
        with pytest.raises(RuntimeError, match="Failed to read"):
            batch_frequency_join(
                [("chr1", 100, "A", "T")],
                "/nonexistent/file.tsv",
            )


@pytest.mark.skipif(not polars_available(), reason="polars not installed")
class TestBatchCoordinateOverlapJoin:
    """Tests for polars-based batch_coordinate_overlap_join."""

    def test_overlapping_regions_found(self) -> None:
        variants = [("chr1", 100, 101)]
        regions = [("chr1", 90, 150, "exon1"), ("chr1", 95, 110, "CDS1")]
        result = batch_coordinate_overlap_join(variants, regions)
        assert set(result[0]) == {"exon1", "CDS1"}

    def test_non_overlapping_variant_gets_empty_list(self) -> None:
        variants = [("chr2", 500, 501)]
        regions = [("chr1", 90, 150, "exon1")]
        result = batch_coordinate_overlap_join(variants, regions)
        assert result[0] == []

    def test_empty_variants_returns_empty(self) -> None:
        regions = [("chr1", 90, 150, "exon1")]
        result = batch_coordinate_overlap_join([], regions)
        assert result == []

    def test_empty_regions_returns_all_empty(self) -> None:
        variants = [("chr1", 100, 101), ("chr2", 200, 201)]
        result = batch_coordinate_overlap_join(variants, [])
        assert result == [[], []]

    def test_multiple_variants_mixed_overlaps(self) -> None:
        variants = [
            ("chr1", 100, 101),
            ("chr1", 250, 251),
            ("chr2", 50, 51),
        ]
        regions = [
            ("chr1", 90, 150, "exon1"),
            ("chr1", 200, 300, "exon2"),
        ]
        result = batch_coordinate_overlap_join(variants, regions)
        assert "exon1" in result[0]
        assert "exon2" in result[1]
        assert result[2] == []


@pytest.mark.skipif(polars_available(), reason="only runs without polars")
class TestPolarsNotAvailable:
    """Tests for polars not installed scenario."""

    def test_batch_frequency_join_raises_import_error(self) -> None:
        with pytest.raises(ImportError):
            batch_frequency_join([("chr1", 100, "A", "T")], "ref.tsv")

    def test_batch_coordinate_overlap_raises_import_error(self) -> None:
        with pytest.raises(ImportError):
            batch_coordinate_overlap_join(
                [("chr1", 100, 101)],
                [("chr1", 90, 150, "exon1")],
            )
