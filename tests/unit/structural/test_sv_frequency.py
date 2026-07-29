"""Unit tests for reciprocal overlap frequency matching."""

from __future__ import annotations

import pytest

from vartriage.structural.annotator import SVAnnotator, SVFrequencyRecord
from vartriage.structural.models import SVType, StructuralVariant


def _make_sv(
    chrom: str = "chr1",
    start: int = 1000,
    end: int = 5000,
    sv_type: SVType = SVType.DEL,
) -> StructuralVariant:
    return StructuralVariant(
        chrom=chrom, start=start, end=end, sv_type=sv_type,
    )


class TestReciprocalOverlapMatching:
    """Tests for the reciprocal overlap frequency lookup algorithm."""

    def test_exact_match_returns_frequency(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 1000, 5000, "DEL", 0.005)
            ]
        }

        sv = _make_sv(start=1000, end=5000)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq == 0.005
        assert unknown is False

    def test_no_overlap_returns_unknown(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 50000, 60000, "DEL", 0.01)
            ]
        }

        sv = _make_sv(start=1000, end=5000)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq is None
        assert unknown is True

    def test_partial_overlap_below_threshold_returns_unknown(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        # Ref SV: 1000-10000 (9001 bp)
        # Query SV: 1000-3000 (2001 bp)
        # Overlap: 1000-3000 (2001 bp)
        # Frac of query: 2001/2001 = 1.0
        # Frac of ref: 2001/9001 = 0.22
        # Reciprocal = min(1.0, 0.22) = 0.22 < 0.5
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 1000, 10000, "DEL", 0.02)
            ]
        }

        sv = _make_sv(start=1000, end=3000)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq is None
        assert unknown is True

    def test_sufficient_reciprocal_overlap_matches(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        # Ref: 1000-6000 (5001 bp)
        # Query: 2000-7000 (5001 bp)
        # Overlap: 2000-6000 (4001 bp)
        # Frac of query: 4001/5001 = 0.80
        # Frac of ref: 4001/5001 = 0.80
        # Reciprocal = min(0.80, 0.80) = 0.80 >= 0.5
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 1000, 6000, "DEL", 0.003)
            ]
        }

        sv = _make_sv(start=2000, end=7000)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq == 0.003
        assert unknown is False

    def test_type_mismatch_does_not_match(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        # Same coordinates but different type
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 1000, 5000, "DUP", 0.01)
            ]
        }

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.DEL)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq is None
        assert unknown is True

    def test_empty_database_returns_unknown(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        annotator._sv_database = {}

        sv = _make_sv(start=1000, end=5000)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq is None
        assert unknown is True

    def test_best_match_selected_among_multiple(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 900, 5100, "DEL", 0.001),
                SVFrequencyRecord("chr1", 1000, 5000, "DEL", 0.008),
            ]
        }

        sv = _make_sv(start=1000, end=5000)
        freq, unknown = annotator._lookup_frequency(sv)

        # Exact match (second entry) has better reciprocal overlap
        assert freq == 0.008
        assert unknown is False

    def test_chr_prefix_normalization_in_frequency(self) -> None:
        annotator = SVAnnotator(reciprocal_overlap=0.5)
        annotator._genes = {}
        annotator._dosage = {}
        annotator._sv_database = {
            "chr1": [
                SVFrequencyRecord("chr1", 1000, 5000, "DEL", 0.005)
            ]
        }

        # Query SV without chr prefix should still match
        sv = _make_sv(chrom="1", start=1000, end=5000)
        freq, unknown = annotator._lookup_frequency(sv)

        assert freq == 0.005
        assert unknown is False
