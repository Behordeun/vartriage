"""Unit tests for structural variant data models."""

from __future__ import annotations

import pytest

from vartriage.structural.models import (
    AnnotatedSV,
    Breakpoint,
    ClassifiedSV,
    GeneOverlap,
    ScoredSV,
    SVClassification,
    SVConsequence,
    SV_CONSEQUENCE_SEVERITY,
    SVEvidenceCategory,
    SV_EVIDENCE_POINTS,
    SVType,
    StructuralVariant,
)


class TestSVType:
    def test_all_members_have_string_value(self) -> None:
        for member in SVType:
            assert isinstance(member.value, str)

    def test_del_value(self) -> None:
        assert SVType.DEL.value == "DEL"
        assert SVType.BND.value == "BND"


class TestSVConsequence:
    def test_severity_list_is_complete(self) -> None:
        assert set(SV_CONSEQUENCE_SEVERITY) == set(SVConsequence)

    def test_whole_gene_deletion_is_most_severe(self) -> None:
        assert SV_CONSEQUENCE_SEVERITY[0] == SVConsequence.WHOLE_GENE_DELETION

    def test_intergenic_is_least_severe(self) -> None:
        assert SV_CONSEQUENCE_SEVERITY[-1] == SVConsequence.INTERGENIC


class TestSVEvidenceCategory:
    def test_all_categories_have_point_values(self) -> None:
        for cat in SVEvidenceCategory:
            assert cat in SV_EVIDENCE_POINTS

    def test_pathogenic_overlap_gives_positive_points(self) -> None:
        assert SV_EVIDENCE_POINTS[SVEvidenceCategory.COMPLETE_OVERLAP_PATHOGENIC] > 0

    def test_benign_overlap_gives_negative_points(self) -> None:
        assert SV_EVIDENCE_POINTS[SVEvidenceCategory.CONTAINED_WITHIN_BENIGN] < 0


class TestBreakpoint:
    def test_creation_with_defaults(self) -> None:
        bp = Breakpoint(chrom="chr1", pos=12345)
        assert bp.chrom == "chr1"
        assert bp.pos == 12345
        assert bp.confidence_interval == (0, 0)

    def test_creation_with_ci(self) -> None:
        bp = Breakpoint(chrom="chrX", pos=100, confidence_interval=(-50, 50))
        assert bp.confidence_interval == (-50, 50)

    def test_frozen(self) -> None:
        bp = Breakpoint(chrom="chr1", pos=1)
        with pytest.raises(AttributeError):
            bp.pos = 2  # type: ignore[misc]


class TestGeneOverlap:
    def test_creation_with_dosage_flags(self) -> None:
        go = GeneOverlap(
            gene_symbol="BRCA1",
            gene_chrom="chr17",
            gene_start=1000,
            gene_end=5000,
            overlap_fraction=0.95,
            is_whole_gene=True,
            exons_affected=23,
            total_exons=23,
            is_haploinsufficient=True,
            is_triplosensitive=False,
            hi_score=3.0,
            ts_score=None,
        )
        assert go.is_haploinsufficient is True
        assert go.hi_score == 3.0

    def test_defaults_for_optional_fields(self) -> None:
        go = GeneOverlap(
            gene_symbol="X",
            gene_chrom="chr1",
            gene_start=1,
            gene_end=100,
            overlap_fraction=1.0,
            is_whole_gene=True,
            exons_affected=5,
            total_exons=5,
        )
        assert go.is_haploinsufficient is False
        assert go.hi_score is None


class TestAnnotatedSV:
    def test_defaults(self) -> None:
        sv = StructuralVariant(chrom="chr1", start=100, end=500, sv_type=SVType.DEL)
        annotated = AnnotatedSV(sv=sv, consequence=SVConsequence.INTERGENIC)
        assert annotated.gene_overlaps == ()
        assert annotated.population_frequency is None
        assert annotated.frequency_unknown is True
        assert annotated.genes_affected == 0
        assert annotated.hi_genes_affected == 0


class TestScoredSV:
    def test_defaults(self) -> None:
        sv = StructuralVariant(chrom="chr1", start=100, end=500, sv_type=SVType.DEL)
        annotated = AnnotatedSV(sv=sv, consequence=SVConsequence.INTERGENIC)
        scored = ScoredSV(annotated=annotated)
        assert scored.pathogenicity_score is None
        assert scored.dosage_score is None
        assert scored.size_score == 0.0
        assert scored.frequency_score == 1.0


class TestClassifiedSV:
    def test_defaults(self) -> None:
        sv = StructuralVariant(chrom="chr1", start=100, end=500, sv_type=SVType.DEL)
        annotated = AnnotatedSV(sv=sv, consequence=SVConsequence.INTERGENIC)
        scored = ScoredSV(annotated=annotated)
        classified = ClassifiedSV(scored=scored)
        assert classified.classification == SVClassification.VUS
        assert classified.evidence_categories == frozenset()
        assert classified.evidence_score == 0.0
        assert classified.missing_data_sources == frozenset()
        assert classified.syndrome_name is None
