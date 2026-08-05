"""Unit tests for ClinGen SV classifier."""

from __future__ import annotations

from vartriage.structural.classifier import SVClassifier
from vartriage.structural.models import (
    AnnotatedSV,
    GeneOverlap,
    ScoredSV,
    StructuralVariant,
    SVClassification,
    SVConsequence,
    SVEvidenceCategory,
    SVType,
)


def _make_scored_sv(
    chrom: str = "chr22",
    start: int = 18916842,
    end: int = 21465659,
    sv_type: SVType = SVType.DEL,
    consequence: SVConsequence = SVConsequence.WHOLE_GENE_DELETION,
    genes_affected: int = 1,
    hi_genes_affected: int = 1,
    population_frequency: float | None = None,
    frequency_unknown: bool = True,
    gene_overlaps: tuple[GeneOverlap, ...] = (),
    pathogenicity_score: float | None = 0.8,
) -> ScoredSV:
    sv = StructuralVariant(
        chrom=chrom,
        start=start,
        end=end,
        sv_type=sv_type,
    )
    annotated = AnnotatedSV(
        sv=sv,
        consequence=consequence,
        gene_overlaps=gene_overlaps,
        population_frequency=population_frequency,
        frequency_unknown=frequency_unknown,
        genes_affected=genes_affected,
        hi_genes_affected=hi_genes_affected,
    )
    return ScoredSV(
        annotated=annotated,
        pathogenicity_score=pathogenicity_score,
    )


class TestClassifierHIGeneDeletion:
    """Tests for pathogenic classification of HI gene deletions."""

    def test_hi_gene_whole_deletion_scores_high(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="TBX1",
            gene_chrom="chr22",
            gene_start=19744226,
            gene_end=19771115,
            overlap_fraction=1.0,
            is_whole_gene=True,
            exons_affected=9,
            total_exons=9,
            is_haploinsufficient=True,
            is_triplosensitive=False,
            hi_score=3.0,
            ts_score=None,
        )

        scored = _make_scored_sv(
            gene_overlaps=(overlap,),
            hi_genes_affected=1,
        )

        classifier = SVClassifier()
        result = list(classifier.classify(iter([scored])))[0]

        # HI gene fully deleted with Section 1B + Section 3A evidence
        # should get >= 0.90 (Likely Pathogenic or Pathogenic)
        assert result.classification in (
            SVClassification.PATHOGENIC,
            SVClassification.LIKELY_PATHOGENIC,
        )
        assert (
            SVEvidenceCategory.CONTAINS_ESTABLISHED_HI_GENE
            in result.evidence_categories
        )
        assert SVEvidenceCategory.GENE_FULLY_CONTAINED in result.evidence_categories

    def test_hi_gene_partial_deletion_is_likely_pathogenic_or_vus(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="BRCA1",
            gene_chrom="chr17",
            gene_start=43044295,
            gene_end=43170245,
            overlap_fraction=0.3,
            is_whole_gene=False,
            exons_affected=7,
            total_exons=23,
            is_haploinsufficient=True,
            is_triplosensitive=False,
            hi_score=3.0,
            ts_score=None,
        )

        scored = _make_scored_sv(
            chrom="chr17",
            start=43044295,
            end=43082000,
            consequence=SVConsequence.PARTIAL_GENE_DELETION,
            gene_overlaps=(overlap,),
        )

        classifier = SVClassifier()
        result = list(classifier.classify(iter([scored])))[0]

        assert result.classification in (
            SVClassification.LIKELY_PATHOGENIC,
            SVClassification.VUS,
        )
        assert SVEvidenceCategory.GENE_PARTIALLY_DELETED in result.evidence_categories


class TestClassifierCommonSVBenign:
    """Tests for benign classification of common SVs."""

    def test_common_sv_with_frequency_gets_benign_evidence(self) -> None:
        scored = _make_scored_sv(
            consequence=SVConsequence.INTERGENIC,
            genes_affected=0,
            hi_genes_affected=0,
            population_frequency=0.02,
            frequency_unknown=False,
            gene_overlaps=(),
        )

        classifier = SVClassifier()
        result = list(classifier.classify(iter([scored])))[0]

        # High AF pushes score negative
        assert result.evidence_score < 0.0

    def test_sv_fully_contained_in_benign_region(self) -> None:
        scored = _make_scored_sv(
            chrom="chr15",
            start=20200000,
            end=20400000,
            consequence=SVConsequence.INTERGENIC,
            genes_affected=0,
            hi_genes_affected=0,
            gene_overlaps=(),
        )

        benign_regions = [("chr15", 20143000, 20570000)]

        classifier = SVClassifier(benign_regions=benign_regions)
        result = list(classifier.classify(iter([scored])))[0]

        assert SVEvidenceCategory.CONTAINED_WITHIN_BENIGN in result.evidence_categories
        assert result.evidence_score < 0.0


class TestClassifierPathogenicRegionMatch:
    """Tests for pathogenic region overlap classification."""

    def test_complete_overlap_with_pathogenic_region(self) -> None:
        # 22q11.2 deletion syndrome coordinates
        scored = _make_scored_sv(
            chrom="chr22",
            start=18916842,
            end=21465659,
            consequence=SVConsequence.WHOLE_GENE_DELETION,
            gene_overlaps=(),
            genes_affected=0,
            hi_genes_affected=0,
        )

        pathogenic_regions = [("chr22", 18916842, 21465659)]
        region_names = {
            ("chr22", 18916842, 21465659): "22q11.2 deletion syndrome (DiGeorge)"
        }

        classifier = SVClassifier(
            pathogenic_regions=pathogenic_regions,
            pathogenic_region_names=region_names,
        )
        result = list(classifier.classify(iter([scored])))[0]

        assert (
            SVEvidenceCategory.COMPLETE_OVERLAP_PATHOGENIC in result.evidence_categories
        )
        assert result.syndrome_name == "22q11.2 deletion syndrome (DiGeorge)"
        assert result.evidence_score >= 0.45


class TestClassifierDuplicationSpecific:
    """Tests for duplication-specific classification (Section 4)."""

    def test_ts_gene_duplication_scored(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="PMP22",
            gene_chrom="chr17",
            gene_start=15138905,
            gene_end=15172379,
            overlap_fraction=1.0,
            is_whole_gene=True,
            exons_affected=5,
            total_exons=5,
            is_haploinsufficient=True,
            is_triplosensitive=True,
            hi_score=3.0,
            ts_score=3.0,
        )

        scored = _make_scored_sv(
            chrom="chr17",
            start=15100000,
            end=15200000,
            sv_type=SVType.DUP,
            consequence=SVConsequence.WHOLE_GENE_DUPLICATION,
            gene_overlaps=(overlap,),
            genes_affected=1,
            hi_genes_affected=1,
        )

        classifier = SVClassifier()
        result = list(classifier.classify(iter([scored])))[0]

        assert SVEvidenceCategory.DUP_TS_GENE_CONTAINED in result.evidence_categories


class TestClassifierMissingData:
    """Tests for missing data source tracking."""

    def test_missing_regions_tracked(self) -> None:
        scored = _make_scored_sv(
            gene_overlaps=(),
            genes_affected=0,
            hi_genes_affected=0,
        )

        # No pathogenic or benign regions provided
        classifier = SVClassifier()
        result = list(classifier.classify(iter([scored])))[0]

        assert "pathogenic_regions" in result.missing_data_sources
        assert "benign_regions" in result.missing_data_sources
