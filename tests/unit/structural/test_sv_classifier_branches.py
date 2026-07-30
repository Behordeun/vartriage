"""Tests for classifier branches not covered by test_sv_classifier.py.

Focuses on Section 2 partial overlap scoring, Section 4 duplication
region matching, frequency evidence tiers, and region overlap helpers.
"""

from __future__ import annotations

import pytest

from vartriage.structural.classifier import SVClassifier
from vartriage.structural.models import (
    AnnotatedSV,
    GeneOverlap,
    ScoredSV,
    SVClassification,
    SVConsequence,
    SVEvidenceCategory,
    SVType,
    StructuralVariant,
)


def _make_scored(
    chrom: str = "chr1",
    start: int = 1000,
    end: int = 5000,
    sv_type: SVType = SVType.DEL,
    consequence: SVConsequence = SVConsequence.WHOLE_GENE_DELETION,
    genes_affected: int = 0,
    hi_genes_affected: int = 0,
    gene_overlaps: tuple[GeneOverlap, ...] = (),
    population_frequency: float | None = None,
    frequency_unknown: bool = True,
    copy_number: int | None = None,
) -> ScoredSV:
    sv = StructuralVariant(
        chrom=chrom, start=start, end=end, sv_type=sv_type,
        copy_number=copy_number,
    )
    annotated = AnnotatedSV(
        sv=sv, consequence=consequence,
        gene_overlaps=gene_overlaps,
        population_frequency=population_frequency,
        frequency_unknown=frequency_unknown,
        genes_affected=genes_affected,
        hi_genes_affected=hi_genes_affected,
    )
    return ScoredSV(annotated=annotated, pathogenicity_score=0.5)


class TestScoreToClassification:
    def test_vus_range(self) -> None:
        classifier = SVClassifier()
        assert classifier._score_to_classification(0.5) == SVClassification.VUS
        assert classifier._score_to_classification(0.0) == SVClassification.VUS
        assert classifier._score_to_classification(-0.5) == SVClassification.VUS

    def test_pathogenic_threshold(self) -> None:
        classifier = SVClassifier()
        assert classifier._score_to_classification(0.99) == SVClassification.PATHOGENIC
        assert classifier._score_to_classification(1.5) == SVClassification.PATHOGENIC

    def test_likely_pathogenic(self) -> None:
        classifier = SVClassifier()
        assert classifier._score_to_classification(0.90) == SVClassification.LIKELY_PATHOGENIC
        assert classifier._score_to_classification(0.95) == SVClassification.LIKELY_PATHOGENIC

    def test_likely_benign(self) -> None:
        classifier = SVClassifier()
        assert classifier._score_to_classification(-0.90) == SVClassification.LIKELY_BENIGN
        assert classifier._score_to_classification(-0.95) == SVClassification.LIKELY_BENIGN

    def test_benign_threshold(self) -> None:
        classifier = SVClassifier()
        assert classifier._score_to_classification(-0.99) == SVClassification.BENIGN
        assert classifier._score_to_classification(-2.0) == SVClassification.BENIGN


class TestIsLoss:
    def test_del_is_loss(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(sv_type=SVType.DEL)
        assert classifier._is_loss(scored.annotated) is True

    def test_dup_is_not_loss(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(sv_type=SVType.DUP)
        assert classifier._is_loss(scored.annotated) is False

    def test_cnv_with_low_copy_number_is_loss(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(sv_type=SVType.CNV, copy_number=1)
        assert classifier._is_loss(scored.annotated) is True

    def test_cnv_with_high_copy_number_is_not_loss(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(sv_type=SVType.CNV, copy_number=3)
        assert classifier._is_loss(scored.annotated) is False

    def test_cnv_without_copy_number_is_not_loss(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(sv_type=SVType.CNV, copy_number=None)
        assert classifier._is_loss(scored.annotated) is False


class TestFrequencyEvidence:
    def test_common_sv_gets_strong_benign_evidence(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(
            population_frequency=0.02, frequency_unknown=False,
        )
        result = list(classifier.classify(iter([scored])))[0]
        # AF >= 0.01 gives -0.60
        assert result.evidence_score < -0.5

    def test_moderately_common_gets_moderate_benign(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(
            population_frequency=0.007, frequency_unknown=False,
        )
        result = list(classifier.classify(iter([scored])))[0]
        assert result.evidence_score < -0.2

    def test_slightly_common_gets_weak_benign(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(
            population_frequency=0.002, frequency_unknown=False,
        )
        result = list(classifier.classify(iter([scored])))[0]
        assert result.evidence_score < 0.0

    def test_rare_sv_no_frequency_penalty(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(
            population_frequency=0.0005, frequency_unknown=False,
        )
        result = list(classifier.classify(iter([scored])))[0]
        # AF < 0.001 returns 0.0 from frequency evidence
        freq_contribution = classifier._evaluate_frequency_evidence(scored.annotated, set())
        assert freq_contribution == 0.0

    def test_unknown_frequency_no_penalty(self) -> None:
        classifier = SVClassifier()
        scored = _make_scored(frequency_unknown=True)
        freq_contribution = classifier._evaluate_frequency_evidence(scored.annotated, set())
        assert freq_contribution == 0.0


class TestSection3GeneEvaluation:
    def test_non_hi_gene_breakpoint_scores_low(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="GENEX", gene_chrom="chr1",
            gene_start=2000, gene_end=4000,
            overlap_fraction=0.3, is_whole_gene=False,
            exons_affected=3, total_exons=10,
            is_haploinsufficient=False, hi_score=None, ts_score=None,
        )
        classifier = SVClassifier()
        scored = _make_scored(
            sv_type=SVType.DEL,
            gene_overlaps=(overlap,),
            genes_affected=1,
        )
        result = list(classifier.classify(iter([scored])))[0]
        assert SVEvidenceCategory.BREAKPOINT_WITHIN_GENE in result.evidence_categories


class TestSection4Duplication:
    def test_partial_dup_gene_disruption(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="GENE1", gene_chrom="chr1",
            gene_start=2000, gene_end=4000,
            overlap_fraction=0.5, is_whole_gene=False,
            exons_affected=3, total_exons=10,
            is_haploinsufficient=False, is_triplosensitive=False,
            hi_score=None, ts_score=None,
        )
        classifier = SVClassifier()
        scored = _make_scored(
            sv_type=SVType.DUP,
            consequence=SVConsequence.PARTIAL_GENE_DUPLICATION,
            gene_overlaps=(overlap,),
            genes_affected=1,
        )
        result = list(classifier.classify(iter([scored])))[0]
        assert SVEvidenceCategory.DUP_GENE_DISRUPTED in result.evidence_categories

    def test_dup_identical_to_pathogenic_region(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="G", gene_chrom="chr1",
            gene_start=1000, gene_end=5000,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=5, total_exons=5,
            hi_score=None, ts_score=None,
        )
        pathogenic_regions = [("chr1", 1000, 5000)]
        classifier = SVClassifier(pathogenic_regions=pathogenic_regions)
        scored = _make_scored(
            sv_type=SVType.DUP,
            consequence=SVConsequence.WHOLE_GENE_DUPLICATION,
            gene_overlaps=(overlap,),
            genes_affected=1,
        )
        result = list(classifier.classify(iter([scored])))[0]
        assert SVEvidenceCategory.DUP_IDENTICAL_TO_PATHOGENIC in result.evidence_categories

    def test_dup_complete_overlap_pathogenic_region(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="G", gene_chrom="chr1",
            gene_start=1500, gene_end=4500,
            overlap_fraction=0.8, is_whole_gene=True,
            exons_affected=5, total_exons=5,
            hi_score=None, ts_score=None,
        )
        # SV from 1000-5000 covers region 1500-4500 completely
        pathogenic_regions = [("chr1", 1500, 4500)]
        classifier = SVClassifier(pathogenic_regions=pathogenic_regions)
        scored = _make_scored(
            sv_type=SVType.DUP,
            consequence=SVConsequence.WHOLE_GENE_DUPLICATION,
            gene_overlaps=(overlap,),
            genes_affected=1,
        )
        result = list(classifier.classify(iter([scored])))[0]
        assert SVEvidenceCategory.DUP_COMPLETE_OVERLAP_PATHOGENIC in result.evidence_categories


class TestRegionOverlapHelper:
    def test_no_overlap_returns_none(self) -> None:
        classifier = SVClassifier()
        result = classifier._best_region_overlap(
            "chr1", 1000, 2000, [("chr2", 1000, 2000)]
        )
        assert result is None

    def test_chr_prefix_normalization(self) -> None:
        classifier = SVClassifier()
        # SV on "1" should match region on "chr1"
        result = classifier._best_region_overlap(
            "1", 1000, 5000, [("chr1", 1000, 5000)]
        )
        assert result is not None
        assert result[0] >= 0.9  # query fraction
        assert result[1] >= 0.9  # ref fraction

    def test_partial_overlap_fractions(self) -> None:
        classifier = SVClassifier()
        # SV: 1000-3000 (2001 bp), Region: 2000-5000 (3001 bp)
        # Overlap: 2000-3000 (1001 bp)
        result = classifier._best_region_overlap(
            "chr1", 1000, 3000, [("chr1", 2000, 5000)]
        )
        assert result is not None
        frac_query, frac_ref = result
        assert 0.4 < frac_query < 0.6  # 1001/2001 ≈ 0.50
        assert 0.3 < frac_ref < 0.4    # 1001/3001 ≈ 0.33

    def test_best_region_overlap_with_key_returns_key(self) -> None:
        classifier = SVClassifier()
        regions = [("chr1", 500, 6000)]
        key, fracs = classifier._best_region_overlap_with_key(
            "chr1", 1000, 5000, regions
        )
        assert key == ("chr1", 500, 6000)
        assert fracs is not None
