"""Unit tests for SV pathogenicity scoring engine."""

from __future__ import annotations

import pytest

from vartriage.structural.models import (
    AnnotatedSV,
    GeneOverlap,
    ScoredSV,
    SVConsequence,
    SVType,
    StructuralVariant,
)
from vartriage.structural.scoring import SVScorer


def _make_annotated(
    sv_type: SVType = SVType.DEL,
    consequence: SVConsequence = SVConsequence.WHOLE_GENE_DELETION,
    genes_affected: int = 1,
    hi_genes_affected: int = 0,
    population_frequency: float | None = None,
    frequency_unknown: bool = True,
    gene_overlaps: tuple[GeneOverlap, ...] = (),
    length: int = 50_000,
    copy_number: int | None = None,
) -> AnnotatedSV:
    sv = StructuralVariant(
        chrom="chr1", start=1000, end=1000 + length - 1,
        sv_type=sv_type, copy_number=copy_number,
    )
    return AnnotatedSV(
        sv=sv,
        consequence=consequence,
        gene_overlaps=gene_overlaps,
        population_frequency=population_frequency,
        frequency_unknown=frequency_unknown,
        genes_affected=genes_affected,
        hi_genes_affected=hi_genes_affected,
    )


class TestFrequencyFilter:
    def test_common_sv_filtered_out(self) -> None:
        scorer = SVScorer(max_allele_frequency=0.01)
        annotated = _make_annotated(population_frequency=0.05, frequency_unknown=False)
        results = list(scorer.score(iter([annotated])))
        assert len(results) == 0

    def test_rare_sv_passes(self) -> None:
        scorer = SVScorer(max_allele_frequency=0.01)
        annotated = _make_annotated(population_frequency=0.005, frequency_unknown=False)
        results = list(scorer.score(iter([annotated])))
        assert len(results) == 1

    def test_unknown_frequency_passes(self) -> None:
        scorer = SVScorer(max_allele_frequency=0.01)
        annotated = _make_annotated(population_frequency=None, frequency_unknown=True)
        results = list(scorer.score(iter([annotated])))
        assert len(results) == 1


class TestImpactScore:
    def test_whole_gene_deletion_gets_highest_base(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(consequence=SVConsequence.WHOLE_GENE_DELETION)
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.pathogenicity_score is not None
        assert scored.pathogenicity_score > 0.3

    def test_intergenic_with_known_frequency_gets_low_score(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(
            consequence=SVConsequence.INTERGENIC,
            genes_affected=0,
            population_frequency=0.001,
            frequency_unknown=False,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.pathogenicity_score is not None
        assert scored.pathogenicity_score < 0.3

    def test_intergenic_frequency_unknown_returns_none_score(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(
            consequence=SVConsequence.INTERGENIC,
            genes_affected=0,
            frequency_unknown=True,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.pathogenicity_score is None

    def test_multi_gene_boost(self) -> None:
        scorer = SVScorer()
        single = _make_annotated(
            consequence=SVConsequence.PARTIAL_GENE_DELETION, genes_affected=1,
        )
        multi = _make_annotated(
            consequence=SVConsequence.PARTIAL_GENE_DELETION, genes_affected=5,
        )
        s1 = list(scorer.score(iter([single])))[0]
        s5 = list(scorer.score(iter([multi])))[0]
        assert s5.pathogenicity_score > s1.pathogenicity_score


class TestDosageScore:
    def test_hi_gene_loss_scores_high_dosage(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="BRCA1", gene_chrom="chr17",
            gene_start=1000, gene_end=5000,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=23, total_exons=23,
            is_haploinsufficient=True, hi_score=3.0,
            ts_score=None,
        )
        scorer = SVScorer()
        annotated = _make_annotated(
            sv_type=SVType.DEL,
            gene_overlaps=(overlap,),
            genes_affected=1,
            hi_genes_affected=1,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.dosage_score is not None
        assert scored.dosage_score >= 0.9

    def test_ts_gene_gain_scores_dosage(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="PMP22", gene_chrom="chr17",
            gene_start=1000, gene_end=5000,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=5, total_exons=5,
            is_haploinsufficient=False, hi_score=None,
            is_triplosensitive=True, ts_score=3.0,
        )
        scorer = SVScorer()
        annotated = _make_annotated(
            sv_type=SVType.DUP,
            consequence=SVConsequence.WHOLE_GENE_DUPLICATION,
            gene_overlaps=(overlap,),
            genes_affected=1,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.dosage_score >= 0.9

    def test_cnv_loss_uses_hi_score(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="RB1", gene_chrom="chr13",
            gene_start=1000, gene_end=5000,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=27, total_exons=27,
            is_haploinsufficient=True, hi_score=3.0,
            ts_score=None,
        )
        scorer = SVScorer()
        annotated = _make_annotated(
            sv_type=SVType.CNV,
            copy_number=1,
            gene_overlaps=(overlap,),
            genes_affected=1,
            hi_genes_affected=1,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.dosage_score >= 0.9

    def test_no_dosage_data_gives_modest_default(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="GENE1", gene_chrom="chr1",
            gene_start=1000, gene_end=5000,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=5, total_exons=5,
            hi_score=None, ts_score=None,
        )
        scorer = SVScorer()
        annotated = _make_annotated(
            gene_overlaps=(overlap,), genes_affected=1,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        # Default is 0.3 when no dosage data but genes are affected
        assert scored.dosage_score == pytest.approx(0.3)


class TestSizeScore:
    def test_large_sv_scores_high(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(length=2_000_000)
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.size_score == 1.0

    def test_medium_sv_scores_mid(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(length=500_000)
        scored = list(scorer.score(iter([annotated])))[0]
        assert 0.5 < scored.size_score < 1.0

    def test_small_sv_scores_low(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(length=5_000)
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.size_score < 0.3

    def test_tiny_sv_scores_near_zero(self) -> None:
        scorer = SVScorer()
        annotated = _make_annotated(length=100)
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.size_score < 0.05


class TestFrequencyScore:
    def test_absent_from_population_gets_max(self) -> None:
        scorer = SVScorer(max_allele_frequency=0.01)
        annotated = _make_annotated(frequency_unknown=True)
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.frequency_score == 1.0

    def test_at_threshold_boundary_gets_zero(self) -> None:
        scorer = SVScorer(max_allele_frequency=0.01)
        annotated = _make_annotated(
            population_frequency=0.01, frequency_unknown=False,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.frequency_score == pytest.approx(0.0)

    def test_half_threshold_gives_half(self) -> None:
        scorer = SVScorer(max_allele_frequency=0.01)
        annotated = _make_annotated(
            population_frequency=0.005, frequency_unknown=False,
        )
        scored = list(scorer.score(iter([annotated])))[0]
        assert scored.frequency_score == pytest.approx(0.5)
