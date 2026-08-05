"""Unit tests for SV gene overlap annotation."""

from __future__ import annotations

from vartriage.structural.annotator import GeneRecord, SVAnnotator
from vartriage.structural.models import (
    StructuralVariant,
    SVConsequence,
    SVType,
)


def _make_sv(
    chrom: str = "chr1",
    start: int = 1000,
    end: int = 5000,
    sv_type: SVType = SVType.DEL,
    copy_number: int | None = None,
) -> StructuralVariant:
    return StructuralVariant(
        chrom=chrom,
        start=start,
        end=end,
        sv_type=sv_type,
        copy_number=copy_number,
    )


class TestGeneOverlapClassification:
    """Tests for classifying SV-gene overlap type."""

    def test_whole_gene_deletion(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        # Gene is 1000-4000, SV is 500-5000 (covers 100% of gene)
        annotator._genes = {"chr1": [GeneRecord("BRCA1", "chr1", 1000, 4000, "+", 23)]}
        annotator._dosage = {}

        sv = _make_sv(start=500, end=5000, sv_type=SVType.DEL)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.WHOLE_GENE_DELETION
        assert result.genes_affected == 1
        assert result.gene_overlaps[0].is_whole_gene is True
        assert result.gene_overlaps[0].overlap_fraction >= 0.8

    def test_partial_gene_deletion(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        # Gene is 1000-10000, SV is 1000-3000 (covers ~22% of gene)
        annotator._genes = {"chr1": [GeneRecord("TP53", "chr1", 1000, 10000, "+", 11)]}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=3000, sv_type=SVType.DEL)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.PARTIAL_GENE_DELETION
        assert result.gene_overlaps[0].is_whole_gene is False

    def test_whole_gene_duplication(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("MYC", "chr1", 2000, 3000, "+", 3)]}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.DUP)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.WHOLE_GENE_DUPLICATION

    def test_gene_disruption_for_inversion(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("NF1", "chr1", 2000, 4000, "+", 58)]}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.INV)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.GENE_DISRUPTION

    def test_intergenic_when_no_gene_overlap(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("FAR", "chr1", 50000, 60000, "+", 5)]}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.DEL)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.INTERGENIC
        assert result.genes_affected == 0

    def test_intergenic_when_no_genes_loaded(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.DEL)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.INTERGENIC

    def test_cnv_loss_classified_as_deletion(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("RB1", "chr1", 2000, 3000, "+", 27)]}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.CNV, copy_number=1)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.WHOLE_GENE_DELETION

    def test_cnv_gain_classified_as_duplication(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("ERBB2", "chr1", 2000, 3000, "+", 31)]}
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.CNV, copy_number=4)
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.WHOLE_GENE_DUPLICATION


class TestMultiGeneOverlap:
    """Tests for SVs affecting multiple genes."""

    def test_multiple_genes_counted(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {
            "chr1": [
                GeneRecord("GENE1", "chr1", 1500, 2500, "+", 5),
                GeneRecord("GENE2", "chr1", 3000, 4000, "+", 8),
                GeneRecord("GENE3", "chr1", 4500, 4800, "+", 3),
            ]
        }
        annotator._dosage = {}

        sv = _make_sv(start=1000, end=5000, sv_type=SVType.DEL)
        result = annotator._annotate_single(sv)

        assert result.genes_affected == 3

    def test_overlaps_sorted_by_fraction_descending(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {
            "chr1": [
                GeneRecord("SMALL", "chr1", 2000, 2100, "+", 2),
                GeneRecord("LARGE", "chr1", 1000, 50000, "+", 40),
            ]
        }
        annotator._dosage = {}

        sv = _make_sv(start=1500, end=3000, sv_type=SVType.DEL)
        result = annotator._annotate_single(sv)

        # SMALL gene (100bp) is fully covered by the SV
        # LARGE gene (49000bp) is barely touched
        assert result.gene_overlaps[0].gene_symbol == "SMALL"


class TestChromNormalization:
    """Tests for chr prefix normalization in gene lookup."""

    def test_sv_with_chr_matches_genes_with_chr(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("X", "chr1", 2000, 3000, "+", 5)]}
        annotator._dosage = {}

        sv = _make_sv(chrom="chr1", start=1000, end=5000)
        result = annotator._annotate_single(sv)
        assert result.genes_affected == 1

    def test_sv_without_chr_matches_genes_with_chr(self) -> None:
        annotator = SVAnnotator(whole_gene_threshold=0.8)
        annotator._genes = {"chr1": [GeneRecord("X", "chr1", 2000, 3000, "+", 5)]}
        annotator._dosage = {}

        sv = _make_sv(chrom="1", start=1000, end=5000)
        result = annotator._annotate_single(sv)
        assert result.genes_affected == 1
