"""Unit tests for SV clinical report builder."""

from __future__ import annotations

import pytest

from vartriage.structural.models import (
    AnnotatedSV,
    ClassifiedSV,
    GeneOverlap,
    ScoredSV,
    SVClassification,
    SVConsequence,
    SVEvidenceCategory,
    SVType,
    StructuralVariant,
)
from vartriage.structural.report import SVReportBuilder, _format_size


def _make_classified(
    chrom: str = "chr22",
    start: int = 18916842,
    end: int = 21465659,
    sv_type: SVType = SVType.DEL,
    consequence: SVConsequence = SVConsequence.WHOLE_GENE_DELETION,
    classification: SVClassification = SVClassification.LIKELY_PATHOGENIC,
    pathogenicity_score: float | None = 0.92,
    evidence_score: float = 0.90,
    syndrome_name: str | None = None,
    gene_overlaps: tuple[GeneOverlap, ...] | None = None,
    genes_affected: int = 0,
    hi_genes_affected: int = 0,
    population_frequency: float | None = None,
    frequency_unknown: bool = True,
) -> ClassifiedSV:
    if gene_overlaps is None:
        gene_overlaps = ()
    sv = StructuralVariant(
        chrom=chrom, start=start, end=end, sv_type=sv_type,
    )
    annotated = AnnotatedSV(
        sv=sv, consequence=consequence,
        gene_overlaps=gene_overlaps, genes_affected=genes_affected,
        hi_genes_affected=hi_genes_affected,
        population_frequency=population_frequency,
        frequency_unknown=frequency_unknown,
    )
    scored = ScoredSV(annotated=annotated, pathogenicity_score=pathogenicity_score)
    return ClassifiedSV(
        scored=scored,
        classification=classification,
        evidence_score=evidence_score,
        syndrome_name=syndrome_name,
    )


class TestFormatSize:
    def test_megabase(self) -> None:
        assert "Mb" in _format_size(2_500_000)

    def test_kilobase(self) -> None:
        assert "kb" in _format_size(150_000)

    def test_base_pairs(self) -> None:
        assert "bp" in _format_size(500)


class TestSVReportBuilderEmpty:
    def test_empty_input_returns_empty_section(self) -> None:
        builder = SVReportBuilder()
        section = builder.build([])
        assert section.findings == []
        assert section.summary.total_svs == 0
        assert section.narratives == []


class TestSVReportBuilderSummary:
    def test_summary_counts_classifications(self) -> None:
        builder = SVReportBuilder()
        results = [
            _make_classified(classification=SVClassification.PATHOGENIC),
            _make_classified(classification=SVClassification.LIKELY_PATHOGENIC),
            _make_classified(classification=SVClassification.VUS),
            _make_classified(classification=SVClassification.VUS),
        ]
        section = builder.build(results)

        assert section.summary.total_svs == 4
        assert section.summary.pathogenic_count == 1
        assert section.summary.likely_pathogenic_count == 1
        assert section.summary.vus_count == 2

    def test_summary_tracks_syndromes(self) -> None:
        builder = SVReportBuilder()
        results = [
            _make_classified(syndrome_name="22q11.2 deletion syndrome"),
            _make_classified(syndrome_name=None),
        ]
        section = builder.build(results)
        assert "22q11.2 deletion syndrome" in section.summary.syndromes_identified

    def test_summary_tracks_hi_genes(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="TBX1", gene_chrom="chr22",
            gene_start=100, gene_end=500,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=9, total_exons=9,
            is_haploinsufficient=True, hi_score=3.0, ts_score=None,
        )
        builder = SVReportBuilder()
        results = [_make_classified(gene_overlaps=(overlap,), hi_genes_affected=1)]
        section = builder.build(results)
        assert "TBX1" in section.summary.hi_genes_affected


class TestSVReportBuilderFindings:
    def test_findings_sorted_by_severity(self) -> None:
        builder = SVReportBuilder()
        vus = _make_classified(classification=SVClassification.VUS, evidence_score=0.3)
        pathogenic = _make_classified(classification=SVClassification.PATHOGENIC, evidence_score=1.0)
        section = builder.build([vus, pathogenic])

        assert section.findings[0].classification == "Pathogenic"
        assert section.findings[1].classification == "VUS"

    def test_findings_row_fields(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="BRCA1", gene_chrom="chr17",
            gene_start=100, gene_end=500,
            overlap_fraction=0.9, is_whole_gene=True,
            exons_affected=23, total_exons=23,
            hi_score=3.0, ts_score=None,
        )
        builder = SVReportBuilder()
        results = [_make_classified(
            gene_overlaps=(overlap,),
            genes_affected=1,
            syndrome_name="Test syndrome",
        )]
        section = builder.build(results)
        row = section.findings[0]

        assert row.sv_type == "DEL"
        assert "chr22" in row.coordinates
        assert row.consequence == "Whole_Gene_Deletion"
        assert row.syndrome == "Test syndrome"
        assert "BRCA1" in row.genes


class TestSVReportBuilderNarrative:
    def test_narrative_includes_sv_description(self) -> None:
        builder = SVReportBuilder()
        results = [_make_classified()]
        section = builder.build(results)
        assert "DEL" in section.narratives[0]
        assert "chr22" in section.narratives[0]

    def test_narrative_mentions_intergenic(self) -> None:
        builder = SVReportBuilder()
        results = [_make_classified(
            consequence=SVConsequence.INTERGENIC,
            genes_affected=0,
            gene_overlaps=(),
        )]
        section = builder.build(results)
        assert "does not overlap" in section.narratives[0]

    def test_narrative_mentions_single_gene_overlap(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="TP53", gene_chrom="chr17",
            gene_start=100, gene_end=500,
            overlap_fraction=0.65, is_whole_gene=False,
            exons_affected=7, total_exons=11,
            hi_score=None, ts_score=None,
        )
        builder = SVReportBuilder()
        results = [_make_classified(gene_overlaps=(overlap,), genes_affected=1)]
        section = builder.build(results)
        assert "TP53" in section.narratives[0]
        assert "65%" in section.narratives[0]

    def test_narrative_mentions_syndrome(self) -> None:
        builder = SVReportBuilder()
        results = [_make_classified(syndrome_name="DiGeorge")]
        section = builder.build(results)
        assert "DiGeorge" in section.narratives[0]

    def test_narrative_mentions_frequency_absent(self) -> None:
        builder = SVReportBuilder()
        results = [_make_classified(frequency_unknown=True)]
        section = builder.build(results)
        assert "gnomAD-SV" in section.narratives[0]

    def test_narrative_mentions_population_frequency(self) -> None:
        builder = SVReportBuilder()
        results = [_make_classified(
            population_frequency=0.003,
            frequency_unknown=False,
        )]
        section = builder.build(results)
        assert "0.300%" in section.narratives[0]

    def test_narrative_escapes_html(self) -> None:
        overlap = GeneOverlap(
            gene_symbol="<script>alert</script>", gene_chrom="chr1",
            gene_start=100, gene_end=500,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=5, total_exons=5,
            hi_score=None, ts_score=None,
        )
        builder = SVReportBuilder()
        results = [_make_classified(gene_overlaps=(overlap,), genes_affected=1)]
        section = builder.build(results)
        # HTML should be escaped
        assert "<script>" not in section.narratives[0]
        assert "&lt;script&gt;" in section.narratives[0]
