"""Unit tests for SNV + SV combined findings merger."""

from __future__ import annotations

from unittest.mock import MagicMock

from vartriage.models.variant import ACMGClassification, FunctionalConsequence
from vartriage.structural.combiner import (
    _serialize_sv,
    merge_findings,
)
from vartriage.structural.models import (
    AnnotatedSV,
    ClassifiedSV,
    GeneOverlap,
    ScoredSV,
    StructuralVariant,
    SVClassification,
    SVConsequence,
    SVType,
)


def _make_snv_classified(
    classification: ACMGClassification = ACMGClassification.VUS,
    score: float = 0.5,
) -> MagicMock:
    """Build a mock ClassifiedVariant (SNV) with needed attributes."""
    mock = MagicMock()
    mock.classification = classification
    mock.scored.prioritization_score = score
    mock.scored.annotated.variant.chrom = "chr17"
    mock.scored.annotated.variant.pos = 7577120
    mock.scored.annotated.variant.ref = "G"
    mock.scored.annotated.variant.alt = "A"
    mock.scored.annotated.gene_name = "TP53"
    mock.scored.annotated.consequence = FunctionalConsequence.MISSENSE
    mock.evidence_tags = frozenset()
    return mock


def _make_sv_classified(
    classification: SVClassification = SVClassification.VUS,
    score: float = 0.6,
) -> ClassifiedSV:
    sv = StructuralVariant(
        chrom="chr22",
        start=18916842,
        end=21465659,
        sv_type=SVType.DEL,
    )
    overlap = GeneOverlap(
        gene_symbol="TBX1",
        gene_chrom="chr22",
        gene_start=19744226,
        gene_end=19771115,
        overlap_fraction=1.0,
        is_whole_gene=True,
        exons_affected=9,
        total_exons=9,
        hi_score=3.0,
        ts_score=None,
    )
    annotated = AnnotatedSV(
        sv=sv,
        consequence=SVConsequence.WHOLE_GENE_DELETION,
        gene_overlaps=(overlap,),
        genes_affected=1,
        hi_genes_affected=1,
    )
    scored = ScoredSV(annotated=annotated, pathogenicity_score=score)
    return ClassifiedSV(
        scored=scored,
        classification=classification,
        evidence_score=0.9,
    )


class TestMergeFindings:
    def test_empty_inputs_return_empty(self) -> None:
        result = merge_findings([], [])
        assert result == []

    def test_snv_only(self) -> None:
        snvs = [_make_snv_classified(ACMGClassification.PATHOGENIC, 0.95)]
        result = merge_findings(snvs, [])
        assert len(result) == 1
        assert result[0].variant_type == "SNV"
        assert result[0].tier == 0

    def test_sv_only(self) -> None:
        svs = [_make_sv_classified(SVClassification.LIKELY_PATHOGENIC, 0.85)]
        result = merge_findings([], svs)
        assert len(result) == 1
        assert result[0].variant_type == "SV"
        assert result[0].tier == 1

    def test_sorted_by_tier_then_score(self) -> None:
        snvs = [
            _make_snv_classified(ACMGClassification.VUS, 0.4),
            _make_snv_classified(ACMGClassification.PATHOGENIC, 0.99),
        ]
        svs = [
            _make_sv_classified(SVClassification.LIKELY_PATHOGENIC, 0.88),
        ]
        result = merge_findings(snvs, svs)

        # Pathogenic (tier=0) first, then LP (tier=1), then VUS (tier=2)
        assert result[0].tier == 0
        assert result[1].tier == 1
        assert result[2].tier == 2

    def test_within_same_tier_sorted_by_score_desc(self) -> None:
        snvs = [
            _make_snv_classified(ACMGClassification.VUS, 0.3),
            _make_snv_classified(ACMGClassification.VUS, 0.7),
        ]
        result = merge_findings(snvs, [])
        assert result[0].score > result[1].score


class TestSerializeSV:
    def test_sv_serialization_fields(self) -> None:
        classified = _make_sv_classified()
        data = _serialize_sv(classified)
        assert data["type"] == "SV"
        assert data["chrom"] == "chr22"
        assert data["sv_type"] == "DEL"
        assert data["classification"] == "VUS"
        assert "TBX1" in data["genes"]
        assert data["genes_affected"] == 1
