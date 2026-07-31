"""Integration test: synthetic 22q11.2 deletion through full SV pipeline."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from vartriage.structural.annotator import GeneRecord, SVAnnotator
from vartriage.structural.classifier import SVClassifier
from vartriage.structural.models import (
    SVClassification,
    SVConsequence,
    SVEvidenceCategory,
    SVType,
    StructuralVariant,
)
from vartriage.structural.scoring import SVScorer


def _build_22q11_deletion() -> StructuralVariant:
    """Create a synthetic 22q11.2 deletion (~2.5 Mb, DiGeorge region)."""
    return StructuralVariant(
        chrom="chr22",
        start=18916842,
        end=21465659,
        sv_type=SVType.DEL,
        id="DiGeorge_DEL",
        svlen=-2548817,
        qual=999.0,
        filter_status="PASS",
        alt="<DEL>",
    )


class TestFullPipeline22q11Deletion:
    """End-to-end test: 22q11.2 deletion → Pathogenic/Likely_Pathogenic."""

    def test_22q11_deletion_classified_pathogenic(self) -> None:
        """A 2.5Mb 22q11.2 deletion containing TBX1 (HI gene) should
        score as Pathogenic or Likely_Pathogenic via ClinGen framework."""

        sv = _build_22q11_deletion()

        # Set up annotator with known genes in the 22q11.2 region
        annotator = SVAnnotator(
            reciprocal_overlap=0.5,
            whole_gene_threshold=0.8,
        )
        annotator._genes = {
            "chr22": [
                GeneRecord("TBX1", "chr22", 19744226, 19771115, "+", 9),
                GeneRecord("HIRA", "chr22", 19363027, 19465975, "+", 25),
                GeneRecord("COMT", "chr22", 19941375, 19969975, "+", 6),
                GeneRecord("GP1BB", "chr22", 19698770, 19709640, "+", 2),
                GeneRecord("DGCR8", "chr22", 20068567, 20107369, "+", 14),
            ]
        }
        # TBX1 is haploinsufficient
        from vartriage.structural.annotator import DosageEntry
        annotator._dosage = {
            "TBX1": DosageEntry("TBX1", hi_score=3.0, ts_score=None),
            "HIRA": DosageEntry("HIRA", hi_score=2.0, ts_score=None),
        }
        annotator._sv_database = {}

        # Annotate
        annotated = annotator._annotate_single(sv)

        assert annotated.consequence == SVConsequence.WHOLE_GENE_DELETION
        assert annotated.genes_affected == 5
        assert annotated.hi_genes_affected >= 1

        # Score
        scorer = SVScorer(max_allele_frequency=0.01)
        scored = list(scorer.score(iter([annotated])))[0]

        assert scored.pathogenicity_score is not None
        assert scored.pathogenicity_score > 0.5

        # Classify with known pathogenic region
        pathogenic_regions = [("chr22", 18916842, 21465659)]
        region_names = {
            ("chr22", 18916842, 21465659): "22q11.2 deletion syndrome (DiGeorge)"
        }

        classifier = SVClassifier(
            pathogenic_regions=pathogenic_regions,
            pathogenic_region_names=region_names,
        )
        classified = list(classifier.classify(iter([scored])))[0]

        # The combination of HI gene (Section 1B + 3A) + pathogenic region
        # overlap (Section 2A) should push this to Pathogenic
        assert classified.classification in (
            SVClassification.PATHOGENIC,
            SVClassification.LIKELY_PATHOGENIC,
        )
        assert classified.syndrome_name == "22q11.2 deletion syndrome (DiGeorge)"
        assert SVEvidenceCategory.CONTAINS_ESTABLISHED_HI_GENE in classified.evidence_categories
        assert SVEvidenceCategory.COMPLETE_OVERLAP_PATHOGENIC in classified.evidence_categories

    def test_common_intergenic_sv_classified_vus_or_benign(self) -> None:
        """A common intergenic SV should not be classified pathogenic."""

        sv = StructuralVariant(
            chrom="chr8",
            start=39227735,
            end=39377735,
            sv_type=SVType.DEL,
            id="common_del",
            qual=500.0,
            filter_status="PASS",
            alt="<DEL>",
        )

        annotator = SVAnnotator(
            reciprocal_overlap=0.5,
            whole_gene_threshold=0.8,
        )
        annotator._genes = {}
        annotator._dosage = {}
        # SV is common in population
        from vartriage.structural.annotator import SVFrequencyRecord
        annotator._sv_database = {
            "chr8": [
                SVFrequencyRecord("chr8", 39227735, 39377735, "DEL", 0.05)
            ]
        }

        annotated = annotator._annotate_single(sv)
        assert annotated.consequence == SVConsequence.INTERGENIC
        assert annotated.population_frequency == 0.05

        # Common SVs get filtered by scorer (AF > 0.01)
        scorer = SVScorer(max_allele_frequency=0.01)
        scored = list(scorer.score(iter([annotated])))

        # Should be filtered out entirely by the frequency filter
        assert len(scored) == 0

    def test_pipeline_json_output_format(self, tmp_path: Path) -> None:
        """Verify JSON output structure from the pipeline report writer."""
        from vartriage.structural.pipeline import SVTriagePipeline
        from vartriage.structural.config import SVTriageConfig
        from vartriage.structural.models import (
            AnnotatedSV, ClassifiedSV, ScoredSV, GeneOverlap,
        )

        # Instead of running full pipeline (needs VCF), test the report writer
        pipeline = object.__new__(SVTriagePipeline)
        pipeline._config = SVTriageConfig(
            vcf_path=tmp_path / "dummy.vcf",
            output_path=tmp_path / "output.json",
            output_format="json",
        )

        sv = StructuralVariant(
            chrom="chr22", start=18916842, end=21465659,
            sv_type=SVType.DEL, id="test_del", qual=999.0,
            filter_status="PASS", alt="<DEL>",
        )
        overlap = GeneOverlap(
            gene_symbol="TBX1", gene_chrom="chr22",
            gene_start=19744226, gene_end=19771115,
            overlap_fraction=1.0, is_whole_gene=True,
            exons_affected=9, total_exons=9,
            is_haploinsufficient=True, is_triplosensitive=False,
            hi_score=3.0, ts_score=None,
        )
        annotated = AnnotatedSV(
            sv=sv, consequence=SVConsequence.WHOLE_GENE_DELETION,
            gene_overlaps=(overlap,), genes_affected=1,
            hi_genes_affected=1,
        )
        scored = ScoredSV(
            annotated=annotated, pathogenicity_score=0.85,
            dosage_score=1.0, size_score=0.9, frequency_score=1.0,
        )
        classified = ClassifiedSV(
            scored=scored,
            classification=SVClassification.LIKELY_PATHOGENIC,
            evidence_categories=frozenset({
                SVEvidenceCategory.CONTAINS_ESTABLISHED_HI_GENE,
                SVEvidenceCategory.GENE_FULLY_CONTAINED,
            }),
            evidence_score=0.90,
            syndrome_name="22q11.2 deletion syndrome (DiGeorge)",
        )

        pipeline._write_report([classified])

        output_path = tmp_path / "output.json"
        assert output_path.exists()

        with open(output_path) as fh:
            data = json.load(fh)

        assert data["pipeline"] == "structural_variant_triage"
        assert data["total_variants"] == 1
        assert len(data["variants"]) == 1

        variant = data["variants"][0]
        assert variant["chrom"] == "chr22"
        assert variant["sv_type"] == "DEL"
        assert variant["classification"] == "Likely_Pathogenic"
        assert variant["syndrome_name"] == "22q11.2 deletion syndrome (DiGeorge)"
        assert len(variant["gene_details"]) == 1
        assert variant["gene_details"][0]["symbol"] == "TBX1"
