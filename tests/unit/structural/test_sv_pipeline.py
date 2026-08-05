"""Unit tests for SV triage pipeline orchestration."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from vartriage.structural.config import SVTriageConfig
from vartriage.structural.models import (
    AnnotatedSV,
    ClassifiedSV,
    GeneOverlap,
    ScoredSV,
    StructuralVariant,
    SVClassification,
    SVConsequence,
    SVEvidenceCategory,
    SVType,
)
from vartriage.structural.pipeline import SVTriagePipeline


def _make_classified(
    classification: SVClassification = SVClassification.LIKELY_PATHOGENIC,
    pathogenicity_score: float | None = 0.85,
    chrom: str = "chr22",
    start: int = 18916842,
    end: int = 21465659,
) -> ClassifiedSV:
    sv = StructuralVariant(
        chrom=chrom,
        start=start,
        end=end,
        sv_type=SVType.DEL,
        id="test_del",
        qual=999.0,
        filter_status="PASS",
        alt="<DEL>",
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
        is_haploinsufficient=True,
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
    scored = ScoredSV(
        annotated=annotated,
        pathogenicity_score=pathogenicity_score,
        dosage_score=1.0,
        size_score=0.9,
        frequency_score=1.0,
    )
    return ClassifiedSV(
        scored=scored,
        classification=classification,
        evidence_categories=frozenset({SVEvidenceCategory.GENE_FULLY_CONTAINED}),
        evidence_score=0.90,
        syndrome_name="22q11.2 deletion syndrome",
    )


class TestPipelineValidation:
    def test_missing_vcf_raises(self, tmp_path: Path) -> None:
        config = SVTriageConfig(
            vcf_path=tmp_path / "missing.vcf",
            output_path=tmp_path / "out.json",
        )
        with pytest.raises(FileNotFoundError, match="Input VCF not found"):
            SVTriagePipeline(config)

    def test_missing_gene_annotation_raises(self, tmp_path: Path) -> None:
        vcf = tmp_path / "input.vcf"
        vcf.touch()
        config = SVTriageConfig(
            vcf_path=vcf,
            output_path=tmp_path / "out.json",
            gene_annotation_path=tmp_path / "missing.gtf",
        )
        with pytest.raises(FileNotFoundError, match="Gene annotation"):
            SVTriagePipeline(config)

    def test_missing_dosage_file_raises(self, tmp_path: Path) -> None:
        vcf = tmp_path / "input.vcf"
        vcf.touch()
        config = SVTriageConfig(
            vcf_path=vcf,
            output_path=tmp_path / "out.json",
            dosage_sensitivity_path=tmp_path / "missing.tsv",
        )
        with pytest.raises(FileNotFoundError, match="Dosage sensitivity"):
            SVTriagePipeline(config)

    def test_missing_gnomad_file_raises(self, tmp_path: Path) -> None:
        vcf = tmp_path / "input.vcf"
        vcf.touch()
        config = SVTriageConfig(
            vcf_path=vcf,
            output_path=tmp_path / "out.json",
            gnomad_sv_path=tmp_path / "missing.bed",
        )
        with pytest.raises(FileNotFoundError, match="gnomAD-SV"):
            SVTriagePipeline(config)


class TestCollectResults:
    def _make_pipeline(self, tmp_path: Path) -> SVTriagePipeline:
        vcf = tmp_path / "input.vcf"
        vcf.touch()
        config = SVTriageConfig(
            vcf_path=vcf,
            output_path=tmp_path / "out.json",
            include_benign=False,
        )
        pipeline = object.__new__(SVTriagePipeline)
        pipeline._config = config
        return pipeline

    def test_filters_benign_when_include_benign_false(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        variants = [
            _make_classified(SVClassification.LIKELY_PATHOGENIC),
            _make_classified(SVClassification.BENIGN),
            _make_classified(SVClassification.LIKELY_BENIGN),
        ]
        results = pipeline._collect_results(iter(variants))
        assert len(results) == 1
        assert results[0].classification == SVClassification.LIKELY_PATHOGENIC

    def test_sorts_by_score_descending(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        low = _make_classified(SVClassification.VUS, pathogenicity_score=0.3)
        high = _make_classified(SVClassification.VUS, pathogenicity_score=0.9)
        results = pipeline._collect_results(iter([low, high]))
        assert (
            results[0].scored.pathogenicity_score
            > results[1].scored.pathogenicity_score
        )

    def test_none_scores_sorted_last(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        scored = _make_classified(SVClassification.VUS, pathogenicity_score=0.5)
        unscored = _make_classified(SVClassification.VUS, pathogenicity_score=None)
        results = pipeline._collect_results(iter([unscored, scored]))
        assert results[0].scored.pathogenicity_score is not None
        assert results[1].scored.pathogenicity_score is None


class TestWriteCSV:
    def test_csv_output_has_header_and_rows(self, tmp_path: Path) -> None:
        vcf = tmp_path / "input.vcf"
        vcf.touch()
        config = SVTriageConfig(
            vcf_path=vcf,
            output_path=tmp_path / "out.csv",
            output_format="csv",
        )
        pipeline = object.__new__(SVTriagePipeline)
        pipeline._config = config

        results = [_make_classified()]
        pipeline._write_report(results)

        output = tmp_path / "out.csv"
        assert output.exists()

        with open(output) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert rows[0]["chrom"] == "chr22"
        assert rows[0]["sv_type"] == "DEL"
        assert rows[0]["classification"] == "Likely_Pathogenic"
        assert "TBX1" in rows[0]["gene_symbols"]


class TestLoadRegions:
    def _make_pipeline(self, tmp_path: Path) -> SVTriagePipeline:
        vcf = tmp_path / "input.vcf"
        vcf.touch()
        config = SVTriageConfig(vcf_path=vcf, output_path=tmp_path / "out.json")
        pipeline = object.__new__(SVTriagePipeline)
        pipeline._config = config
        return pipeline

    def test_load_regions_from_bed(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        bed = tmp_path / "regions.bed"
        bed.write_text(
            "# header\nchr22\t18916842\t21465659\nchr15\t20143000\t20570000\n"
        )
        regions = pipeline._load_regions(bed)
        assert len(regions) == 2
        assert regions[0] == ("chr22", 18916842, 21465659)

    def test_load_regions_skips_track_lines(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        bed = tmp_path / "regions.bed"
        bed.write_text("track name=test\nchr1\t100\t200\n")
        regions = pipeline._load_regions(bed)
        assert len(regions) == 1

    def test_load_regions_returns_empty_for_none(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        regions = pipeline._load_regions(None)
        assert regions == []

    def test_load_regions_with_names(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        bed = tmp_path / "regions.bed"
        bed.write_text("chr22\t18916842\t21465659\tDiGeorge syndrome\n")
        regions, names = pipeline._load_regions_with_names(bed)
        assert len(regions) == 1
        assert names[("chr22", 18916842, 21465659)] == "DiGeorge syndrome"

    def test_load_regions_with_names_handles_no_name(self, tmp_path: Path) -> None:
        pipeline = self._make_pipeline(tmp_path)
        bed = tmp_path / "regions.bed"
        bed.write_text("chr1\t100\t500\n")
        regions, names = pipeline._load_regions_with_names(bed)
        assert len(regions) == 1
        assert len(names) == 0
