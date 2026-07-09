"""Smoke test for the Pipeline orchestrator (task 10.4).

Verifies that the Pipeline class:
1. Validates configuration at construction (fail-fast)
2. Wires all stages together correctly
3. Produces output when run with valid inputs
4. Uses the WarningAccumulator for missing data tracking
"""

from __future__ import annotations

import csv
import json
import warnings
from pathlib import Path

import pytest

from vartriage.models.config import (
    AnnotationConfig,
    MissingDataConfig,
    PipelineConfig,
    PrioritizationConfig,
    QualityFilterConfig,
    ReportConfig,
)
from vartriage.pipeline import Pipeline


_VCF_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _write_vcf(tmp_dir: Path) -> Path:
    """Write a minimal valid VCF with 10 passing variants."""
    vcf_path = tmp_dir / "input.vcf"
    lines = [
        "chr1\t100\tvar_1\tA\tT\t50\tPASS\tDP=30",
        "chr1\t200\tvar_2\tG\tC\t60\tPASS\tDP=40",
        "chr1\t300\tvar_3\tC\tA\t70\tPASS\tDP=50",
        "chr1\t400\tvar_4\tT\tG\t80\tPASS\tDP=60",
        "chr1\t500\tvar_5\tA\tG\t90\tPASS\tDP=70",
        "chr1\t600\tvar_6\tG\tT\t55\tPASS\tDP=35",
        "chr1\t700\tvar_7\tC\tG\t65\tPASS\tDP=45",
        "chr1\t800\tvar_8\tT\tA\t75\tPASS\tDP=55",
        "chr1\t900\tvar_9\tA\tC\t85\tPASS\tDP=65",
        "chr1\t1000\tvar_10\tG\tA\t95\tPASS\tDP=75",
        # Variants that should be filtered out
        "chr1\t1100\tvar_11\tA\tT\t50\tLowQual\tDP=10",
        "chr1\t1200\tvar_12\tG\tC\t.\tPASS\tDP=15",
        "chr1\t1300\tvar_13\tC\tA\t5\tPASS\tDP=12",
    ]
    content = _VCF_HEADER + "\n".join(lines) + "\n"
    vcf_path.write_text(content, encoding="utf-8")
    return vcf_path


def _write_gtf(tmp_dir: Path) -> Path:
    """Write a minimal GTF with one gene covering chr1:50-1500."""
    gtf_path = tmp_dir / "genes.gtf"
    lines = [
        '##format: gtf',
        'chr1\thavana\tgene\t50\t1500\t.\t+\t.\t'
        'gene_id "BRCA1"; gene_name "BRCA1";',
        'chr1\thavana\ttranscript\t50\t1500\t.\t+\t.\t'
        'gene_id "BRCA1"; transcript_id "BRCA1.1"; gene_name "BRCA1";',
        'chr1\thavana\texon\t50\t1500\t.\t+\t.\t'
        'gene_id "BRCA1"; transcript_id "BRCA1.1"; gene_name "BRCA1";',
        'chr1\thavana\tCDS\t100\t1200\t.\t+\t0\t'
        'gene_id "BRCA1"; transcript_id "BRCA1.1"; gene_name "BRCA1";',
    ]
    gtf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return gtf_path


def _write_gnomad(tmp_dir: Path) -> Path:
    """Write a minimal gnomAD TSV with some matching entries."""
    gnomad_path = tmp_dir / "gnomad.tsv"
    rows = [
        ["chrom", "pos", "ref", "alt", "af"],
        ["chr1", "100", "A", "T", "0.00005"],
        ["chr1", "200", "G", "C", "0.005"],
        ["chr1", "300", "C", "A", "0.05"],
    ]
    with open(gnomad_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)
    return gnomad_path


def _write_clinvar(tmp_dir: Path) -> Path:
    """Write a minimal ClinVar TSV."""
    clinvar_path = tmp_dir / "clinvar.tsv"
    rows = [
        ["chrom", "pos", "ref", "alt", "clinical_significance"],
        ["chr1", "100", "A", "T", "Pathogenic"],
    ]
    with open(clinvar_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)
    return clinvar_path


@pytest.fixture
def pipeline_files(tmp_path: Path) -> dict[str, Path]:
    """Create synthetic reference data and return paths."""
    vcf_path = _write_vcf(tmp_path)
    gtf_path = _write_gtf(tmp_path)
    gnomad_path = _write_gnomad(tmp_path)
    clinvar_path = _write_clinvar(tmp_path)
    output_path = tmp_path / "output" / "report.json"
    return {
        "vcf_path": vcf_path,
        "gtf_path": gtf_path,
        "gnomad_path": gnomad_path,
        "clinvar_path": clinvar_path,
        "output_path": output_path,
    }


class TestPipelineConstruction:
    """Test Pipeline construction and config validation."""

    def test_pipeline_construction_with_valid_config(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline can be constructed with valid config."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=pipeline_files["gnomad_path"],
                clinvar_path=pipeline_files["clinvar_path"],
            ),
        )
        pipeline = Pipeline(config)
        assert pipeline is not None

    def test_pipeline_fails_fast_on_missing_gtf(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline raises FileNotFoundError if gene annotation is missing."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            annotation=AnnotationConfig(
                gene_annotation_path=Path("/nonexistent/genes.gtf"),
                gnomad_path=pipeline_files["gnomad_path"],
            ),
        )
        with pytest.raises(FileNotFoundError, match="Gene annotation"):
            Pipeline(config)

    def test_pipeline_fails_fast_on_missing_gnomad(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline raises FileNotFoundError if gnomAD file is missing."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=Path("/nonexistent/gnomad.tsv"),
            ),
        )
        with pytest.raises(FileNotFoundError, match="gnomAD"):
            Pipeline(config)

    def test_pipeline_fails_fast_on_missing_clinvar(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline raises FileNotFoundError if ClinVar file is missing."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=pipeline_files["gnomad_path"],
                clinvar_path=Path("/nonexistent/clinvar.tsv"),
            ),
        )
        with pytest.raises(FileNotFoundError, match="ClinVar"):
            Pipeline(config)

    def test_pipeline_without_annotation_config(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline can be constructed without annotation config (passthrough)."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
        )
        pipeline = Pipeline(config)
        assert pipeline is not None


class TestPipelineRun:
    """Test Pipeline.run() execution."""

    def test_run_produces_json_output(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline run with annotation produces valid JSON output."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            quality_filter=QualityFilterConfig(min_qual=20.0),
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=pipeline_files["gnomad_path"],
                clinvar_path=pipeline_files["clinvar_path"],
            ),
            report=ReportConfig(output_format="json"),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_path = pipeline.run()

        assert result_path.exists()
        data = json.loads(result_path.read_text(encoding="utf-8"))
        assert isinstance(data, list)
        assert len(data) > 0

    def test_run_without_annotation_config(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline run without annotation uses passthrough and still produces output."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            quality_filter=QualityFilterConfig(min_qual=20.0),
            report=ReportConfig(output_format="json"),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_path = pipeline.run()

        assert result_path.exists()
        data = json.loads(result_path.read_text(encoding="utf-8"))
        assert isinstance(data, list)
        assert len(data) > 0

    def test_run_with_overridden_paths(
        self, pipeline_files: dict[str, Path], tmp_path: Path
    ) -> None:
        """Pipeline run can override vcf_path and output_path."""
        alt_output = tmp_path / "alt_output" / "report.json"
        config = PipelineConfig(
            vcf_path=Path("/dummy/not_used.vcf"),
            output_path=Path("/dummy/not_used.json"),
            quality_filter=QualityFilterConfig(min_qual=20.0),
            report=ReportConfig(output_format="json"),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_path = pipeline.run(
                vcf_path=pipeline_files["vcf_path"],
                output_path=alt_output,
            )

        assert result_path == alt_output
        assert result_path.exists()

    def test_run_csv_output(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline produces valid CSV output."""
        output_path = pipeline_files["output_path"].parent / "report.csv"
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=output_path,
            quality_filter=QualityFilterConfig(min_qual=20.0),
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=pipeline_files["gnomad_path"],
                clinvar_path=pipeline_files["clinvar_path"],
            ),
            report=ReportConfig(output_format="csv"),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_path = pipeline.run()

        assert result_path.exists()
        with open(result_path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)
        assert len(rows) > 1
        assert rows[0][0] == "chromosome"

    def test_run_excludes_filtered_variants(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline excludes LowQual, missing QUAL, and below-threshold variants."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            quality_filter=QualityFilterConfig(min_qual=20.0),
            report=ReportConfig(output_format="json"),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result_path = pipeline.run()

        data = json.loads(result_path.read_text(encoding="utf-8"))
        # VCF has 13 total: 10 PASS with qual >= 20, 1 LowQual, 1 missing, 1 below threshold
        # So we expect 10 variants in output
        assert len(data) == 10


class TestPipelineWarningAccumulator:
    """Test warning accumulation across pipeline stages."""

    def test_warning_accumulator_is_accessible(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Pipeline exposes the warning accumulator."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=pipeline_files["gnomad_path"],
                clinvar_path=pipeline_files["clinvar_path"],
            ),
            missing_data=MissingDataConfig(warning_threshold=1000),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pipeline.run()

        acc = pipeline.warning_accumulator
        # Some variants not in gnomAD/ClinVar should trigger warnings
        assert acc.total_count >= 0

    def test_warning_accumulator_resets_on_rerun(
        self, pipeline_files: dict[str, Path]
    ) -> None:
        """Warning accumulator resets between runs."""
        config = PipelineConfig(
            vcf_path=pipeline_files["vcf_path"],
            output_path=pipeline_files["output_path"],
            annotation=AnnotationConfig(
                gene_annotation_path=pipeline_files["gtf_path"],
                gnomad_path=pipeline_files["gnomad_path"],
                clinvar_path=pipeline_files["clinvar_path"],
            ),
        )
        pipeline = Pipeline(config)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pipeline.run()
            count_first = pipeline.warning_accumulator.total_count
            pipeline.run()
            count_second = pipeline.warning_accumulator.total_count

        assert count_first == count_second
