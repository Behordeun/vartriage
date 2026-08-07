"""Integration tests for Pipeline.run with mitochondrial variant wiring.

Exercises the top-level Pipeline to verify:
- chrM variants are excluded from the nuclear classification stream
- pipeline.mito_results is populated when mito is enabled and chrM present
- JSON output includes mitochondrial_findings when mito variants exist
- CSV output includes mitochondrial section when mito variants exist
- Legacy flat formats are preserved when mito is disabled or no chrM present
"""

from __future__ import annotations

import json
from pathlib import Path

from vartriage.mito.config import MitoConfig
from vartriage.models.config import PipelineConfig, ReportConfig
from vartriage.pipeline import Pipeline

_VCF_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

_VCF_HEADER_WITH_SAMPLE = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""


def _write_mixed_vcf(tmp_dir: Path) -> Path:
    """Write a VCF containing both nuclear (chr1) and mitochondrial (chrM) variants."""
    vcf_path = tmp_dir / "mixed.vcf"
    lines = [
        # Nuclear variants
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=30",
        "chr1\t200\t.\tG\tC\t60\tPASS\tDP=40",
        "chr1\t300\t.\tC\tA\t70\tPASS\tDP=50",
        # Mitochondrial variants
        "chrM\t73\t.\tA\tG\t200\tPASS\tDP=1000",
        "chrM\t11778\t.\tG\tA\t200\tPASS\tDP=1200",
    ]
    content = _VCF_HEADER + "\n".join(lines) + "\n"
    vcf_path.write_text(content, encoding="utf-8")
    return vcf_path


def _write_mixed_vcf_with_sample(tmp_dir: Path) -> Path:
    """Write a VCF with FORMAT/AD fields for heteroplasmy extraction."""
    vcf_path = tmp_dir / "mixed_sample.vcf"
    lines = [
        # Nuclear
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=30\tGT:AD\t0/1:15,15",
        "chr1\t200\t.\tG\tC\t60\tPASS\tDP=40\tGT:AD\t0/1:20,20",
        # chrM with heteroplasmy data via AD
        "chrM\t73\t.\tA\tG\t200\tPASS\tDP=1000\tGT:AD\t0/1:50,950",
        "chrM\t11778\t.\tG\tA\t200\tPASS\tDP=1200\tGT:AD\t0/1:20,980",
    ]
    content = _VCF_HEADER_WITH_SAMPLE + "\n".join(lines) + "\n"
    vcf_path.write_text(content, encoding="utf-8")
    return vcf_path


def _write_nuclear_only_vcf(tmp_dir: Path) -> Path:
    """Write a VCF with nuclear variants only (no chrM)."""
    vcf_path = tmp_dir / "nuclear_only.vcf"
    lines = [
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=30",
        "chr1\t200\t.\tG\tC\t60\tPASS\tDP=40",
        "chr1\t300\t.\tC\tA\t70\tPASS\tDP=50",
    ]
    content = _VCF_HEADER + "\n".join(lines) + "\n"
    vcf_path.write_text(content, encoding="utf-8")
    return vcf_path


class TestPipelineMitoRouting:
    """Verify chrM variants are routed to MitochondrialPipeline."""

    def test_mito_results_populated_when_chrm_present(self, tmp_path: Path):
        vcf_path = _write_mixed_vcf(tmp_path)
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
        )
        pipeline = Pipeline(config)
        pipeline.run()

        assert pipeline.mito_results is not None
        assert len(pipeline.mito_results) > 0

    def test_mito_results_none_when_no_chrm(self, tmp_path: Path):
        vcf_path = _write_nuclear_only_vcf(tmp_path)
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
        )
        pipeline = Pipeline(config)
        pipeline.run()

        # No chrM variants -> mito results should be None or empty
        mito = pipeline.mito_results
        assert mito is None or len(mito) == 0

    def test_mito_results_none_when_disabled(self, tmp_path: Path):
        vcf_path = _write_mixed_vcf(tmp_path)
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            mito=MitoConfig(enabled=False),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        assert pipeline.mito_results is None


class TestJsonOutputWithMito:
    """Verify JSON output format with and without mitochondrial findings."""

    def test_json_includes_mitochondrial_findings(self, tmp_path: Path):
        vcf_path = _write_mixed_vcf(tmp_path)
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="json"),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        report = json.loads(output_path.read_text())

        assert "mitochondrial_findings" in report
        assert isinstance(report["mitochondrial_findings"], list)
        assert len(report["mitochondrial_findings"]) > 0

        assert "variants" in report
        assert isinstance(report["variants"], list)

        # Nuclear variants should not contain chrM
        for var in report["variants"]:
            assert var["chromosome"] != "chrM"

    def test_json_legacy_format_when_no_chrm(self, tmp_path: Path):
        vcf_path = _write_nuclear_only_vcf(tmp_path)
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="json"),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        report = json.loads(output_path.read_text())

        # Legacy flat array format (no mitochondrial_findings key)
        assert isinstance(report, list)
        assert "mitochondrial_findings" not in (
            report if isinstance(report, dict) else {}
        )

    def test_json_legacy_format_when_mito_disabled(self, tmp_path: Path):
        vcf_path = _write_mixed_vcf(tmp_path)
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="json"),
            mito=MitoConfig(enabled=False),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        report = json.loads(output_path.read_text())

        # With mito disabled, output is flat array (all variants, including
        # chrM since they pass through the nuclear path unclassified)
        assert isinstance(report, list)


class TestCsvOutputWithMito:
    """Verify CSV output format with and without mitochondrial findings."""

    def test_csv_includes_mitochondrial_section(self, tmp_path: Path):
        vcf_path = _write_mixed_vcf(tmp_path)
        output_path = tmp_path / "output.csv"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="csv"),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        content = output_path.read_text()
        assert "Mitochondrial" in content

    def test_csv_legacy_format_when_no_chrm(self, tmp_path: Path):
        vcf_path = _write_nuclear_only_vcf(tmp_path)
        output_path = tmp_path / "output.csv"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="csv"),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        content = output_path.read_text()
        assert "Mitochondrial" not in content
        # Should still have the nuclear header row
        assert "chromosome" in content

    def test_csv_legacy_format_when_mito_disabled(self, tmp_path: Path):
        vcf_path = _write_mixed_vcf(tmp_path)
        output_path = tmp_path / "output.csv"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="csv"),
            mito=MitoConfig(enabled=False),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        content = output_path.read_text()
        assert "Mitochondrial" not in content
