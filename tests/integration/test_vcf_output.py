"""Integration test for VCF output format through the pipeline.

Creates a tiny VCF, runs the Pipeline with output_format="vcf",
and verifies the output is correctly bgzipped, indexed, and
annotated with VARTRIAGE_* INFO headers.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from vartriage.models.config import (
    PipelineConfig,
    ReportConfig,
)
from vartriage.pipeline import Pipeline


def _write_tiny_vcf(path: Path) -> Path:
    """Create a minimal bgzipped VCF with two variants using pysam."""
    header = pysam.VariantHeader()
    header.add_sample("SAMPLE")
    header.add_line(
        '##FORMAT=<ID=GT,Number=1,Type=String,'
        'Description="Genotype">'
    )
    header.add_line("##contig=<ID=chr1,length=1000000>")

    with pysam.VariantFile(str(path), "wz", header=header) as vcf:
        rec1 = vcf.new_record(
            contig="chr1",
            start=999,  # 0-based, pos=1000
            alleles=("A", "T"),
        )
        rec1.samples["SAMPLE"]["GT"] = (0, 1)
        vcf.write(rec1)

        rec2 = vcf.new_record(
            contig="chr1",
            start=1999,  # 0-based, pos=2000
            alleles=("G", "C"),
        )
        rec2.samples["SAMPLE"]["GT"] = (0, 1)
        vcf.write(rec2)

    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


class TestVCFOutputIntegration:
    """Pipeline produces valid bgzipped VCF with tabix index."""

    @pytest.fixture()
    def pipeline_vcf(self, tmp_path: Path) -> tuple[Path, Path]:
        """Run the pipeline with VCF output and return paths.

        Returns (output_vcf_path, source_vcf_path).
        """
        source_vcf = _write_tiny_vcf(tmp_path / "input.vcf.gz")
        output_path = tmp_path / "result.vcf.gz"

        config = PipelineConfig(
            vcf_path=source_vcf,
            output_path=output_path,
            report=ReportConfig(output_format="vcf"),
        )
        pipeline = Pipeline(config)
        pipeline.run()

        return output_path, source_vcf

    def test_output_is_bgzipped(
        self, pipeline_vcf: tuple[Path, Path]
    ) -> None:
        """Output file exists at expected path and is bgzipped."""
        output_path, _ = pipeline_vcf
        assert output_path.exists()

        # bgzip files start with the gzip magic bytes 1f 8b
        magic = output_path.read_bytes()[:2]
        assert magic == b"\x1f\x8b"

    def test_tbi_index_exists(
        self, pipeline_vcf: tuple[Path, Path]
    ) -> None:
        """A .tbi tabix index file exists alongside the VCF."""
        output_path, _ = pipeline_vcf
        tbi_path = Path(str(output_path) + ".tbi")
        assert tbi_path.exists()
        assert tbi_path.stat().st_size > 0

    def test_vartriage_headers_present(
        self, pipeline_vcf: tuple[Path, Path]
    ) -> None:
        """Output VCF header declares all VARTRIAGE_* INFO fields."""
        output_path, _ = pipeline_vcf

        with pysam.VariantFile(str(output_path)) as vcf:
            info_ids = set(vcf.header.info)

        expected = {
            "VARTRIAGE_CONSEQUENCE",
            "VARTRIAGE_AF",
            "VARTRIAGE_RANK",
            "VARTRIAGE_ACMG",
            "VARTRIAGE_TAGS",
        }
        assert expected.issubset(info_ids)

    def test_opens_with_pysam(
        self, pipeline_vcf: tuple[Path, Path]
    ) -> None:
        """Output can be opened and iterated with pysam."""
        output_path, _ = pipeline_vcf

        with pysam.VariantFile(str(output_path)) as vcf:
            records = list(vcf)

        # Should have the same two records we put into the source
        assert len(records) == 2

    def test_records_have_correct_positions(
        self, pipeline_vcf: tuple[Path, Path]
    ) -> None:
        """Output records preserve the original genomic positions."""
        output_path, _ = pipeline_vcf

        with pysam.VariantFile(str(output_path)) as vcf:
            records = list(vcf)

        positions = [rec.pos for rec in records]
        assert 1000 in positions
        assert 2000 in positions
