"""Unit tests for VCFParser."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import pytest

from vartriage.io.exceptions import ParseError
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.variant import Variant


MINIMAL_VCF_CONTENT = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30.0\tPASS\tDP=50
chr1\t200\trs123\tC\tT\t45.5\t.\tDP=100
chr2\t300\t.\tG\tA\t.\tLowQual\tDP=10
"""


def _write_vcf(tmp_path: Path, content: str) -> Path:
    """Write VCF content to a plain .vcf file."""
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(content)
    return vcf_path


def _write_vcf_gz(tmp_path: Path, content: str) -> tuple[Path, Path]:
    """Write bgzipped VCF and create a .tbi index via pysam."""
    import pysam

    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(content)

    gz_path = tmp_path / "test.vcf.gz"
    pysam.tabix_compress(str(vcf_path), str(gz_path))
    pysam.tabix_index(str(gz_path), preset="vcf")

    tbi_path = Path(str(gz_path) + ".tbi")
    return gz_path, tbi_path


class TestVCFParserInit:
    """Tests for VCFParser initialization and validation."""

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        """Raises FileNotFoundError when the file doesn't exist."""
        fake_path = tmp_path / "nonexistent.vcf"
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            VCFParser(fake_path)

    def test_missing_tbi_index_raises(self, tmp_path: Path) -> None:
        """Raises FileNotFoundError when .tbi is missing for .vcf.gz."""
        gz_path = tmp_path / "test.vcf.gz"
        gz_path.write_bytes(b"fake content")
        with pytest.raises(FileNotFoundError, match="Tabix index not found"):
            VCFParser(gz_path)

    def test_valid_vcf_opens_without_error(self, tmp_path: Path) -> None:
        """Opens a valid .vcf file without raising."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        parser = VCFParser(vcf_path)
        parser.close()


class TestVCFParserIteration:
    """Tests for iterating over VCF records."""

    def test_yields_correct_number_of_variants(self, tmp_path: Path) -> None:
        """Yields one Variant per data line."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        with VCFParser(vcf_path) as parser:
            variants = list(parser)
        assert len(variants) == 3

    def test_variant_fields_match_vcf_content(self, tmp_path: Path) -> None:
        """Parsed Variant fields match the source data line values."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        with VCFParser(vcf_path) as parser:
            variants = list(parser)

        first = variants[0]
        assert isinstance(first, Variant)
        assert first.chrom == "chr1"
        assert first.pos == 100
        assert first.id is None
        assert first.ref == "A"
        assert first.alt == "G"
        assert first.qual == pytest.approx(30.0)
        assert first.filter_status == "PASS"
        assert first.info.get("DP") == 50

    def test_variant_id_parsed(self, tmp_path: Path) -> None:
        """Parses variant ID when present."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        with VCFParser(vcf_path) as parser:
            variants = list(parser)

        second = variants[1]
        assert second.id == "rs123"
        assert second.qual == pytest.approx(45.5)

    def test_missing_qual_yields_none(self, tmp_path: Path) -> None:
        """Missing QUAL field (.) produces qual=None."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        with VCFParser(vcf_path) as parser:
            variants = list(parser)

        third = variants[2]
        assert third.qual is None

    def test_filter_status_preserved(self, tmp_path: Path) -> None:
        """FILTER field value is correctly extracted."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        with VCFParser(vcf_path) as parser:
            variants = list(parser)

        assert variants[0].filter_status == "PASS"
        assert variants[1].filter_status == "."
        assert variants[2].filter_status == "LowQual"

    def test_compressed_vcf_with_index(self, tmp_path: Path) -> None:
        """Parses .vcf.gz files with tabix index."""
        gz_path, _ = _write_vcf_gz(tmp_path, MINIMAL_VCF_CONTENT)
        with VCFParser(gz_path) as parser:
            variants = list(parser)
        assert len(variants) == 3
        assert variants[0].chrom == "chr1"


class TestVCFParserContextManager:
    """Tests for context manager protocol."""

    def test_context_manager_closes_file(self, tmp_path: Path) -> None:
        """File handle is closed after exiting context."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        parser = VCFParser(vcf_path)
        with parser:
            _ = list(parser)
        assert parser._closed is True

    def test_close_idempotent(self, tmp_path: Path) -> None:
        """Calling close multiple times is safe."""
        vcf_path = _write_vcf(tmp_path, MINIMAL_VCF_CONTENT)
        parser = VCFParser(vcf_path)
        parser.close()
        parser.close()
        assert parser._closed is True


class TestVCFParserErrorHandling:
    """Tests for error handling on malformed VCF data."""

    def test_missing_fileformat_raises_parse_error(self, tmp_path: Path) -> None:
        """Raises ParseError when ##fileformat is missing."""
        bad_content = """\
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50
"""
        vcf_path = _write_vcf(tmp_path, bad_content)
        with pytest.raises(ParseError, match="fileformat"):
            VCFParser(vcf_path)

    def test_empty_file_raises_parse_error(self, tmp_path: Path) -> None:
        """Raises ParseError for an empty file."""
        vcf_path = tmp_path / "empty.vcf"
        vcf_path.write_text("")
        with pytest.raises((ParseError, FileNotFoundError)):
            VCFParser(vcf_path)
