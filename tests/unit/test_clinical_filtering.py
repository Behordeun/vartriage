"""Unit tests for RegionFilter, SampleExtractor, and CLI integration."""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path
from unittest.mock import patch

import pytest

from vartriage.filter.region_filter import RegionFilter
from vartriage.filter.sample_extractor import SampleExtractor
from vartriage.io.exceptions import ParseError
from vartriage.models.config import (
    RegionFilterConfig,
    SampleConfig,
)
from vartriage.models.variant import Variant


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    info: dict | None = None,
) -> Variant:
    """Create a minimal test variant."""
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
        info=info or {},
    )


def _write_bed(content: str) -> Path:
    """Write content to a temp BED file, return path."""
    f = tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    )
    f.write(content)
    f.close()
    return Path(f.name)


def _make_sample_variant(
    sample_name: str,
    gt: tuple,
    gq: int | None = 50,
    chrom: str = "chr1",
    pos: int = 100,
) -> Variant:
    """Create a variant with _pysam_samples data attached."""
    sample_data: dict = {"GT": gt}
    if gq is not None:
        sample_data["GQ"] = gq
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
        info={"_pysam_samples": {sample_name: sample_data}},
    )


# ===========================================================================
# RegionFilter Tests
# ===========================================================================


class TestRegionFilterBEDLoading:
    """BED file loading behavior."""

    def test_four_column_bed(self) -> None:
        """4-column BED (chrom, start, end, name) loads fine."""
        content = "chr1\t100\t200\tgene_A\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            assert "chr1" in rf._intervals
            assert len(rf._intervals["chr1"]) == 1
            assert rf._intervals["chr1"][0] == (100, 200)
        finally:
            bed_path.unlink(missing_ok=True)

    def test_empty_bed_file(self) -> None:
        """Empty BED file gives 0 intervals, filters everything."""
        bed_path = _write_bed("")
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            assert rf._intervals == {}
            v = _make_variant()
            result = list(rf.apply(iter([v])))
            assert result == []
        finally:
            bed_path.unlink(missing_ok=True)

    def test_overlap_at_start_boundary(self) -> None:
        """VCF pos at start+1 (0-based start) overlaps."""
        content = "chr1\t100\t200\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            v = _make_variant(pos=101)
            result = list(rf.apply(iter([v])))
            assert result == [v]
        finally:
            bed_path.unlink(missing_ok=True)

    def test_overlap_at_end_minus_one(self) -> None:
        """VCF pos at end (0-based end-1) overlaps [start, end)."""
        content = "chr1\t100\t200\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            v = _make_variant(pos=200)
            result = list(rf.apply(iter([v])))
            assert result == [v]
        finally:
            bed_path.unlink(missing_ok=True)

    def test_no_overlap_one_past_end(self) -> None:
        """VCF pos at end+1 (0-based = end) does NOT overlap."""
        content = "chr1\t100\t200\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            v = _make_variant(pos=201)
            result = list(rf.apply(iter([v])))
            assert result == []
        finally:
            bed_path.unlink(missing_ok=True)

    def test_missing_file_raises_file_not_found(self) -> None:
        """Non-existent BED path raises FileNotFoundError."""
        fake_path = Path("/tmp/nonexistent_bed_file_xyz.bed")
        config = RegionFilterConfig(bed_path=fake_path)
        with pytest.raises(FileNotFoundError):
            RegionFilter(config)

    def test_malformed_line_raises_parse_error_with_line_number(
        self,
    ) -> None:
        """Malformed BED line raises ParseError with correct line number."""
        content = "chr1\t100\t200\nchr2\tnot_a_number\t300\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            with pytest.raises(ParseError) as exc_info:
                RegionFilter(config)
            assert exc_info.value.line_number == 2
        finally:
            bed_path.unlink(missing_ok=True)

    def test_too_few_columns_raises_parse_error(self) -> None:
        """Line with < 3 columns raises ParseError."""
        content = "chr1\t100\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            with pytest.raises(ParseError) as exc_info:
                RegionFilter(config)
            assert exc_info.value.line_number == 1
        finally:
            bed_path.unlink(missing_ok=True)

    def test_comment_lines_skipped(self) -> None:
        """Lines starting with # are skipped."""
        content = "# header comment\nchr1\t100\t200\n"
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            assert len(rf._intervals.get("chr1", [])) == 1
        finally:
            bed_path.unlink(missing_ok=True)

    def test_track_and_browser_lines_skipped(self) -> None:
        """Track and browser header lines are skipped."""
        content = (
            "browser position chr1:1-1000\n"
            "track name=test\n"
            "chr1\t50\t150\n"
        )
        bed_path = _write_bed(content)
        try:
            config = RegionFilterConfig(bed_path=bed_path)
            rf = RegionFilter(config)
            assert len(rf._intervals.get("chr1", [])) == 1
        finally:
            bed_path.unlink(missing_ok=True)


# ===========================================================================
# SampleExtractor Tests
# ===========================================================================


class TestSampleExtractorGenotype:
    """Genotype-based filtering."""

    def test_het_included(self) -> None:
        """Het (0/1) passes."""
        sample = "S1"
        config = SampleConfig(sample_name=sample)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1))
        result = list(ext.apply(iter([v])))
        assert len(result) == 1

    def test_hom_ref_excluded(self) -> None:
        """Hom ref (0/0) excluded."""
        sample = "S1"
        config = SampleConfig(sample_name=sample)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 0))
        result = list(ext.apply(iter([v])))
        assert len(result) == 0

    def test_no_call_excluded(self) -> None:
        """No-call (./.) excluded."""
        sample = "S1"
        config = SampleConfig(sample_name=sample)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(None, None))
        result = list(ext.apply(iter([v])))
        assert len(result) == 0

    def test_hom_alt_included(self) -> None:
        """Hom alt (1/1) passes."""
        sample = "S1"
        config = SampleConfig(sample_name=sample)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(1, 1))
        result = list(ext.apply(iter([v])))
        assert len(result) == 1


class TestSampleExtractorGQ:
    """Genotype quality filtering."""

    def test_gq_below_threshold_excluded(self) -> None:
        """GQ below min_gq excludes."""
        sample = "S1"
        config = SampleConfig(sample_name=sample, min_gq=30)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1), gq=20)
        result = list(ext.apply(iter([v])))
        assert len(result) == 0

    def test_gq_above_threshold_included(self) -> None:
        """GQ above min_gq passes."""
        sample = "S1"
        config = SampleConfig(sample_name=sample, min_gq=30)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1), gq=50)
        result = list(ext.apply(iter([v])))
        assert len(result) == 1

    def test_gq_at_threshold_included(self) -> None:
        """GQ exactly at min_gq passes (>= comparison)."""
        sample = "S1"
        config = SampleConfig(sample_name=sample, min_gq=30)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1), gq=30)
        result = list(ext.apply(iter([v])))
        assert len(result) == 1

    def test_missing_gq_emits_warning_and_excludes(self) -> None:
        """Missing GQ with min_gq set emits UserWarning, excludes."""
        sample = "S1"
        config = SampleConfig(sample_name=sample, min_gq=30)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1), gq=None)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = list(ext.apply(iter([v])))

        assert len(result) == 0, "Missing GQ should exclude"
        assert len(w) >= 1, "Should emit a warning"
        assert issubclass(w[0].category, UserWarning)


class TestSampleExtractorValidation:
    """Sample validation and output structure."""

    def test_sample_not_in_header_raises_value_error(self) -> None:
        """Unknown sample raises ValueError."""
        config = SampleConfig(sample_name="MISSING")
        with pytest.raises(ValueError, match="MISSING"):
            SampleExtractor(config, ["S1", "S2", "S3"])

    def test_available_samples_in_error_message(self) -> None:
        """ValueError lists available sample names."""
        config = SampleConfig(sample_name="MISSING")
        with pytest.raises(ValueError, match="S1.*S2.*S3|Available"):
            SampleExtractor(config, ["S1", "S2", "S3"])

    def test_output_has_sample_gt_gq_name(self) -> None:
        """Output info contains sample_gt, sample_gq, sample_name."""
        sample = "PROBAND"
        config = SampleConfig(sample_name=sample)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1), gq=45)
        result = list(ext.apply(iter([v])))

        assert len(result) == 1
        info = result[0].info
        assert "sample_gt" in info
        assert info["sample_gt"] == "0/1"
        assert "sample_name" in info
        assert info["sample_name"] == sample
        assert "sample_gq" in info
        assert info["sample_gq"] == 45

    def test_pysam_samples_stripped_from_output(self) -> None:
        """_pysam_samples removed from output info."""
        sample = "S1"
        config = SampleConfig(sample_name=sample)
        ext = SampleExtractor(config, [sample])
        v = _make_sample_variant(sample, gt=(0, 1), gq=50)
        result = list(ext.apply(iter([v])))

        assert len(result) == 1
        assert "_pysam_samples" not in result[0].info


# ===========================================================================
# CLI Tests
# ===========================================================================


class TestCLIParsing:
    """CLI argument parsing for clinical filtering options."""

    def _parse_args(self, argv: list[str]) -> object:
        """Parse CLI args using the internal parser builder."""
        from vartriage.cli import _build_parser

        parser = _build_parser()
        return parser.parse_args(argv)

    def test_regions_parsed(self) -> None:
        """--regions parsed as Path."""
        args = self._parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
            "--regions", "/path/to/regions.bed",
        ])
        assert args.regions == Path("/path/to/regions.bed")

    def test_sample_parsed(self) -> None:
        """--sample parsed as string."""
        args = self._parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
            "--sample", "NA12878",
        ])
        assert args.sample == "NA12878"

    def test_min_gq_without_sample_prints_error(
        self, capsys: pytest.CaptureFixture
    ) -> None:
        """--min-gq without --sample triggers error exit."""
        import sys

        args = self._parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
            "--min-gq", "30",
        ])
        assert args.min_gq == 30
        assert args.sample is None

        from unittest.mock import MagicMock

        with pytest.raises(SystemExit) as exc_info:
            from vartriage.cli import _run_pipeline

            with tempfile.NamedTemporaryFile(
                suffix=".vcf", delete=False
            ) as f:
                f.write(b"##fileformat=VCFv4.2\n")
                vcf_path = Path(f.name)

            try:
                _run_pipeline(args, vcf_path)
            finally:
                vcf_path.unlink(missing_ok=True)

        assert exc_info.value.code == 2
        captured = capsys.readouterr()
        assert (
            "min-gq" in captured.err.lower()
            or "sample" in captured.err.lower()
        )

    def test_all_absent_produces_none(self) -> None:
        """No --regions/--sample/--min-gq means all None."""
        args = self._parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
        ])
        assert args.regions is None
        assert args.sample is None
        assert args.min_gq is None

    def test_min_gq_parsed_as_int(self) -> None:
        """--min-gq parsed as integer."""
        args = self._parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
            "--sample", "S1",
            "--min-gq", "25",
        ])
        assert args.min_gq == 25
        assert isinstance(args.min_gq, int)
