"""Unit tests for SVTriageConfig validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.structural.config import SVTriageConfig


class TestSVTriageConfigValidation:
    def test_valid_defaults(self, tmp_path: Path) -> None:
        config = SVTriageConfig(
            vcf_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.json",
        )
        assert config.min_sv_size == 50
        assert config.max_allele_frequency == 0.01
        assert config.output_format == "json"

    def test_min_sv_size_too_small(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="min_sv_size"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                min_sv_size=0,
            )

    def test_min_sv_size_too_large(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="min_sv_size"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                min_sv_size=20_000_000,
            )

    def test_max_sv_size_negative(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="max_sv_size"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                max_sv_size=-1,
            )

    def test_max_sv_size_less_than_min(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="max_sv_size"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                min_sv_size=100,
                max_sv_size=50,
            )

    def test_max_allele_frequency_out_of_range(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="max_allele_frequency"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                max_allele_frequency=1.5,
            )

    def test_reciprocal_overlap_out_of_range(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="reciprocal_overlap"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                reciprocal_overlap=-0.1,
            )

    def test_whole_gene_threshold_too_low(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="whole_gene_threshold"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                whole_gene_threshold=0.3,
            )

    def test_min_quality_negative(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="min_quality"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                min_quality=-1.0,
            )

    def test_batch_size_too_small(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="batch_size"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                batch_size=50,
            )

    def test_batch_size_too_large(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="batch_size"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                batch_size=200_000,
            )

    def test_invalid_output_format(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="output_format"):
            SVTriageConfig(
                vcf_path=tmp_path / "x.vcf",
                output_path=tmp_path / "out.json",
                output_format="xml",
            )

    def test_csv_format_accepted(self, tmp_path: Path) -> None:
        config = SVTriageConfig(
            vcf_path=tmp_path / "x.vcf",
            output_path=tmp_path / "out.csv",
            output_format="csv",
        )
        assert config.output_format == "csv"
