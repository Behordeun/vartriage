"""Unit tests for CLI --gene-list argument integration."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch, MagicMock

from vartriage.cli import _build_parser, _run_pipeline
from vartriage.models.config import GeneFilterConfig


class TestGeneListArgumentParsing:
    """Verify --gene-list is registered and parses correctly."""

    def test_gene_list_argument_accepted(self) -> None:
        parser = _build_parser()
        args = parser.parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
            "--gene-list", "/tmp/genes.txt",
        ])
        assert args.gene_list == Path("/tmp/genes.txt")

    def test_gene_list_defaults_to_none(self) -> None:
        parser = _build_parser()
        args = parser.parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
        ])
        assert args.gene_list is None

    def test_gene_list_is_path_type(self) -> None:
        parser = _build_parser()
        args = parser.parse_args([
            "--vcf", "input.vcf",
            "--output", "out.json",
            "--gene-list", "my_panel.txt",
        ])
        assert isinstance(args.gene_list, Path)


class TestGeneListConfigConstruction:
    """Verify CLI wires --gene-list into PipelineConfig correctly."""

    @patch("vartriage.pipeline.Pipeline")
    def test_gene_list_produces_gene_filter_config(
        self, mock_pipeline_cls: MagicMock, tmp_path: Path
    ) -> None:
        vcf_file = tmp_path / "input.vcf"
        vcf_file.touch()
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\n")

        mock_pipeline = MagicMock()
        mock_pipeline.run.return_value = tmp_path / "out.json"
        mock_pipeline_cls.return_value = mock_pipeline

        parser = _build_parser()
        args = parser.parse_args([
            "--vcf", str(vcf_file),
            "--output", str(tmp_path / "out.json"),
            "--gene-list", str(gene_file),
        ])

        _run_pipeline(args, vcf_file)

        config = mock_pipeline_cls.call_args[0][0]
        assert config.gene_filter is not None
        assert isinstance(config.gene_filter, GeneFilterConfig)
        assert config.gene_filter.gene_list_path == gene_file

    @patch("vartriage.pipeline.Pipeline")
    def test_no_gene_list_sets_gene_filter_to_none(
        self, mock_pipeline_cls: MagicMock, tmp_path: Path
    ) -> None:
        vcf_file = tmp_path / "input.vcf"
        vcf_file.touch()

        mock_pipeline = MagicMock()
        mock_pipeline.run.return_value = tmp_path / "out.json"
        mock_pipeline_cls.return_value = mock_pipeline

        parser = _build_parser()
        args = parser.parse_args([
            "--vcf", str(vcf_file),
            "--output", str(tmp_path / "out.json"),
        ])

        _run_pipeline(args, vcf_file)

        config = mock_pipeline_cls.call_args[0][0]
        assert config.gene_filter is None
