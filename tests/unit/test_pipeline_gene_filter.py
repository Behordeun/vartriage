"""Unit tests for pipeline integration with gene filtering."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vartriage.models.config import (
    GeneFilterConfig,
    PipelineConfig,
)
from vartriage.models.variant import (
    AnnotatedVariant,
    FunctionalConsequence,
    Variant,
)
from vartriage.pipeline import Pipeline


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    gene_name: str | None = "BRCA1",
) -> AnnotatedVariant:
    return AnnotatedVariant(
        variant=Variant(
            chrom=chrom,
            pos=pos,
            id=None,
            ref="A",
            alt="T",
            qual=30.0,
            filter_status="PASS",
        ),
        consequence=FunctionalConsequence.MISSENSE,
        allele_frequency=0.001,
        gene_name=gene_name,
    )


class TestPipelineGeneFilterApplied:
    """Pipeline applies GeneFilter when gene_filter config is set."""

    def test_gene_filter_applied_after_annotation(
        self, tmp_path: Path
    ) -> None:
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\n")

        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            gene_filter=GeneFilterConfig(gene_list_path=gene_file),
        )

        brca1_variant = _make_variant(gene_name="BRCA1")
        egfr_variant = _make_variant(gene_name="EGFR")
        tp53_variant = _make_variant(gene_name="TP53")

        annotated_stream = iter(
            [brca1_variant, egfr_variant, tp53_variant]
        )

        with (
            patch.object(Pipeline, "_validate_config"),
            patch(
                "vartriage.pipeline.VCFParser"
            ) as mock_parser_cls,
            patch(
                "vartriage.pipeline.QualityFilter"
            ) as mock_qf_cls,
            patch(
                "vartriage.pipeline.AnnotationEngine"
            ),
            patch(
                "vartriage.pipeline.PrioritizationEngine"
            ) as mock_pri_cls,
            patch(
                "vartriage.pipeline.ACMGClassifier"
            ) as mock_acmg_cls,
            patch(
                "vartriage.pipeline.ReportGenerator"
            ) as mock_report_cls,
        ):
            pipeline = Pipeline(config)

            mock_parser = MagicMock()
            mock_parser.__enter__ = MagicMock(
                return_value=mock_parser
            )
            mock_parser.__exit__ = MagicMock(return_value=False)
            mock_parser.__iter__ = MagicMock(return_value=iter([]))
            mock_parser_cls.return_value = mock_parser

            mock_qf_cls.return_value.apply.return_value = iter([])

            pipeline._config = config

            mock_pri = mock_pri_cls.return_value
            captured_input: list[AnnotatedVariant] = []

            def capture_prioritize(variants):
                captured_input.extend(variants)
                return iter([])

            mock_pri.prioritize.side_effect = capture_prioritize

            mock_acmg_cls.return_value.classify.return_value = (
                iter([])
            )
            mock_report_cls.return_value.generate.return_value = (
                output_path
            )

            with patch.object(
                pipeline,
                "_passthrough_annotation",
                return_value=annotated_stream,
            ):
                pipeline.run()

            assert len(captured_input) == 2
            genes_passed = [
                v.gene_name for v in captured_input
            ]
            assert "BRCA1" in genes_passed
            assert "TP53" in genes_passed
            assert "EGFR" not in genes_passed


class TestPipelineNoGeneFilter:
    """Pipeline passes stream through unchanged without gene_filter."""

    def test_no_gene_filter_passes_all_variants(
        self, tmp_path: Path
    ) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            gene_filter=None,
        )

        brca1_variant = _make_variant(gene_name="BRCA1")
        egfr_variant = _make_variant(gene_name="EGFR")
        tp53_variant = _make_variant(gene_name="TP53")

        annotated_stream = iter(
            [brca1_variant, egfr_variant, tp53_variant]
        )

        with (
            patch.object(Pipeline, "_validate_config"),
            patch(
                "vartriage.pipeline.VCFParser"
            ) as mock_parser_cls,
            patch(
                "vartriage.pipeline.QualityFilter"
            ) as mock_qf_cls,
            patch(
                "vartriage.pipeline.AnnotationEngine"
            ),
            patch(
                "vartriage.pipeline.PrioritizationEngine"
            ) as mock_pri_cls,
            patch(
                "vartriage.pipeline.ACMGClassifier"
            ) as mock_acmg_cls,
            patch(
                "vartriage.pipeline.ReportGenerator"
            ) as mock_report_cls,
        ):
            pipeline = Pipeline(config)

            mock_parser = MagicMock()
            mock_parser.__enter__ = MagicMock(
                return_value=mock_parser
            )
            mock_parser.__exit__ = MagicMock(return_value=False)
            mock_parser.__iter__ = MagicMock(return_value=iter([]))
            mock_parser_cls.return_value = mock_parser

            mock_qf_cls.return_value.apply.return_value = iter([])

            mock_pri = mock_pri_cls.return_value
            captured_input: list[AnnotatedVariant] = []

            def capture_prioritize(variants):
                captured_input.extend(variants)
                return iter([])

            mock_pri.prioritize.side_effect = capture_prioritize

            mock_acmg_cls.return_value.classify.return_value = (
                iter([])
            )
            mock_report_cls.return_value.generate.return_value = (
                output_path
            )

            with patch.object(
                pipeline,
                "_passthrough_annotation",
                return_value=annotated_stream,
            ):
                pipeline.run()

            assert len(captured_input) == 3
            genes_passed = [
                v.gene_name for v in captured_input
            ]
            assert "BRCA1" in genes_passed
            assert "EGFR" in genes_passed
            assert "TP53" in genes_passed


class TestPipelineConfigValidation:
    """Pipeline config validation for gene_filter path."""

    def test_raises_file_not_found_for_missing_gene_list(
        self, tmp_path: Path
    ) -> None:
        missing_path = tmp_path / "nonexistent_genes.txt"
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            gene_filter=GeneFilterConfig(
                gene_list_path=missing_path
            ),
        )

        with pytest.raises(FileNotFoundError, match="Gene list"):
            Pipeline(config)
