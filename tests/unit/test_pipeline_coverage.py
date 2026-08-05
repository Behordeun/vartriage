"""Unit tests covering pipeline.py logic branches not hit by existing tests."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vartriage.models.config import (
    AnnotationConfig,
    InheritanceConfig,
    PipelineConfig,
    PrioritizationConfig,
    RegionFilterConfig,
)
from vartriage.models.variant import (
    AnnotatedVariant,
    ClassifiedVariant,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.pipeline import Pipeline


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    gene_name: str | None = "BRCA1",
    info: dict | None = None,
) -> AnnotatedVariant:
    return AnnotatedVariant(
        variant=Variant(
            chrom=chrom,
            pos=pos,
            id=None,
            ref=ref,
            alt=alt,
            qual=30.0,
            filter_status="PASS",
            info=info or {},
        ),
        consequence=FunctionalConsequence.MISSENSE,
        allele_frequency=0.001,
        gene_name=gene_name,
    )


def _make_classified(scored: ScoredVariant) -> ClassifiedVariant:
    return ClassifiedVariant(scored=scored)


def _pipeline_with_mocked_validation(config: PipelineConfig) -> Pipeline:
    with patch.object(Pipeline, "_validate_config"):
        return Pipeline(config)


class TestPassthroughAnnotation:
    """_passthrough_annotation wraps raw Variants as AnnotatedVariant."""

    def test_wraps_variants_with_intergenic_consequence(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(vcf_path=vcf_path, output_path=output_path)

        pipeline = _pipeline_with_mocked_validation(config)

        raw = Variant(
            chrom="chr7",
            pos=500,
            id=None,
            ref="G",
            alt="C",
            qual=40.0,
            filter_status="PASS",
        )
        results = list(pipeline._passthrough_annotation(iter([raw])))

        assert len(results) == 1
        av = results[0]
        assert av.consequence == FunctionalConsequence.INTERGENIC
        assert av.allele_frequency is None
        assert av.clinvar_assertion is None
        assert av.frequency_unknown is True
        assert av.clinvar_unknown is True
        assert av.variant.chrom == "chr7"

    def test_passthrough_preserves_variant_info(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(vcf_path=vcf_path, output_path=output_path)
        pipeline = _pipeline_with_mocked_validation(config)

        raw = Variant(
            chrom="chr1",
            pos=10,
            id="rs123",
            ref="T",
            alt="A",
            qual=99.0,
            filter_status="PASS",
            info={"DP": 50},
        )
        results = list(pipeline._passthrough_annotation(iter([raw])))
        assert results[0].variant.info == {"DP": 50}


class TestVariantsWithGeneInfo:
    """_variants_with_gene_info copies gene_name into variant info dict."""

    def test_copies_gene_name_into_info(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(vcf_path=vcf_path, output_path=output_path)
        pipeline = _pipeline_with_mocked_validation(config)

        av = _make_variant(gene_name="TP53")
        results = list(pipeline._variants_with_gene_info(iter([av])))

        assert len(results) == 1
        assert results[0].info["gene"] == "TP53"
        assert results[0].chrom == "chr1"

    def test_skips_gene_key_when_none(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(vcf_path=vcf_path, output_path=output_path)
        pipeline = _pipeline_with_mocked_validation(config)

        av = _make_variant(gene_name=None)
        results = list(pipeline._variants_with_gene_info(iter([av])))

        assert "gene" not in results[0].info


class TestReattachAnnotations:
    """_reattach_annotations maps inherited variants back to annotation data."""

    def test_reattaches_matching_annotations(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(vcf_path=vcf_path, output_path=output_path)
        pipeline = _pipeline_with_mocked_validation(config)

        av = _make_variant(chrom="chr2", pos=200, ref="G", alt="A", gene_name="EGFR")
        inherited_variant = Variant(
            chrom="chr2",
            pos=200,
            id=None,
            ref="G",
            alt="A",
            qual=30.0,
            filter_status="PASS",
            info={"inheritance_pattern": "de_novo"},
        )

        results = pipeline._reattach_annotations([inherited_variant], [av])

        assert len(results) == 1
        assert results[0].gene_name == "EGFR"
        assert results[0].variant.info["inheritance_pattern"] == "de_novo"

    def test_falls_back_to_intergenic_when_no_match(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(vcf_path=vcf_path, output_path=output_path)
        pipeline = _pipeline_with_mocked_validation(config)

        # No annotations match chr3:300
        av = _make_variant(chrom="chr2", pos=200, ref="G", alt="A")
        orphan = Variant(
            chrom="chr3",
            pos=300,
            id=None,
            ref="C",
            alt="T",
            qual=20.0,
            filter_status="PASS",
        )

        results = pipeline._reattach_annotations([orphan], [av])
        assert results[0].consequence == FunctionalConsequence.INTERGENIC
        assert results[0].frequency_unknown is True
        assert results[0].clinvar_unknown is True


class TestCollectReferencePaths:
    """_collect_reference_paths gathers paths from annotation + prioritization configs."""

    def test_collects_annotation_and_prioritization_paths(self, tmp_path: Path) -> None:
        gene_file = tmp_path / "genes.gtf"
        gene_file.write_text("")
        gnomad_file = tmp_path / "gnomad.vcf.gz"
        gnomad_file.write_text("")
        clinvar_file = tmp_path / "clinvar.vcf.gz"
        clinvar_file.write_text("")
        cadd_file = tmp_path / "cadd.tsv"
        cadd_file.write_text("")

        config = PipelineConfig(
            vcf_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.json",
            annotation=AnnotationConfig(
                gene_annotation_path=gene_file,
                gnomad_path=gnomad_file,
                clinvar_path=clinvar_file,
            ),
            prioritization=PrioritizationConfig(cadd_scores_path=cadd_file),
        )
        with patch.object(Pipeline, "_validate_config"):
            pipeline = Pipeline(config)

        paths = pipeline._collect_reference_paths()
        assert gene_file in paths
        assert gnomad_file in paths
        assert clinvar_file in paths
        assert cadd_file in paths

    def test_skips_none_clinvar_path(self, tmp_path: Path) -> None:
        gene_file = tmp_path / "genes.gtf"
        gene_file.write_text("")
        gnomad_file = tmp_path / "gnomad.vcf.gz"
        gnomad_file.write_text("")

        config = PipelineConfig(
            vcf_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.json",
            annotation=AnnotationConfig(
                gene_annotation_path=gene_file,
                gnomad_path=gnomad_file,
                clinvar_path=None,
            ),
            prioritization=PrioritizationConfig(),
        )
        with patch.object(Pipeline, "_validate_config"):
            pipeline = Pipeline(config)

        paths = pipeline._collect_reference_paths()
        assert gene_file in paths
        assert gnomad_file in paths
        # clinvar_path=None means it's excluded from the list
        assert len(paths) == 2


class TestValidateConfig:
    """_validate_config fail-fast paths."""

    def test_rejects_missing_annotation_gene_file(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        gnomad_file = tmp_path / "gnomad.vcf.gz"
        gnomad_file.write_text("")
        missing_gene = tmp_path / "missing_genes.gtf"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=tmp_path / "output.json",
            annotation=AnnotationConfig(
                gene_annotation_path=missing_gene,
                gnomad_path=gnomad_file,
            ),
        )
        with pytest.raises(FileNotFoundError, match="Gene annotation"):
            Pipeline(config)

    def test_rejects_missing_gnomad_file(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        gene_file = tmp_path / "genes.gtf"
        gene_file.write_text("")
        missing_gnomad = tmp_path / "missing_gnomad.vcf.gz"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=tmp_path / "output.json",
            annotation=AnnotationConfig(
                gene_annotation_path=gene_file,
                gnomad_path=missing_gnomad,
            ),
        )
        with pytest.raises(FileNotFoundError, match="gnomAD"):
            Pipeline(config)

    def test_rejects_missing_cadd_scores_file(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=tmp_path / "output.json",
            prioritization=PrioritizationConfig(
                cadd_scores_path=tmp_path / "missing_cadd.tsv"
            ),
        )
        with pytest.raises(FileNotFoundError, match="CADD"):
            Pipeline(config)

    def test_rejects_compound_het_without_annotation(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=tmp_path / "output.json",
            inheritance=InheritanceConfig(
                proband="child",
                mother="mom",
                father="dad",
                patterns=["compound_het"],
            ),
            annotation=None,
        )
        with pytest.raises(ValueError, match="compound_het"):
            Pipeline(config)

    def test_rejects_missing_region_bed_file(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=tmp_path / "output.json",
            region_filter=RegionFilterConfig(bed_path=tmp_path / "missing.bed"),
        )
        with pytest.raises(FileNotFoundError, match="BED"):
            Pipeline(config)


class TestGenerateReport:
    """_generate_report dispatches VCF format with source path."""

    def test_vcf_format_passes_source_path(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.vcf"
        from vartriage.models.config import ReportConfig

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="vcf"),
        )
        pipeline = _pipeline_with_mocked_validation(config)

        mock_gen = MagicMock()
        mock_gen.generate.return_value = output_path

        result = pipeline._generate_report(mock_gen, iter([]), output_path, vcf_path)
        mock_gen.generate.assert_called_once_with(
            mock_gen.generate.call_args[0][0], output_path, vcf_path
        )
        assert result == output_path

    def test_non_vcf_format_omits_source_path(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        from vartriage.models.config import ReportConfig

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            report=ReportConfig(output_format="json"),
        )
        pipeline = _pipeline_with_mocked_validation(config)

        mock_gen = MagicMock()
        mock_gen.generate.return_value = output_path

        result = pipeline._generate_report(mock_gen, iter([]), output_path, vcf_path)
        # Non-VCF format: no source_vcf_path argument
        call_args = mock_gen.generate.call_args
        assert len(call_args[0]) == 2
        assert result == output_path


class TestRunWithSv:
    """run_with_sv triggers SV pipeline when sv_vcf_path is configured."""

    def test_skips_sv_when_path_is_none(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            sv_vcf_path=None,
        )
        pipeline = _pipeline_with_mocked_validation(config)

        with patch.object(pipeline, "run") as mock_run:
            mock_run.return_value = output_path
            pipeline.run_with_sv()
            mock_run.assert_called_once()

    def test_skips_sv_when_path_does_not_exist(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        sv_path = tmp_path / "nonexistent_sv.vcf"

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            sv_vcf_path=sv_path,
        )
        pipeline = _pipeline_with_mocked_validation(config)

        with patch.object(pipeline, "run") as mock_run:
            mock_run.return_value = output_path
            pipeline.run_with_sv()
            mock_run.assert_called_once()

    def test_runs_sv_pipeline_when_path_exists(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text("")
        output_path = tmp_path / "output.json"
        sv_path = tmp_path / "sv.vcf"
        sv_path.write_text("")

        config = PipelineConfig(
            vcf_path=vcf_path,
            output_path=output_path,
            sv_vcf_path=sv_path,
        )
        pipeline = _pipeline_with_mocked_validation(config)

        with (
            patch.object(pipeline, "run") as mock_run,
            patch("vartriage.structural.pipeline.SVTriagePipeline") as mock_sv_cls,
        ):
            mock_run.return_value = output_path
            mock_sv_pipeline = MagicMock()
            mock_sv_cls.return_value = mock_sv_pipeline

            pipeline.run_with_sv()

            mock_sv_pipeline.run.assert_called_once()
