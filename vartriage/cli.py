"""CLI entry point for the vartriage pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal, Optional, cast

from vartriage.models.config import ClinicalReportConfig, InheritanceConfig


def _get_version() -> str:
    """Return the installed package version, falling back to __version__."""
    try:
        from importlib.metadata import version

        return version("vartriage")
    except Exception:
        from vartriage import __version__

        return __version__


def _build_parser() -> argparse.ArgumentParser:
    """Build the argparse parser with all CLI options."""
    parser = argparse.ArgumentParser(
        prog="vartriage",
        description=(
            "Variant prioritization and ACMG classification pipeline. "
            "Reads a VCF file, applies quality filtering, annotation, "
            "prioritization, and ACMG classification, then writes a "
            "structured report."
        ),
    )

    parser.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="Path to the input VCF file (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path where the output report will be written",
    )
    parser.add_argument(
        "--output-format",
        choices=[
            "json",
            "csv",
            "pdf",
            "vcf",
            "clinical-pdf",
            "clinical-html",
            "clinical-docx",
        ],
        default="json",
        help="Output report format (default: json)",
    )
    parser.add_argument(
        "--gene-annotation",
        type=Path,
        default=None,
        help="Path to GTF/GFF gene annotation reference file",
    )
    parser.add_argument(
        "--gnomad",
        type=Path,
        default=None,
        help="Path to gnomAD population frequency reference file",
    )
    parser.add_argument(
        "--clinvar",
        type=Path,
        default=None,
        help="Path to ClinVar clinical significance reference file",
    )
    parser.add_argument(
        "--cadd-scores",
        type=Path,
        default=None,
        help="Path to CADD Phred score TSV reference file",
    )
    parser.add_argument(
        "--revel-scores",
        type=Path,
        default=None,
        help="Path to REVEL score TSV reference file",
    )
    parser.add_argument(
        "--spliceai-scores",
        type=Path,
        default=None,
        help="Path to SpliceAI score TSV reference file",
    )
    parser.add_argument(
        "--gene-list",
        type=Path,
        default=None,
        help="Path to a gene list file for gene-based filtering",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {_get_version()}",
    )

    parser.add_argument(
        "--regions",
        type=Path,
        default=None,
        help="Path to BED file for region-based variant filtering",
    )
    parser.add_argument(
        "--sample",
        type=str,
        default=None,
        help="Sample name to extract from multi-sample VCF",
    )
    parser.add_argument(
        "--min-gq",
        type=int,
        default=None,
        help="Minimum genotype quality threshold (0-99). Requires --sample",
    )
    parser.add_argument(
        "--proband",
        type=str,
        default=None,
        help="Proband sample name for trio inheritance analysis",
    )
    parser.add_argument(
        "--mother",
        type=str,
        default=None,
        help="Mother sample name for trio inheritance analysis",
    )
    parser.add_argument(
        "--father",
        type=str,
        default=None,
        help="Father sample name for trio inheritance analysis",
    )
    parser.add_argument(
        "--inheritance-pattern",
        type=str,
        action="append",
        default=None,
        help=(
            "Inheritance pattern to evaluate (may be specified "
            "multiple times). Supported: de_novo, dominant, "
            "recessive, compound_het, x_linked. Defaults to all "
            "patterns when trio is active."
        ),
    )
    parser.add_argument(
        "--patient-id",
        type=str,
        default=None,
        help="Patient identifier for clinical reports",
    )
    parser.add_argument(
        "--panel-name",
        type=str,
        default=None,
        help="Gene panel name for clinical reports",
    )

    return parser


def main(argv: Optional[list[str]] = None) -> None:
    """Run the vartriage CLI.

    Parameters
    ----------
    argv : list[str], optional
        Arguments to parse. Uses sys.argv[1:] when None.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Validate clinical format requirements before VCF check.
    # _build_clinical_config is the single source of truth for
    # this validation and will sys.exit(2) on missing flags.
    output_fmt: str = args.output_format
    clinical_config = _build_clinical_config(args, output_fmt)

    vcf_path: Path = args.vcf
    if not vcf_path.exists():
        print(
            f"Error: VCF file not found: {vcf_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        result_path = _run_pipeline(args, vcf_path, clinical_config)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    except OSError as exc:
        print(
            f"Error: report generation failed: {exc}",
            file=sys.stderr,
        )
        sys.exit(1)
    except Exception as exc:
        _handle_unexpected_error(exc)

    print(str(result_path))
    sys.exit(0)


def _handle_unexpected_error(exc: Exception) -> None:
    """Print an appropriate error message and exit."""
    from vartriage.io.exceptions import VariantPrioritizationError

    if isinstance(exc, VariantPrioritizationError):
        print(f"Error: pipeline failed: {exc}", file=sys.stderr)
    else:
        print(f"Error: unexpected failure: {exc}", file=sys.stderr)
    sys.exit(1)


def _run_pipeline(
    args: argparse.Namespace,
    vcf_path: Path,
    clinical_config: Optional[ClinicalReportConfig] = None,
) -> Path:
    """Assemble pipeline config from parsed args and run it.

    Returns
    -------
    Path
        Path to the generated report.
    """
    from vartriage.models.config import (AnnotationConfig, GeneFilterConfig,
                                         PipelineConfig, PrioritizationConfig,
                                         RegionFilterConfig, ReportConfig,
                                         SampleConfig)
    from vartriage.pipeline import Pipeline

    output_format: str = args.output_format

    inheritance_config = _build_inheritance_config(args)

    annotation_config: Optional[AnnotationConfig] = None
    gene_annotation: Optional[Path] = args.gene_annotation
    gnomad: Optional[Path] = args.gnomad

    if gene_annotation is not None and gnomad is not None:
        annotation_config = AnnotationConfig(
            gene_annotation_path=gene_annotation,
            gnomad_path=gnomad,
            clinvar_path=args.clinvar,
        )

    prioritization_config = PrioritizationConfig(
        cadd_scores_path=args.cadd_scores,
        revel_scores_path=args.revel_scores,
        spliceai_scores_path=args.spliceai_scores,
    )

    report_fmt = cast(
        Literal[
            "json",
            "csv",
            "pdf",
            "vcf",
            "clinical-pdf",
            "clinical-html",
            "clinical-docx",
        ],
        output_format,
    )
    report_config = ReportConfig(
        output_format=report_fmt,
    )

    gene_filter_config: Optional[GeneFilterConfig] = None
    if args.gene_list is not None:
        gene_filter_config = GeneFilterConfig(
            gene_list_path=args.gene_list,
        )

    region_filter_config: Optional[RegionFilterConfig] = None
    if args.regions is not None:
        region_filter_config = RegionFilterConfig(
            bed_path=args.regions,
        )

    sample_config: Optional[SampleConfig] = None
    if args.sample is not None:
        sample_config = SampleConfig(
            sample_name=args.sample,
            min_gq=args.min_gq,
        )

    if args.min_gq is not None and args.sample is None:
        print(
            "Error: --min-gq requires --sample",
            file=sys.stderr,
        )
        sys.exit(2)

    pipeline_config = PipelineConfig(
        vcf_path=vcf_path,
        output_path=args.output,
        annotation=annotation_config,
        prioritization=prioritization_config,
        report=report_config,
        inheritance=inheritance_config,
        gene_filter=gene_filter_config,
        region_filter=region_filter_config,
        sample=sample_config,
        clinical_report=clinical_config,
    )

    pipeline = Pipeline(pipeline_config)
    return pipeline.run()


def _build_clinical_config(
    args: argparse.Namespace,
    output_format: str,
) -> Optional[ClinicalReportConfig]:
    """Build ClinicalReportConfig if clinical format is requested."""
    if not output_format.startswith("clinical-"):
        return None

    from vartriage.models.config import ClinicalReportConfig as _CRC

    """Build ClinicalReportConfig if clinical format is requested."""
    if not output_format.startswith("clinical-"):
        return None

    from vartriage.models.config import ClinicalReportConfig as _CRC

    clinical_fmt = cast(
        Literal["clinical-pdf", "clinical-html", "clinical-docx"],
        output_format,
    )
    return _CRC(
        patient_id=args.patient_id,
        panel_name=args.panel_name,
        output_format=clinical_fmt,
    )


def _build_inheritance_config(
    args: argparse.Namespace,
) -> Optional[InheritanceConfig]:
    """Build InheritanceConfig from trio arguments."""
    from vartriage.models.config import InheritanceConfig as _IC

    proband: Optional[str] = args.proband
    mother: Optional[str] = args.mother
    father: Optional[str] = args.father
    sample: Optional[str] = getattr(args, "sample", None)

    trio_args = [proband, mother, father]
    trio_provided = [a for a in trio_args if a is not None]

    if trio_provided and len(trio_provided) < 3:
        print(
            "Error: --proband, --mother, and --father must all "
            "be provided together for trio analysis.",
            file=sys.stderr,
        )
        sys.exit(1)

    if trio_provided and sample is not None:
        print(
            "Error: --sample and trio arguments (--proband, "
            "--mother, --father) are mutually exclusive.",
            file=sys.stderr,
        )
        sys.exit(1)

    if len(trio_provided) != 3:
        return None

    patterns = args.inheritance_pattern
    if patterns is None:
        patterns = list(_IC.SUPPORTED_PATTERNS)
    return _IC(
        proband=proband,  # type: ignore[arg-type]
        mother=mother,  # type: ignore[arg-type]
        father=father,  # type: ignore[arg-type]
        patterns=patterns,
    )
