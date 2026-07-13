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
    parser.add_argument(
        "--use-bundles",
        action="store_true",
        default=False,
        help=(
            "Auto-resolve reference file paths from installed bundles "
            "(~/.vartriage/bundles/). Paths explicitly passed via "
            "--gnomad, --clinvar, etc. take precedence."
        ),
    )
    parser.add_argument(
        "--genome-build",
        type=str,
        default="grch38",
        help="Genome build for bundle resolution (default: grch38)",
    )

    return parser


def main(argv: Optional[list[str]] = None) -> None:
    """Run the vartriage CLI.

    Parameters
    ----------
    argv : list[str], optional
        Arguments to parse. Uses sys.argv[1:] when None.
    """
    # Intercept 'bundle' subcommand before main parser
    effective_argv = argv if argv is not None else sys.argv[1:]
    if effective_argv and effective_argv[0] == "bundle":
        _run_bundle_cli(effective_argv[1:])
        return

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

    # Bundle auto-resolution: fill in missing paths from installed bundles
    use_bundles: bool = getattr(args, "use_bundles", False)
    genome_build: str = getattr(args, "genome_build", "grch38")

    annotation_config: Optional[AnnotationConfig] = None
    gene_annotation: Optional[Path] = args.gene_annotation
    gnomad: Optional[Path] = args.gnomad
    clinvar: Optional[Path] = args.clinvar
    cadd_scores: Optional[Path] = args.cadd_scores
    revel_scores: Optional[Path] = args.revel_scores
    spliceai_scores: Optional[Path] = args.spliceai_scores

    if use_bundles:
        from vartriage.bundle.storage import BundleStorage
        storage = BundleStorage()

        if gene_annotation is None:
            resolved = storage.resolve_path(genome_build, "gencode")
            if resolved:
                gene_annotation = resolved

        if gnomad is None:
            resolved = storage.resolve_path(genome_build, "gnomad-exomes-chr22")
            if resolved:
                gnomad = resolved

        if clinvar is None:
            resolved = storage.resolve_path(genome_build, "clinvar")
            if resolved:
                clinvar = resolved

        if cadd_scores is None:
            resolved = storage.resolve_path(genome_build, "cadd")
            if resolved:
                cadd_scores = resolved

        if revel_scores is None:
            resolved = storage.resolve_path(genome_build, "revel")
            if resolved:
                revel_scores = resolved

        if spliceai_scores is None:
            resolved = storage.resolve_path(genome_build, "spliceai")
            if resolved:
                spliceai_scores = resolved

    if gene_annotation is not None and gnomad is not None:
        annotation_config = AnnotationConfig(
            gene_annotation_path=gene_annotation,
            gnomad_path=gnomad,
            clinvar_path=clinvar,
        )

    prioritization_config = PrioritizationConfig(
        cadd_scores_path=cadd_scores,
        revel_scores_path=revel_scores,
        spliceai_scores_path=spliceai_scores,
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
        use_bundles=use_bundles,
        genome_build=genome_build,
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


def _run_bundle_cli(argv: list[str]) -> None:
    """Handle the 'vartriage bundle' subcommand."""
    from vartriage.bundle.cli import add_bundle_subcommands, run_bundle_command

    parser = argparse.ArgumentParser(prog="vartriage bundle")
    subparsers = parser.add_subparsers(dest="bundle_command")
    add_bundle_subcommands(subparsers)

    # Re-parse with the bundle-specific parser
    # add_bundle_subcommands adds to a sub-parser, but we need
    # to parse directly here since we intercepted early.
    bundle_parser = argparse.ArgumentParser(prog="vartriage bundle")
    bundle_sub = bundle_parser.add_subparsers(dest="bundle_command")

    dl = bundle_sub.add_parser("download", help="Download a reference bundle")
    dl.add_argument("--bundle", required=True)
    dl.add_argument("--build", default=None)
    dl.add_argument("--dest", default=None)
    dl.add_argument("--no-transform", action="store_true")
    dl.add_argument("--no-progress", action="store_true")

    ls = bundle_sub.add_parser("list", help="List bundles")
    ls.add_argument("--build", default=None)
    ls.add_argument("--json", action="store_true", dest="json_output")

    vf = bundle_sub.add_parser("verify", help="Verify bundles")
    vf.add_argument("--bundle", default=None)
    vf.add_argument("--build", default=None)

    bundle_sub.add_parser("status", help="Show status")
    bundle_sub.add_parser("update-registry", help="Update registry")

    args = bundle_parser.parse_args(argv)
    exit_code = run_bundle_command(args)
    sys.exit(exit_code)
