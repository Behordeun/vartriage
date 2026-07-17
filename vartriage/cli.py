"""CLI entry point for the vartriage pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Literal, Optional, TypeVar, cast

from vartriage.models.config import ClinicalReportConfig, InheritanceConfig

if TYPE_CHECKING:
    from vartriage.models.config import SampleConfig


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
        "--reference-fasta",
        type=Path,
        default=None,
        help=(
            "Path to indexed reference genome FASTA (.fa + .fai). "
            "Enables codon-level consequence calling and variant normalization. "
            "Without this, the pipeline uses a positional heuristic for consequences."
        ),
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
    parser.add_argument(
        "--mode",
        choices=["local", "api", "hybrid"],
        default="local",
        help=(
            "Annotation mode: 'local' uses file-based backends (default), "
            "'api' queries remote services (Ensembl VEP, ClinVar), "
            "'hybrid' uses local files where available and API for gaps. "
            "API mode requires: pip install vartriage[api]"
        ),
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=None,
        help="NCBI API key for higher ClinVar rate limits (also reads NCBI_API_KEY env var)",
    )
    parser.add_argument(
        "--no-confirm",
        action="store_true",
        default=False,
        help="Skip confirmation prompts for large API-mode runs (>1000 variants)",
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

    use_bundles: bool = getattr(args, "use_bundles", False)
    genome_build: str = getattr(args, "genome_build", "grch38")
    mode: str = getattr(args, "mode", "local")

    # Build API config if mode requires it
    api_config = _build_api_config(args, mode, genome_build)

    paths = _resolve_reference_paths(args, use_bundles, genome_build)

    annotation_config: Optional[AnnotationConfig] = None
    if paths["gene_annotation"] is not None and paths["gnomad"] is not None:
        annotation_config = AnnotationConfig(
            gene_annotation_path=paths["gene_annotation"],
            gnomad_path=paths["gnomad"],
            clinvar_path=paths["clinvar"],
            reference_fasta_path=getattr(args, "reference_fasta", None),
        )

    prioritization_config = PrioritizationConfig(
        cadd_scores_path=paths["cadd_scores"],
        revel_scores_path=paths["revel_scores"],
        spliceai_scores_path=paths["spliceai_scores"],
    )

    report_config = ReportConfig(
        output_format=cast(
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
        ),
    )

    gene_filter_config = _build_optional_config(
        args.gene_list, lambda p: GeneFilterConfig(gene_list_path=p)
    )
    region_filter_config = _build_optional_config(
        args.regions, lambda p: RegionFilterConfig(bed_path=p)
    )
    sample_config = _build_sample_config(args)

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
        api=api_config,
    )

    pipeline = Pipeline(pipeline_config)
    return pipeline.run()


def _build_api_config(
    args: argparse.Namespace,
    mode: str,
    genome_build: str,
) -> "Optional[object]":
    """Build APIConfig if mode is api or hybrid. Returns None for local mode."""
    if mode == "local":
        return None

    try:
        from vartriage.api.config import APIConfig
    except ImportError:
        print(
            "Error: API mode requires the 'httpx' package.\n"
            "Install with: pip install vartriage[api]",
            file=sys.stderr,
        )
        sys.exit(1)

    api_key: Optional[str] = getattr(args, "api_key", None)

    return APIConfig.load(
        mode=mode,
        genome_build=genome_build,
        ncbi_api_key=api_key,
    )


def _resolve_reference_paths(
    args: argparse.Namespace, use_bundles: bool, genome_build: str
) -> dict[str, Optional[Path]]:
    """Resolve reference file paths, filling from bundles if enabled."""
    paths: dict[str, Optional[Path]] = {
        "gene_annotation": args.gene_annotation,
        "gnomad": args.gnomad,
        "clinvar": args.clinvar,
        "cadd_scores": args.cadd_scores,
        "revel_scores": args.revel_scores,
        "spliceai_scores": args.spliceai_scores,
    }

    if not use_bundles:
        return paths

    from vartriage.bundle.storage import BundleStorage

    storage = BundleStorage()

    bundle_names = {
        "gene_annotation": "gencode",
        "gnomad": "gnomad-exomes-chr22",
        "clinvar": "clinvar",
        "cadd_scores": "cadd",
        "revel_scores": "revel",
        "spliceai_scores": "spliceai",
    }

    for key, bundle_name in bundle_names.items():
        if paths[key] is None:
            resolved = storage.resolve_path(genome_build, bundle_name)
            if resolved:
                paths[key] = resolved

    return paths


ConfigT = TypeVar("ConfigT")


def _build_optional_config(
    value: Optional[Path], factory: Callable[[Path], ConfigT]
) -> Optional[ConfigT]:
    """Build an optional config if value is not None."""
    return factory(value) if value is not None else None


def _build_sample_config(
    args: argparse.Namespace,
) -> Optional["SampleConfig"]:
    """Build SampleConfig from args, validating --min-gq requires --sample."""
    from vartriage.models.config import SampleConfig

    if args.min_gq is not None and args.sample is None:
        print("Error: --min-gq requires --sample", file=sys.stderr)
        sys.exit(2)

    if args.sample is None:
        return None

    return SampleConfig(sample_name=args.sample, min_gq=args.min_gq)


def _build_clinical_config(
    args: argparse.Namespace,
    output_format: str,
) -> Optional[ClinicalReportConfig]:
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
    args = parser.parse_args(argv)
    exit_code = run_bundle_command(args)
    sys.exit(exit_code)
