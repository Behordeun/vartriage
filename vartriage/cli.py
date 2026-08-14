"""CLI entry point for the vartriage pipeline."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeVar, cast

from vartriage._internal.path_safety import resolve_path
from vartriage.models.config import ClinicalReportConfig, InheritanceConfig

if TYPE_CHECKING:
    from vartriage.knowledge.config import KnowledgeBaseConfig
    from vartriage.mito.config import MitoConfig
    from vartriage.models.config import SampleConfig
    from vartriage.remote.config import RemoteTabixConfig
    from vartriage.remote.presets import PresetEntry


def _get_version() -> str:
    """Return the installed package version, falling back to __version__."""
    try:
        from importlib.metadata import version

        return version("vartriage")
    except Exception:
        from vartriage import __version__

        return __version__


def _validated_heteroplasmy(value: str) -> float:
    """Argparse type validator for --mt-min-heteroplasmy."""
    try:
        f = float(value)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"invalid float value: '{value}'") from err
    if not (0.0 <= f <= 100.0):
        raise argparse.ArgumentTypeError(f"must be between 0.0 and 100.0, got {f}")
    return f


def _add_reference_arguments(parser: argparse.ArgumentParser) -> None:
    """Add shared reference file and bundle arguments to a parser.

    Used by both the main parser and the cohort subcommand to keep
    argument definitions, defaults, and help text in sync.
    """
    parser.add_argument(
        "--gene-annotation",
        type=Path,
        default=None,
        help=(
            "Path to GTF/GFF gene annotation reference file. "
            "Required together with --gnomad for annotation."
        ),
    )
    parser.add_argument(
        "--gnomad",
        type=Path,
        default=None,
        help=(
            "Path to gnomAD population frequency reference file. "
            "Required together with --gene-annotation for annotation."
        ),
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
    _add_reference_arguments(parser)

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
    parser.add_argument(
        "--secondary-findings",
        action="store_true",
        default=False,
        help=(
            "Flag variants in ACMG Secondary Findings (SF v3.2) genes "
            "regardless of primary gene panel filter. Adds a dedicated "
            "section to clinical reports."
        ),
    )

    # Gene-disease linkage arguments
    parser.add_argument(
        "--hpo-terms",
        type=str,
        default=None,
        help=(
            "Comma-separated HPO term IDs for phenotype-driven prioritization "
            "(e.g., HP:0001250,HP:0001249). Genes matching patient phenotype "
            "receive a ranking boost."
        ),
    )
    parser.add_argument(
        "--inheritance-mode",
        type=str,
        default=None,
        help=(
            "Filter variants to genes matching this inheritance mode "
            "(AD, AR, XL, XLD, XLR, MT). Genes without OMIM data or "
            "intergenic variants pass through unfiltered."
        ),
    )
    parser.add_argument(
        "--flag-actionable",
        action="store_true",
        default=False,
        help=(
            "Flag variants in medically actionable genes (ClinGen "
            "actionability curations). Adds is_actionable to output."
        ),
    )
    parser.add_argument(
        "--knowledge-dir",
        type=Path,
        default=None,
        help=(
            "Path to custom gene knowledge data directory. "
            "Defaults to bundled package data."
        ),
    )

    # Mitochondrial analysis options
    parser.add_argument(
        "--skip-mito",
        action="store_true",
        default=False,
        help=(
            "Skip mitochondrial variant analysis. Use for targeted panels "
            "without mtDNA capture where chrM variants are noise."
        ),
    )
    parser.add_argument(
        "--mt-min-heteroplasmy",
        type=_validated_heteroplasmy,
        default=1.0,
        help=(
            "Minimum heteroplasmy percentage for reporting mtDNA variants "
            "(0.0-100.0). Variants below this threshold are filtered as "
            "sub-threshold noise. Default: 1.0%%"
        ),
    )

    # Structural variant analysis (integrated mode)
    parser.add_argument(
        "--sv-vcf",
        type=Path,
        default=None,
        help=(
            "Path to a VCF file containing structural variant calls. "
            "When provided alongside --vcf, both SNV and SV pipelines "
            "run and produce separate output files."
        ),
    )

    # Remote tabix score backend
    parser.add_argument(
        "--cadd-remote",
        type=str,
        default=None,
        help=(
            "Remote CADD score source: a named preset (e.g., cadd-v1.7-grch38) "
            "or a full URL to a bgzipped/tabix-indexed CADD TSV. Queries use "
            "HTTP byte-range requests — no local download needed. "
            "Ignored when --cadd-scores is provided (local file takes priority)."
        ),
    )
    parser.add_argument(
        "--gnomad-remote",
        type=str,
        default=None,
        help=(
            "Remote gnomAD frequency source: a named preset "
            "(e.g., gnomad-exomes-v4-grch38) or a URL template with {chrom} "
            "placeholder. Queries per-chromosome VCFs via HTTP byte-range. "
            "Ignored when --gnomad is provided (local file takes priority)."
        ),
    )
    parser.add_argument(
        "--remote-cache-ttl",
        type=int,
        default=30,
        help=(
            "Remote score cache TTL in days. Scores fetched from remote "
            "tabix are cached locally to avoid redundant network requests. "
            "Use -1 for pinned mode (never expire, clinical reproducibility). "
            "Default: 30"
        ),
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the vartriage CLI.

    Parameters
    ----------
    argv : list[str], optional
        Arguments to parse. Uses sys.argv[1:] when None.
    """
    # Intercept subcommands before main parser
    effective_argv = argv if argv is not None else sys.argv[1:]
    if effective_argv and effective_argv[0] == "bundle":
        _run_bundle_cli(effective_argv[1:])
        return
    if effective_argv and effective_argv[0] == "cohort":
        _run_cohort_cli(effective_argv[1:])
        return
    if effective_argv and effective_argv[0] == "sv":
        _run_sv_cli(effective_argv[1:])
        return
    if effective_argv and effective_argv[0] == "remote":
        _run_remote_cli(effective_argv[1:])
        return

    parser = _build_parser()
    args = parser.parse_args(argv)

    # Validate clinical format requirements before VCF check.
    # _build_clinical_config is the single source of truth for
    # this validation and will sys.exit(2) on missing flags.
    output_fmt: str = args.output_format
    clinical_config = _build_clinical_config(args, output_fmt)

    vcf_path: Path = resolve_path(args.vcf)
    if not vcf_path.exists():
        print(
            f"Error: VCF file not found: {vcf_path}",  # nosec: output to stderr, not rendered in browser
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        result_path = _run_pipeline(args, vcf_path, clinical_config)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)  # nosec: output to stderr, not rendered in browser
        sys.exit(1)
    except OSError as exc:
        print(
            f"Error: report generation failed: {exc}",  # nosec: output to stderr, not rendered in browser
            file=sys.stderr,
        )
        sys.exit(1)
    except Exception as exc:
        _handle_unexpected_error(exc)

    print(result_path)
    sys.exit(0)


def _handle_unexpected_error(exc: Exception) -> None:
    """Print an appropriate error message and exit."""
    from vartriage.io.exceptions import VariantPrioritizationError

    if isinstance(exc, VariantPrioritizationError):
        print(f"Error: pipeline failed: {exc}", file=sys.stderr)  # nosec: output to stderr, not rendered in browser
    else:
        print(f"Error: unexpected failure: {exc}", file=sys.stderr)  # nosec: output to stderr, not rendered in browser
    sys.exit(1)


def _run_pipeline(
    args: argparse.Namespace,
    vcf_path: Path,
    clinical_config: ClinicalReportConfig | None = None,
) -> Path:
    """Assemble pipeline config from parsed args and run it.

    Returns
    -------
    Path
        Path to the generated report.
    """
    from vartriage.models.config import (
        AnnotationConfig,
        GeneFilterConfig,
        PipelineConfig,
        PrioritizationConfig,
        RegionFilterConfig,
        ReportConfig,
    )
    from vartriage.pipeline import Pipeline

    output_format: str = args.output_format
    inheritance_config = _build_inheritance_config(args)

    use_bundles: bool = getattr(args, "use_bundles", False)
    genome_build: str = getattr(args, "genome_build", "grch38")
    mode: str = getattr(args, "mode", "local")

    # Build API config if mode requires it
    api_config = _build_api_config(args, mode, genome_build)

    paths = _resolve_reference_paths(args, use_bundles, genome_build)

    annotation_config: AnnotationConfig | None = None
    if paths["gene_annotation"] is not None and paths["gnomad"] is not None:
        annotation_config = AnnotationConfig(
            gene_annotation_path=paths["gene_annotation"],
            gnomad_path=paths["gnomad"],
            clinvar_path=paths["clinvar"],
            reference_fasta_path=getattr(args, "reference_fasta", None),
        )
    elif paths["gene_annotation"] is not None:
        # gene annotation available but no local gnomAD — still allow
        # annotation if remote gnomAD is configured (frequency comes
        # from the remote backend injected by the pipeline)
        annotation_config = AnnotationConfig(
            gene_annotation_path=paths["gene_annotation"],
            gnomad_path=None,
            clinvar_path=paths["clinvar"],
            reference_fasta_path=getattr(args, "reference_fasta", None),
        )

    # Build remote tabix config if any remote flags are set
    remote_config = _build_remote_config(args)

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

    knowledge_config = _build_knowledge_config(args)

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
        knowledge=knowledge_config,
        sv_vcf_path=getattr(args, "sv_vcf", None),
        mito=_build_mito_config(args),
        remote=remote_config,
    )

    pipeline = Pipeline(pipeline_config)
    if pipeline_config.sv_vcf_path is not None:
        pipeline.run_with_sv()
        return pipeline_config.output_path
    return pipeline.run()


def _build_api_config(
    args: argparse.Namespace,
    mode: str,
    genome_build: str,
) -> object | None:
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

    api_key: str | None = getattr(args, "api_key", None)

    return APIConfig.load(
        mode=mode,
        genome_build=genome_build,
        ncbi_api_key=api_key,
    )


def _resolve_reference_paths(
    args: argparse.Namespace, use_bundles: bool, genome_build: str
) -> dict[str, Path | None]:
    """Resolve reference file paths, filling from bundles if enabled."""
    paths: dict[str, Path | None] = {
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
        if paths[key] is None and (
            resolved := storage.resolve_path(genome_build, bundle_name)
        ):
            paths[key] = resolved

    return paths


ConfigT = TypeVar("ConfigT")


def _build_optional_config(
    value: Path | None, factory: Callable[[Path], ConfigT]
) -> ConfigT | None:
    """Build an optional config if value is not None."""
    return factory(value) if value is not None else None


def _build_sample_config(
    args: argparse.Namespace,
) -> SampleConfig | None:
    """Build SampleConfig from args, validating --min-gq requires --sample."""
    from vartriage.models.config import SampleConfig

    if args.min_gq is not None and args.sample is None:
        print("Error: --min-gq requires --sample", file=sys.stderr)
        sys.exit(2)

    if args.sample is None:
        return None

    return SampleConfig(sample_name=args.sample, min_gq=args.min_gq)


def _build_knowledge_config(
    args: argparse.Namespace,
) -> KnowledgeBaseConfig | None:
    """Build KnowledgeBaseConfig if gene-disease linkage features are requested."""
    from vartriage.knowledge.config import KnowledgeBaseConfig

    hpo_terms_raw: str | None = getattr(args, "hpo_terms", None)
    knowledge_dir: Path | None = getattr(args, "knowledge_dir", None)
    flag_actionable: bool = getattr(args, "flag_actionable", False)
    inheritance_mode: str | None = getattr(args, "inheritance_mode", None)

    # Only create config when at least one gene-knowledge feature is active
    has_knowledge_request = (
        hpo_terms_raw is not None
        or knowledge_dir is not None
        or flag_actionable
        or inheritance_mode is not None
    )

    if not has_knowledge_request:
        return None

    hpo_terms: frozenset[str] = frozenset()
    if hpo_terms_raw:
        parsed = [t.strip() for t in hpo_terms_raw.split(",") if t.strip()]
        hpo_terms = frozenset(parsed)

    return KnowledgeBaseConfig(
        data_dir=knowledge_dir,
        hpo_terms=hpo_terms,
        inheritance_mode=inheritance_mode,
        flag_actionable=flag_actionable,
    )


def _build_mito_config(
    args: argparse.Namespace,
) -> MitoConfig | None:
    """Build MitoConfig from CLI arguments.

    Returns None when --skip-mito is not set and defaults are acceptable
    (auto-detection handles the rest). Returns a MitoConfig with
    enabled=False when --skip-mito is passed.
    """
    from vartriage.mito.config import MitoConfig

    if getattr(args, "skip_mito", False):
        return MitoConfig(enabled=False)

    # Only create explicit config if non-default threshold was set
    if abs((threshold := getattr(args, "mt_min_heteroplasmy", 1.0)) - 1.0) > 1e-9:
        return MitoConfig(min_heteroplasmy=threshold)

    return None


def _build_clinical_config(
    args: argparse.Namespace,
    output_format: str,
) -> ClinicalReportConfig | None:
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
) -> InheritanceConfig | None:
    """Build InheritanceConfig from trio arguments."""
    from vartriage.models.config import InheritanceConfig as _IC

    proband: str | None = args.proband
    mother: str | None = args.mother
    father: str | None = args.father
    sample: str | None = getattr(args, "sample", None)

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


def _run_cohort_cli(argv: list[str]) -> None:
    """Handle the 'vartriage cohort' subcommand.

    Provides multi-sample cohort analysis: processes multiple VCF files,
    aggregates variants across samples, computes recurrence frequencies
    and per-gene burden, then writes cohort-level reports.
    """
    parser = argparse.ArgumentParser(
        prog="vartriage cohort",
        description=(
            "Multi-sample cohort analysis. Processes multiple VCF files "
            "through the standard pipeline, then aggregates variants "
            "across samples to identify shared mutations, compute "
            "recurrence frequencies, and generate cohort-level reports."
        ),
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=(
            "Path to a manifest file listing sample VCFs (one path per line). "
            "Lines starting with '#' are comments. Optional second column "
            "(tab-separated) provides a sample label."
        ),
    )
    input_group.add_argument(
        "--vcf",
        type=Path,
        nargs="+",
        default=None,
        help="Two or more VCF file paths to include in the cohort",
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for cohort report files",
    )
    parser.add_argument(
        "--cohort-name",
        type=str,
        default="cohort",
        help="Cohort identifier used in output filenames (default: cohort)",
    )
    parser.add_argument(
        "--output-format",
        choices=["json", "csv"],
        default="json",
        help="Output format for cohort reports (default: json)",
    )
    parser.add_argument(
        "--min-recurrence",
        type=int,
        default=2,
        help=(
            "Minimum number of samples a variant must appear in "
            "to be included in the output (default: 2). Variants "
            "below this threshold are excluded."
        ),
    )
    parser.add_argument(
        "--max-af",
        type=float,
        default=0.05,
        help=(
            "Maximum population allele frequency threshold for "
            "cohort inclusion (default: 0.05)"
        ),
    )
    parser.add_argument(
        "--no-singletons",
        action="store_true",
        default=False,
        help="Exclude singleton variants (appearing in only one sample) from output",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        default=False,
        help="Process samples concurrently using a thread pool",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Maximum parallel workers when --parallel is set (default: 4)",
    )

    # Reference file and bundle options (shared definition)
    _add_reference_arguments(parser)

    args = parser.parse_args(argv)

    # Resolve sample VCF list
    from vartriage.cohort.runner import (
        CohortCLIConfig,
        parse_cohort_manifest,
        run_cohort,
    )

    sample_vcfs: list[Path]
    sample_labels: dict[str, str] | None = None

    if args.manifest is not None:
        if not args.manifest.resolve().exists():
            print(
                f"Error: manifest file not found: {args.manifest}",
                file=sys.stderr,
            )
            sys.exit(1)
        sample_vcfs, sample_labels = parse_cohort_manifest(args.manifest.resolve())
    else:
        sample_vcfs = args.vcf

    if len(sample_vcfs) < 2:
        print(
            "Error: cohort analysis requires at least 2 sample VCFs",
            file=sys.stderr,
        )
        sys.exit(2)

    for vcf_path in sample_vcfs:
        if not vcf_path.exists():
            print(f"Error: VCF file not found: {vcf_path}", file=sys.stderr)
            sys.exit(1)

    config = CohortCLIConfig(
        sample_vcfs=sample_vcfs,
        output=args.output,
        cohort_name=args.cohort_name,
        output_format=args.output_format,
        min_recurrence=args.min_recurrence,
        max_af=args.max_af,
        include_singletons=not args.no_singletons,
        parallel=args.parallel,
        max_workers=args.max_workers,
        use_bundles=getattr(args, "use_bundles", False),
        genome_build=getattr(args, "genome_build", "grch38"),
        gene_list=args.gene_list,
        gene_annotation=args.gene_annotation,
        gnomad=args.gnomad,
        clinvar=args.clinvar,
        cadd_scores=args.cadd_scores,
        revel_scores=args.revel_scores,
        spliceai_scores=args.spliceai_scores,
        sample_labels=sample_labels,
    )

    try:
        report_paths = run_cohort(config)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(2)
    except Exception as exc:
        print(f"Error: cohort analysis failed: {exc}", file=sys.stderr)
        sys.exit(1)

    for p in report_paths:
        print(p)
    sys.exit(0)


def _run_sv_cli(argv: list[str]) -> None:
    """Handle the 'vartriage sv' subcommand for structural variant triage.

    Runs the SV triage pipeline: parse SV VCF, annotate gene overlap,
    score by dosage sensitivity, classify via ClinGen framework, and
    write a triage report.
    """
    parser = argparse.ArgumentParser(
        prog="vartriage sv",
        description=(
            "Structural variant triage pipeline. Parses SV calls from VCF, "
            "annotates gene overlap and dosage sensitivity, scores pathogenicity, "
            "classifies via ClinGen 2020 framework, and writes a prioritized report."
        ),
    )

    parser.add_argument(
        "--sv-vcf",
        type=Path,
        required=True,
        help="Path to a VCF file containing structural variant calls",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path where the SV triage report will be written",
    )
    parser.add_argument(
        "--output-format",
        choices=["json", "csv"],
        default="json",
        help="Output report format (default: json)",
    )
    parser.add_argument(
        "--gene-annotation",
        type=Path,
        default=None,
        help="Path to GTF/GFF gene annotation for gene overlap assessment",
    )
    parser.add_argument(
        "--dosage-sensitivity",
        type=Path,
        default=None,
        help="Path to ClinGen dosage sensitivity TSV (gene_symbol, hi_score, ts_score)",
    )
    parser.add_argument(
        "--gnomad-sv",
        type=Path,
        default=None,
        help="Path to gnomAD-SV reference (BED/TSV: chrom, start, end, sv_type, af)",
    )
    parser.add_argument(
        "--pathogenic-regions",
        type=Path,
        default=None,
        help="BED file of known pathogenic CNV regions for classification",
    )
    parser.add_argument(
        "--benign-regions",
        type=Path,
        default=None,
        help="BED file of known benign CNV regions for classification",
    )
    parser.add_argument(
        "--min-sv-size",
        type=int,
        default=50,
        help="Minimum SV size in bp to include (default: 50)",
    )
    parser.add_argument(
        "--max-sv-size",
        type=int,
        default=0,
        help="Maximum SV size in bp to include (default: 0, no limit)",
    )
    parser.add_argument(
        "--sv-types",
        type=str,
        default=None,
        help=(
            "Comma-separated list of SV types to include "
            "(DEL,DUP,INV,INS,BND,CNV). Default: all types."
        ),
    )
    parser.add_argument(
        "--max-af",
        type=float,
        default=0.01,
        help="Maximum population frequency for SV inclusion (default: 0.01)",
    )
    parser.add_argument(
        "--min-quality",
        type=float,
        default=20.0,
        help="Minimum QUAL score for SV calls (default: 20.0)",
    )
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        default=0.5,
        help="Minimum reciprocal overlap for frequency matching (default: 0.5)",
    )
    parser.add_argument(
        "--whole-gene-threshold",
        type=float,
        default=0.8,
        help="Gene overlap fraction to classify as whole-gene event (default: 0.8)",
    )
    parser.add_argument(
        "--include-benign",
        action="store_true",
        default=False,
        help="Include Benign/Likely_Benign SVs in output (excluded by default)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {_get_version()}",
    )

    args = parser.parse_args(argv)

    sv_vcf: Path = args.sv_vcf
    if not sv_vcf.exists():
        print(f"Error: SV VCF file not found: {sv_vcf}", file=sys.stderr)
        sys.exit(1)

    from vartriage.structural.config import SVTriageConfig
    from vartriage.structural.pipeline import SVTriagePipeline

    try:
        config = SVTriageConfig(
            vcf_path=sv_vcf,
            output_path=args.output,
            gene_annotation_path=args.gene_annotation,
            dosage_sensitivity_path=args.dosage_sensitivity,
            gnomad_sv_path=args.gnomad_sv,
            pathogenic_regions_path=args.pathogenic_regions,
            benign_regions_path=args.benign_regions,
            min_sv_size=args.min_sv_size,
            max_allele_frequency=args.max_af,
            reciprocal_overlap=args.reciprocal_overlap,
            whole_gene_threshold=args.whole_gene_threshold,
            min_quality=args.min_quality,
            output_format=args.output_format,
            include_benign=args.include_benign,
        )
    except ValueError as exc:
        print(f"Error: invalid configuration: {exc}", file=sys.stderr)
        sys.exit(2)

    try:
        pipeline = SVTriagePipeline(config)
        result_path = pipeline.run()
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f"Error: SV triage failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print(result_path)
    sys.exit(0)


def _build_remote_config(
    args: argparse.Namespace,
) -> RemoteTabixConfig | None:
    """Build RemoteTabixConfig from CLI arguments.

    Returns None when no remote flags are set, preserving existing
    behavior for users who don't use remote tabix.
    """
    from vartriage.remote.config import RemoteTabixConfig

    cadd_remote: str | None = getattr(args, "cadd_remote", None)
    gnomad_remote: str | None = getattr(args, "gnomad_remote", None)
    cache_ttl: int = getattr(args, "remote_cache_ttl", 30)

    if cadd_remote is None and gnomad_remote is None:
        return None

    # Local file takes priority over remote — clear remote if local exists
    if getattr(args, "cadd_scores", None) is not None and cadd_remote is not None:
        _log_remote_override(
            "--cadd-scores takes priority over --cadd-remote; remote CADD disabled"
        )
        cadd_remote = None

    if getattr(args, "gnomad", None) is not None and gnomad_remote is not None:
        _log_remote_override(
            "--gnomad takes priority over --gnomad-remote; remote gnomAD disabled"
        )
        gnomad_remote = None

    if cadd_remote is None and gnomad_remote is None:
        return None

    return RemoteTabixConfig(
        cadd_remote_url=cadd_remote,
        gnomad_remote_url=gnomad_remote,
        cache_ttl_days=cache_ttl,
    )


def _log_remote_override(message: str) -> None:
    import logging

    logging.getLogger(__name__).info(message)


def _print_presets(presets: list[PresetEntry]) -> None:
    """Print remote tabix presets to stdout."""
    if not presets:
        print("No presets found.", file=sys.stderr)
        return
    print(f"{'Name':<30} {'Source':<8} {'Build':<8} Description")
    print("-" * 80)
    for p in presets:
        print(f"{p.name:<30} {p.source:<8} {p.genome_build:<8} {p.description}")


def _run_remote_cli(argv: list[str]) -> None:
    """Handle the 'vartriage remote' subcommand.

    Currently supports:
        vartriage remote list-presets [--source cadd|gnomad]
    """
    from vartriage.remote.presets import list_presets

    parser = argparse.ArgumentParser(
        prog="vartriage remote",
        description="Manage remote tabix score backends.",
    )
    subparsers = parser.add_subparsers(dest="remote_command")

    list_parser = subparsers.add_parser(
        "list-presets",
        help="List available named presets for remote tabix databases",
    )
    list_parser.add_argument(
        "--source",
        type=str,
        choices=["cadd", "gnomad"],
        default=None,
        help="Filter presets by source type",
    )

    args = parser.parse_args(argv)

    if args.remote_command is None:
        parser.print_help()
        sys.exit(0)

    if args.remote_command == "list-presets":
        _print_presets(list_presets(source=args.source))
        sys.exit(0)
