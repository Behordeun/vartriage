"""CLI entry point for the vartriage pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional


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
        choices=["json", "csv", "pdf", "vcf"],
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
        "--version",
        action="version",
        version=f"%(prog)s {_get_version()}",
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

    vcf_path: Path = args.vcf
    if not vcf_path.exists():
        print(
            f"Error: VCF file not found: {vcf_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        result_path = _run_pipeline(args, vcf_path)
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
        from vartriage.io.exceptions import (
            VariantPrioritizationError,
        )

        if isinstance(exc, VariantPrioritizationError):
            print(
                f"Error: pipeline failed: {exc}",
                file=sys.stderr,
            )
        else:
            print(
                f"Error: unexpected failure: {exc}",
                file=sys.stderr,
            )
        sys.exit(1)

    print(str(result_path))
    sys.exit(0)


def _run_pipeline(args: argparse.Namespace, vcf_path: Path) -> Path:
    """Assemble pipeline config from parsed args and run it.

    Returns
    -------
    Path
        Path to the generated report.
    """
    from vartriage.models.config import (
        AnnotationConfig,
        PipelineConfig,
        PrioritizationConfig,
        ReportConfig,
    )
    from vartriage.pipeline import Pipeline

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
    )

    report_config = ReportConfig(
        output_format=args.output_format,
    )

    pipeline_config = PipelineConfig(
        vcf_path=vcf_path,
        output_path=args.output,
        annotation=annotation_config,
        prioritization=prioritization_config,
        report=report_config,
    )

    pipeline = Pipeline(pipeline_config)
    return pipeline.run()
