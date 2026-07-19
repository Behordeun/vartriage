"""Cohort analysis orchestration and manifest parsing.

Separated from CLI to keep argument parsing thin and to allow
direct programmatic invocation and unit testing without going
through argparse or sys.exit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

from vartriage.cohort.pipeline import CohortPipeline
from vartriage.models.cohort import CohortConfig
from vartriage.models.config import (
    AnnotationConfig,
    GeneFilterConfig,
    PipelineConfig,
    PrioritizationConfig,
    ReportConfig,
)


@dataclass(frozen=True)
class CohortCLIConfig:
    """Typed configuration for a cohort CLI invocation.

    Captures all parsed arguments in a structured form so that
    orchestration logic can be tested without argparse.
    """

    sample_vcfs: list[Path]
    output: Path
    cohort_name: str = "cohort"
    output_format: Literal["json", "csv"] = "json"
    min_recurrence: int = 2
    max_af: float = 0.05
    include_singletons: bool = True
    parallel: bool = False
    max_workers: int = 4
    use_bundles: bool = False
    genome_build: str = "grch38"
    gene_list: Optional[Path] = None
    gene_annotation: Optional[Path] = None
    gnomad: Optional[Path] = None
    clinvar: Optional[Path] = None
    cadd_scores: Optional[Path] = None
    revel_scores: Optional[Path] = None
    spliceai_scores: Optional[Path] = None
    sample_labels: Optional[dict[str, str]] = field(default=None)


def run_cohort(config: CohortCLIConfig) -> list[Path]:
    """Build cohort pipeline from CLI config and execute.

    Parameters
    ----------
    config : CohortCLIConfig
        Typed configuration from CLI argument parsing.

    Returns
    -------
    list[Path]
        Paths to generated report files.

    Raises
    ------
    FileNotFoundError
        If a reference file or VCF does not exist.
    ValueError
        If configuration is invalid.
    """
    paths = _resolve_paths(config)

    annotation_config: Optional[AnnotationConfig] = None
    if paths["gene_annotation"] is not None and paths["gnomad"] is not None:
        annotation_config = AnnotationConfig(
            gene_annotation_path=paths["gene_annotation"],
            gnomad_path=paths["gnomad"],
            clinvar_path=paths["clinvar"],
        )

    prioritization_config = PrioritizationConfig(
        cadd_scores_path=paths["cadd_scores"],
        revel_scores_path=paths["revel_scores"],
        spliceai_scores_path=paths["spliceai_scores"],
    )

    gene_filter_config = None
    if config.gene_list is not None:
        gene_filter_config = GeneFilterConfig(gene_list_path=config.gene_list)

    base_pipeline_config = PipelineConfig(
        vcf_path=config.sample_vcfs[0],
        output_path=config.output / "tmp.json",
        annotation=annotation_config,
        prioritization=prioritization_config,
        report=ReportConfig(output_format="json"),
        gene_filter=gene_filter_config,
    )

    cohort_config = CohortConfig(
        sample_vcfs=config.sample_vcfs,
        output_path=config.output,
        cohort_name=config.cohort_name,
        min_recurrence=config.min_recurrence,
        output_format=config.output_format,
        max_af_threshold=config.max_af,
        include_singletons=config.include_singletons,
        sample_labels=config.sample_labels,
        parallel=config.parallel,
        max_workers=config.max_workers,
    )

    pipeline = CohortPipeline(
        cohort_config=cohort_config,
        pipeline_config=base_pipeline_config,
        annotation_config=annotation_config,
        prioritization_config=prioritization_config,
    )
    return pipeline.run()


def parse_cohort_manifest(
    manifest_path: Path,
) -> tuple[list[Path], dict[str, str] | None]:
    """Parse a cohort manifest file.

    Format: one VCF path per line. Optional tab-separated second column
    provides a human-readable sample label. Lines starting with '#' are
    comments. Blank lines are skipped.

    Parameters
    ----------
    manifest_path : Path
        Path to the manifest file.

    Returns
    -------
    tuple[list[Path], dict[str, str] | None]
        (list of VCF paths, optional label mapping keyed by file stem)

    Raises
    ------
    FileNotFoundError
        If the manifest file does not exist.
    """
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

    vcf_paths: list[Path] = []
    labels: dict[str, str] = {}
    has_labels = False

    with open(manifest_path, encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            vcf_path = Path(parts[0].strip())

            if not vcf_path.is_absolute():
                vcf_path = manifest_path.parent / vcf_path

            vcf_paths.append(vcf_path)

            if len(parts) >= 2:
                label = parts[1].strip()
                if label:
                    stem = _stem_from_path(vcf_path)
                    labels[stem] = label
                    has_labels = True

    return vcf_paths, labels if has_labels else None


def _stem_from_path(vcf_path: Path) -> str:
    """Extract the base stem from a VCF path, stripping .vcf suffix."""
    stem = vcf_path.stem
    if stem.endswith(".vcf"):
        stem = stem[:-4]
    return stem


def _resolve_paths(config: CohortCLIConfig) -> dict[str, Optional[Path]]:
    """Resolve reference file paths, filling from bundles if enabled."""
    paths: dict[str, Optional[Path]] = {
        "gene_annotation": config.gene_annotation,
        "gnomad": config.gnomad,
        "clinvar": config.clinvar,
        "cadd_scores": config.cadd_scores,
        "revel_scores": config.revel_scores,
        "spliceai_scores": config.spliceai_scores,
    }

    if not config.use_bundles:
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
            resolved = storage.resolve_path(config.genome_build, bundle_name)
            if resolved:
                paths[key] = resolved

    return paths
