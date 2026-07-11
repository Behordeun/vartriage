"""Pipeline configuration models with startup validation.

All configuration dataclasses are frozen (immutable after creation) and validate
their parameters in ``__post_init__``. Invalid values raise ``ValueError`` with
a message specifying the valid range.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional


@dataclass(frozen=True)
class QualityFilterConfig:
    """Configuration for quality-based variant filtering.

    Parameters
    ----------
    min_qual : float
        Minimum QUAL score threshold. Variants with a QUAL score below this
        value are excluded from downstream analysis. Must be in the range
        [0, 1_000_000]. Default is 20.0.

    Raises
    ------
    ValueError
        If ``min_qual`` is outside the range [0, 1_000_000].
    """

    min_qual: float = 20.0

    def __post_init__(self) -> None:
        if not (0 <= self.min_qual <= 1_000_000):
            raise ValueError(
                f"min_qual must be between 0 and 1000000, got {self.min_qual}"
            )


@dataclass(frozen=True)
class AnnotationConfig:
    """Configuration for the annotation engine.

    Parameters
    ----------
    gene_annotation_path : Path
        Path to a GTF/GFF gene annotation reference file used for functional
        consequence assignment via coordinate overlap.
    gnomad_path : Path
        Path to a local gnomAD reference file for population allele frequency
        lookups.
    clinvar_path : Optional[Path]
        Path to a ClinVar reference file for clinical significance lookups.
        When None, ClinVar annotation is skipped and variants receive a null
        clinical significance value.
    batch_size : int
        Number of variants processed per batch during vectorized annotation
        operations. Must be in the range [1_000, 100_000]. Default is 10_000.

    Raises
    ------
    ValueError
        If ``batch_size`` is outside the range [1_000, 100_000].
    """

    gene_annotation_path: Path
    gnomad_path: Path
    clinvar_path: Optional[Path] = None
    batch_size: int = 10_000

    def __post_init__(self) -> None:
        if not (1_000 <= self.batch_size <= 100_000):
            raise ValueError(
                f"batch_size must be between 1000 and 100000, got {self.batch_size}"
            )


@dataclass(frozen=True)
class PrioritizationConfig:
    """Configuration for the prioritization engine.

    Parameters
    ----------
    max_allele_frequency : float
        Maximum allele frequency threshold. Variants with a population
        frequency strictly above this value are excluded (unless they carry
        the ``frequency_unknown`` flag). Must be in the range [0.0, 1.0].
        Default is 0.01.
    cadd_scores_path : Optional[Path]
        Path to a CADD Phred score reference file. When None, CADD scores are
        not incorporated into composite ranking.
    revel_scores_path : Optional[Path]
        Path to a REVEL score reference file. When None, REVEL scores are not
        incorporated into composite ranking.
    spliceai_scores_path : Optional[Path]
        Path to a SpliceAI score TSV reference file. When None, SpliceAI
        scores are not incorporated into composite ranking.
    batch_size : int
        Number of variants processed per batch during vectorized score
        normalization. Must be in the range [1_000, 100_000]. Default is
        10_000.

    Raises
    ------
    ValueError
        If ``max_allele_frequency`` is outside the range [0.0, 1.0].
    ValueError
        If ``batch_size`` is outside the range [1_000, 100_000].
    """

    max_allele_frequency: float = 0.01
    cadd_scores_path: Optional[Path] = None
    revel_scores_path: Optional[Path] = None
    spliceai_scores_path: Optional[Path] = None
    batch_size: int = 10_000

    def __post_init__(self) -> None:
        if not (0.0 <= self.max_allele_frequency <= 1.0):
            raise ValueError(
                f"max_allele_frequency must be between 0.0 and 1.0, "
                f"got {self.max_allele_frequency}"
            )
        if not (1_000 <= self.batch_size <= 100_000):
            raise ValueError(
                f"batch_size must be between 1000 and 100000, got {self.batch_size}"
            )


@dataclass(frozen=True)
class ReportConfig:
    """Configuration for report generation.

    Parameters
    ----------
    output_format : Literal["json", "csv", "pdf", "vcf"]
        Desired output format for the final clinical report. Default is
        ``"json"``.
    """

    output_format: Literal["json", "csv", "pdf", "vcf"] = "json"


@dataclass(frozen=True)
class MissingDataConfig:
    """Configuration for missing data handling behavior.

    Parameters
    ----------
    warning_threshold : int
        Maximum number of ``MissingDataWarning`` events allowed before the
        pipeline emits a summary warning. The summary includes the total count
        of missing-data events and the reference sources that contributed.
        Default is 1000.
    """

    warning_threshold: int = 1000


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level pipeline configuration aggregating all sub-configs.

    Parameters
    ----------
    vcf_path : Path
        Path to the input VCF file (``.vcf`` or ``.vcf.gz``).
    output_path : Path
        Path where the output report file will be written.
    quality_filter : QualityFilterConfig
        Quality filtering settings. Defaults to standard thresholds.
    annotation : AnnotationConfig
        Annotation engine settings including reference file paths.
    prioritization : PrioritizationConfig
        Prioritization engine settings for frequency filtering and scoring.
    report : ReportConfig
        Report generation format settings.
    missing_data : MissingDataConfig
        Missing data handling and warning threshold settings.
    """

    vcf_path: Path
    output_path: Path
    quality_filter: QualityFilterConfig = field(default_factory=QualityFilterConfig)
    annotation: Optional[AnnotationConfig] = None
    prioritization: PrioritizationConfig = field(default_factory=PrioritizationConfig)
    report: ReportConfig = field(default_factory=ReportConfig)
    missing_data: MissingDataConfig = field(default_factory=MissingDataConfig)
