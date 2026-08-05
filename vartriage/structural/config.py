"""Configuration for the structural variant triage pipeline.

Frozen dataclasses with startup validation. Invalid values raise
ValueError with a message naming the valid range.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SVTriageConfig:
    """Top-level configuration for structural variant triage.

    Parameters
    ----------
    vcf_path : Path
        Input VCF containing structural variant calls.
    output_path : Path
        Destination for the triage report.
    gene_annotation_path : Optional[Path]
        GTF/GFF gene annotation for overlap assessment. Required for
        gene-level impact classification.
    dosage_sensitivity_path : Optional[Path]
        ClinGen dosage sensitivity map (TSV with gene, HI score, TS score).
        When None, dosage scoring is skipped.
    gnomad_sv_path : Optional[Path]
        gnomAD-SV reference for population frequency matching.
        When None, all SVs are treated as frequency-unknown.
    pathogenic_regions_path : Optional[Path]
        BED file of established pathogenic CNV regions for Section 2
        overlap scoring. When None, Section 2 evidence is skipped.
    benign_regions_path : Optional[Path]
        BED file of established benign CNV regions. When None, benign
        region overlap is not evaluated.
    min_sv_size : int
        Minimum SV length in bp to process. Smaller events are filtered
        out as likely indels handled by the point-variant pipeline.
        Must be in range [1, 10_000_000]. Default is 50.
    max_allele_frequency : float
        Maximum gnomAD-SV frequency for an SV to pass filtering.
        Must be in range [0.0, 1.0]. Default is 0.01.
    reciprocal_overlap : float
        Minimum reciprocal overlap fraction for matching SVs against
        the population frequency database. Must be in range [0.0, 1.0].
        Default is 0.5 (50%), standard in SV analysis.
    whole_gene_threshold : float
        Fraction of gene length that must be overlapped to classify
        as a whole-gene event. Must be in range [0.5, 1.0].
        Default is 0.8 (80%).
    min_quality : float
        Minimum QUAL score for SV calls. SVs below this threshold
        are excluded. Must be in range [0.0, 1_000_000]. Default is 20.0.
    batch_size : int
        Number of SVs processed per batch. Must be in range
        [100, 100_000]. Default is 5_000.
    output_format : str
        Report format. One of "json", "csv". Default is "json".
    include_benign : bool
        Whether to include Benign/Likely_Benign SVs in output.
        Default is False (only VUS and above).

    Raises
    ------
    ValueError
        If any parameter is outside its valid range.
    """

    vcf_path: Path
    output_path: Path
    gene_annotation_path: Path | None = None
    dosage_sensitivity_path: Path | None = None
    gnomad_sv_path: Path | None = None
    pathogenic_regions_path: Path | None = None
    benign_regions_path: Path | None = None
    min_sv_size: int = 50
    max_sv_size: int = 0
    max_allele_frequency: float = 0.01
    reciprocal_overlap: float = 0.5
    whole_gene_threshold: float = 0.8
    min_quality: float = 20.0
    batch_size: int = 5_000
    output_format: str = "json"
    include_benign: bool = False

    def __post_init__(self) -> None:
        if not (1 <= self.min_sv_size <= 10_000_000):
            raise ValueError(
                f"min_sv_size must be between 1 and 10000000, got {self.min_sv_size}"
            )
        if self.max_sv_size < 0:
            raise ValueError(
                f"max_sv_size must be >= 0 (0 means no limit), got {self.max_sv_size}"
            )
        if self.max_sv_size > 0 and self.max_sv_size < self.min_sv_size:
            raise ValueError(
                f"max_sv_size ({self.max_sv_size}) must be >= "
                f"min_sv_size ({self.min_sv_size})"
            )
        if not (0.0 <= self.max_allele_frequency <= 1.0):
            raise ValueError(
                f"max_allele_frequency must be between 0.0 and 1.0, "
                f"got {self.max_allele_frequency}"
            )
        if not (0.0 <= self.reciprocal_overlap <= 1.0):
            raise ValueError(
                f"reciprocal_overlap must be between 0.0 and 1.0, "
                f"got {self.reciprocal_overlap}"
            )
        if not (0.5 <= self.whole_gene_threshold <= 1.0):
            raise ValueError(
                f"whole_gene_threshold must be between 0.5 and 1.0, "
                f"got {self.whole_gene_threshold}"
            )
        if not (0.0 <= self.min_quality <= 1_000_000):
            raise ValueError(
                f"min_quality must be between 0.0 and 1000000, got {self.min_quality}"
            )
        if not (100 <= self.batch_size <= 100_000):
            raise ValueError(
                f"batch_size must be between 100 and 100000, got {self.batch_size}"
            )
        valid_formats = ("json", "csv")
        if self.output_format not in valid_formats:
            raise ValueError(
                f"output_format must be one of {valid_formats}, "
                f"got '{self.output_format}'"
            )
