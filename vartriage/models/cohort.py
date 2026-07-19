"""Data models for multi-sample cohort analysis.

Defines the configuration, per-variant cohort records, sample-level
occurrence tracking, and cohort-wide summary statistics. All models
are frozen dataclasses with startup validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

from vartriage.models.variant import (
    ACMGClassification,
    ClassifiedVariant,
    EvidenceTag,
    FunctionalConsequence,
)


@dataclass(frozen=True)
class CohortConfig:
    """Configuration for multi-sample cohort analysis.

    Parameters
    ----------
    sample_vcfs : list[Path]
        Paths to individual sample VCF files. Each file is processed
        independently through the standard pipeline, then results are
        merged. At least 2 samples are required for meaningful cohort
        analysis.
    output_path : Path
        Directory where cohort output files are written.
    cohort_name : str
        Human-readable cohort identifier used in report headers.
    min_recurrence : int
        Minimum number of samples a variant must appear in to be
        included in the cohort report. Must be >= 1.
        Default is 2 (shared by at least two samples).
    output_format : str
        Report output format. Default is "json".
    max_af_threshold : float
        Maximum population allele frequency for cohort inclusion.
        Variants exceeding this in gnomAD are excluded from cohort
        aggregation. Range [0.0, 1.0]. Default is 0.05.
    include_singletons : bool
        When True, variants appearing in only one sample are included
        in the full output (marked as singletons). When False, only
        recurrent variants (count >= min_recurrence) appear.
        Default is True.
    sample_labels : dict[str, str] | None
        Optional mapping of VCF file stems to human-readable sample
        labels for reports. Keys are Path.stem values.
    parallel : bool
        When True, process samples concurrently using a thread pool.
        Default is False (sequential processing).
    max_workers : int
        Maximum parallel workers when parallel=True. Must be >= 1.
        Default is 4.

    Raises
    ------
    ValueError
        If fewer than 2 sample VCFs are provided, min_recurrence < 1,
        max_af_threshold outside [0.0, 1.0], or max_workers < 1.
    """

    sample_vcfs: list[Path]
    output_path: Path
    cohort_name: str = "cohort"
    min_recurrence: int = 2
    output_format: Literal["json", "csv"] = "json"
    max_af_threshold: float = 0.05
    include_singletons: bool = True
    sample_labels: dict[str, str] | None = None
    parallel: bool = False
    max_workers: int = 4

    def __post_init__(self) -> None:
        if len(self.sample_vcfs) < 2:
            raise ValueError(
                f"Cohort analysis requires at least 2 samples, "
                f"got {len(self.sample_vcfs)}"
            )
        if self.min_recurrence < 1:
            raise ValueError(
                f"min_recurrence must be >= 1, got {self.min_recurrence}"
            )
        if not (0.0 <= self.max_af_threshold <= 1.0):
            raise ValueError(
                f"max_af_threshold must be between 0.0 and 1.0, "
                f"got {self.max_af_threshold}"
            )
        if self.max_workers < 1:
            raise ValueError(
                f"max_workers must be >= 1, got {self.max_workers}"
            )

    @property
    def sample_count(self) -> int:
        """Total number of samples in the cohort."""
        return len(self.sample_vcfs)

    def label_for(self, vcf_path: Path) -> str:
        """Resolve the display label for a sample VCF.

        Uses sample_labels mapping if available, otherwise falls
        back to the file stem.
        """
        stem = vcf_path.stem
        # Strip .vcf from stems like "sample.vcf.gz" -> "sample.vcf" -> "sample"
        if stem.endswith(".vcf"):
            stem = stem[:-4]
        if self.sample_labels and stem in self.sample_labels:
            return self.sample_labels[stem]
        return stem


# Canonical variant coordinate used as the merge key across samples
VariantKey = tuple[str, int, str, str]  # (chrom, pos, ref, alt)


@dataclass(frozen=True, slots=True)
class SampleOccurrence:
    """Record of a variant's appearance in a single sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier (file stem or label).
    vcf_path : Path
        Path to the source VCF file.
    classified : ClassifiedVariant
        Full classification result from the standard pipeline.
    """

    sample_id: str
    vcf_path: Path
    classified: ClassifiedVariant


@dataclass(frozen=True, slots=True)
class CohortVariant:
    """A variant aggregated across multiple samples in a cohort.

    Parameters
    ----------
    chrom : str
        Chromosome.
    pos : int
        1-based genomic position.
    ref : str
        Reference allele.
    alt : str
        Alternate allele.
    gene_name : str | None
        Gene symbol if annotated.
    consequence : FunctionalConsequence
        Most severe consequence observed across samples.
    sample_count : int
        Number of samples carrying this variant.
    total_samples : int
        Total samples in the cohort (denominator for frequency).
    occurrences : tuple[SampleOccurrence, ...]
        Per-sample classification details.
    max_classification : ACMGClassification
        Most severe ACMG classification across all samples.
    all_evidence_tags : frozenset[EvidenceTag]
        Union of all evidence tags assigned across samples.
    allele_frequency : float | None
        Population allele frequency from gnomAD (consensus across samples).
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene_name: Optional[str]
    consequence: FunctionalConsequence
    sample_count: int
    total_samples: int
    occurrences: tuple[SampleOccurrence, ...]
    max_classification: ACMGClassification
    all_evidence_tags: frozenset[EvidenceTag]
    allele_frequency: Optional[float] = None

    @property
    def key(self) -> VariantKey:
        """Canonical coordinate key for this variant."""
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def cohort_frequency(self) -> float:
        """Fraction of cohort samples carrying this variant."""
        if self.total_samples == 0:
            return 0.0
        return self.sample_count / self.total_samples

    @property
    def is_singleton(self) -> bool:
        """True if variant appears in exactly one sample."""
        return self.sample_count == 1

    @property
    def is_universal(self) -> bool:
        """True if variant appears in all cohort samples."""
        return self.sample_count == self.total_samples

    @property
    def sample_ids(self) -> list[str]:
        """List of sample identifiers carrying this variant."""
        return [occ.sample_id for occ in self.occurrences]


@dataclass(frozen=True, slots=True)
class GeneBurden:
    """Per-gene variant burden across the cohort.

    Parameters
    ----------
    gene_name : str
        Gene symbol.
    total_variants : int
        Total distinct variants in this gene across all samples.
    pathogenic_count : int
        Variants classified as Pathogenic or Likely_Pathogenic.
    samples_affected : int
        Number of samples with at least one variant in this gene.
    total_samples : int
        Total samples in cohort.
    most_severe : ACMGClassification
        Most severe classification in this gene.
    """

    gene_name: str
    total_variants: int
    pathogenic_count: int
    samples_affected: int
    total_samples: int
    most_severe: ACMGClassification

    @property
    def penetrance(self) -> float:
        """Fraction of cohort samples affected in this gene."""
        if self.total_samples == 0:
            return 0.0
        return self.samples_affected / self.total_samples


@dataclass(frozen=True, slots=True)
class CohortSummary:
    """Aggregate statistics for a completed cohort analysis.

    Parameters
    ----------
    cohort_name : str
        Cohort identifier from config.
    total_samples : int
        Number of samples analyzed.
    total_variants : int
        Total distinct variants across all samples.
    shared_variants : int
        Variants appearing in >= 2 samples.
    singleton_variants : int
        Variants appearing in exactly 1 sample.
    universal_variants : int
        Variants appearing in all samples.
    pathogenic_variants : int
        Variants classified Pathogenic in at least one sample.
    likely_pathogenic_variants : int
        Variants classified Likely_Pathogenic in at least one sample.
    genes_affected : int
        Total distinct genes with at least one variant.
    top_recurrent_genes : tuple[str, ...]
        Gene symbols with highest recurrence (top 10).
    samples_processed : tuple[str, ...]
        Sample identifiers in processing order.
    """

    cohort_name: str
    total_samples: int
    total_variants: int
    shared_variants: int
    singleton_variants: int
    universal_variants: int
    pathogenic_variants: int
    likely_pathogenic_variants: int
    genes_affected: int
    top_recurrent_genes: tuple[str, ...]
    samples_processed: tuple[str, ...]
