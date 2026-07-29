"""Data models for structural variant representation and classification.

Frozen dataclasses and enums representing SVs from VCF parsing through
ClinGen-based interpretation. Separate from the point-variant models
because SVs carry fundamentally different attributes (spans, breakpoints,
copy number state) and use a different classification framework.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional


class SVType(Enum):
    """Structural variant type from the VCF SVTYPE field.

    Covers the standard SV types defined in VCF 4.3+ specification.
    """

    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    BND = "BND"
    CNV = "CNV"


class SVConsequence(Enum):
    """Predicted gene-level impact of a structural variant.

    Ordered by clinical severity (highest first). When an SV overlaps
    multiple genes with different impact levels, the pipeline assigns
    the most severe value.
    """

    WHOLE_GENE_DELETION = "Whole_Gene_Deletion"
    PARTIAL_GENE_DELETION = "Partial_Gene_Deletion"
    WHOLE_GENE_DUPLICATION = "Whole_Gene_Duplication"
    PARTIAL_GENE_DUPLICATION = "Partial_Gene_Duplication"
    GENE_DISRUPTION = "Gene_Disruption"
    INTRONIC = "Intronic"
    REGULATORY = "Regulatory"
    INTERGENIC = "Intergenic"


# Severity ordering for SVConsequence: index 0 is most severe.
SV_CONSEQUENCE_SEVERITY: list[SVConsequence] = [
    SVConsequence.WHOLE_GENE_DELETION,
    SVConsequence.PARTIAL_GENE_DELETION,
    SVConsequence.WHOLE_GENE_DUPLICATION,
    SVConsequence.PARTIAL_GENE_DUPLICATION,
    SVConsequence.GENE_DISRUPTION,
    SVConsequence.INTRONIC,
    SVConsequence.REGULATORY,
    SVConsequence.INTERGENIC,
]


class SVClassification(Enum):
    """Final ClinGen-based classification for a structural variant.

    Uses the 5-tier system from the ACMG/ClinGen Technical Standards
    for interpretation of copy-number gains and losses (2020).
    """

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_Pathogenic"
    VUS = "VUS"
    LIKELY_BENIGN = "Likely_Benign"
    BENIGN = "Benign"


class SVEvidenceCategory(Enum):
    """ClinGen SV evidence categories for scoring.

    Based on Riggs et al. 2020 framework. Categories 1-3 cover
    losses (deletions), categories 4+ cover gains (duplications).
    """

    # Section 1: Initial assessment of genomic content
    CONTAINS_PROTEIN_CODING = "1A"
    CONTAINS_ESTABLISHED_HI_GENE = "1B"

    # Section 2: Overlap with established pathogenic/benign regions
    COMPLETE_OVERLAP_PATHOGENIC = "2A"
    PARTIAL_OVERLAP_PATHOGENIC = "2B"
    OVERLAP_SMALLER_THAN_PATHOGENIC = "2C"
    CONTAINED_WITHIN_BENIGN = "2D"
    PARTIAL_OVERLAP_BENIGN = "2E"
    COMPLETELY_CONTAINS_BENIGN = "2F"
    NO_OVERLAP_KNOWN = "2G"
    HI_GENE_COUNT = "2H"

    # Section 3: Gene-level evaluation
    GENE_FULLY_CONTAINED = "3A"
    GENE_PARTIALLY_DELETED = "3B"
    BREAKPOINT_WITHIN_GENE = "3C"

    # Section 4: Duplication-specific evaluation
    DUP_COMPLETE_OVERLAP_PATHOGENIC = "4A"
    DUP_IDENTICAL_TO_PATHOGENIC = "4B"
    DUP_SMALLER_THAN_PATHOGENIC = "4C"
    DUP_CONTAINED_WITHIN_BENIGN = "4D"
    DUP_NO_OVERLAP = "4E"
    DUP_TS_GENE_CONTAINED = "4F"
    DUP_GENE_DISRUPTED = "4G"
    DUP_INTRAGENIC_NO_DISRUPTION = "4H"

    # Section 5: Family and clinical evaluation
    SEGREGATION_EVIDENCE = "5A"
    DE_NOVO = "5B"
    PHENOTYPE_CONSISTENT = "5C"


# Point values for each evidence category (ClinGen 2020 scoring)
SV_EVIDENCE_POINTS: dict[SVEvidenceCategory, float] = {
    SVEvidenceCategory.CONTAINS_PROTEIN_CODING: 0.0,
    SVEvidenceCategory.CONTAINS_ESTABLISHED_HI_GENE: 0.0,
    SVEvidenceCategory.COMPLETE_OVERLAP_PATHOGENIC: 1.0,
    SVEvidenceCategory.PARTIAL_OVERLAP_PATHOGENIC: 0.5,
    SVEvidenceCategory.OVERLAP_SMALLER_THAN_PATHOGENIC: 0.0,
    SVEvidenceCategory.CONTAINED_WITHIN_BENIGN: -1.0,
    SVEvidenceCategory.PARTIAL_OVERLAP_BENIGN: -0.5,
    SVEvidenceCategory.COMPLETELY_CONTAINS_BENIGN: 0.0,
    SVEvidenceCategory.NO_OVERLAP_KNOWN: 0.0,
    SVEvidenceCategory.HI_GENE_COUNT: 0.0,
    SVEvidenceCategory.GENE_FULLY_CONTAINED: 0.5,
    SVEvidenceCategory.GENE_PARTIALLY_DELETED: 0.25,
    SVEvidenceCategory.BREAKPOINT_WITHIN_GENE: 0.25,
    SVEvidenceCategory.DUP_COMPLETE_OVERLAP_PATHOGENIC: 1.0,
    SVEvidenceCategory.DUP_IDENTICAL_TO_PATHOGENIC: 0.75,
    SVEvidenceCategory.DUP_SMALLER_THAN_PATHOGENIC: 0.0,
    SVEvidenceCategory.DUP_CONTAINED_WITHIN_BENIGN: -1.0,
    SVEvidenceCategory.DUP_NO_OVERLAP: 0.0,
    SVEvidenceCategory.DUP_TS_GENE_CONTAINED: 0.5,
    SVEvidenceCategory.DUP_GENE_DISRUPTED: 0.5,
    SVEvidenceCategory.DUP_INTRAGENIC_NO_DISRUPTION: 0.0,
    SVEvidenceCategory.SEGREGATION_EVIDENCE: 0.15,
    SVEvidenceCategory.DE_NOVO: 0.3,
    SVEvidenceCategory.PHENOTYPE_CONSISTENT: 0.15,
}


@dataclass(frozen=True, slots=True)
class Breakpoint:
    """One end of a structural variant breakpoint pair.

    For BND (translocation) records, both breakpoints define the
    rearrangement. For simple SVs (DEL/DUP/INV), the start breakpoint
    is the POS and the end is END/POS+SVLEN.
    """

    chrom: str
    pos: int
    confidence_interval: tuple[int, int] = (0, 0)


@dataclass(frozen=True, slots=True)
class StructuralVariant:
    """Parsed structural variant record from VCF.

    Captures SV-specific attributes beyond what the point-variant
    Variant dataclass holds. The original VCF line fields (CHROM, POS,
    etc.) are stored directly rather than wrapping a Variant, since
    REF/ALT semantics differ for SVs (ALT is symbolic like <DEL>).

    Parameters
    ----------
    chrom : str
        Chromosome of the start breakpoint.
    start : int
        1-based start position (VCF POS field).
    end : int
        1-based end position (inclusive). From INFO/END or POS+|SVLEN|.
    sv_type : SVType
        Structural variant class.
    id : Optional[str]
        Variant ID from the VCF ID column.
    svlen : Optional[int]
        SV length. Negative for deletions per VCF spec.
    qual : Optional[float]
        Phred-scaled quality.
    filter_status : str
        VCF FILTER field value.
    alt : str
        ALT field (symbolic allele like <DEL> or BND notation).
    copy_number : Optional[int]
        Copy number state from CN INFO field.
    start_ci : tuple[int, int]
        Confidence interval around start position (CIPOS).
    end_ci : tuple[int, int]
        Confidence interval around end position (CIEND).
    mate_id : Optional[str]
        MATEID for BND records linking paired breakends.
    info : dict[str, Any]
        Raw INFO fields for pass-through to reports.
    """

    chrom: str
    start: int
    end: int
    sv_type: SVType
    id: Optional[str] = None
    svlen: Optional[int] = None
    qual: Optional[float] = None
    filter_status: str = "PASS"
    alt: str = ""
    copy_number: Optional[int] = None
    start_ci: tuple[int, int] = (0, 0)
    end_ci: tuple[int, int] = (0, 0)
    mate_id: Optional[str] = None
    info: dict[str, Any] = field(default_factory=dict)

    @property
    def length(self) -> int:
        """Span of the SV in base pairs."""
        if self.svlen is not None:
            return abs(self.svlen)
        return self.end - self.start + 1

    @property
    def is_intrachromosomal(self) -> bool:
        """True for SVs contained within a single chromosome."""
        return self.sv_type != SVType.BND


@dataclass(frozen=True, slots=True)
class GeneOverlap:
    """Result of overlapping an SV with a single gene.

    Captures the relationship between the SV boundaries and the
    gene boundaries to determine functional impact.
    """

    gene_symbol: str
    gene_chrom: str
    gene_start: int
    gene_end: int
    overlap_fraction: float
    is_whole_gene: bool
    exons_affected: int
    total_exons: int
    is_haploinsufficient: bool = False
    is_triplosensitive: bool = False
    hi_score: Optional[float] = None
    ts_score: Optional[float] = None


@dataclass(frozen=True, slots=True)
class AnnotatedSV:
    """Structural variant enriched with gene overlap and frequency data.

    Parameters
    ----------
    sv : StructuralVariant
        Original parsed SV record.
    consequence : SVConsequence
        Most severe gene-level impact.
    gene_overlaps : tuple[GeneOverlap, ...]
        Per-gene overlap details, sorted by severity.
    population_frequency : Optional[float]
        gnomAD-SV frequency based on reciprocal overlap matching.
    frequency_unknown : bool
        True if no matching SV found in reference database.
    genes_affected : int
        Total number of protein-coding genes overlapped.
    hi_genes_affected : int
        Number of haploinsufficient genes overlapped.
    """

    sv: StructuralVariant
    consequence: SVConsequence
    gene_overlaps: tuple[GeneOverlap, ...] = ()
    population_frequency: Optional[float] = None
    frequency_unknown: bool = True
    genes_affected: int = 0
    hi_genes_affected: int = 0


@dataclass(frozen=True, slots=True)
class ScoredSV:
    """Annotated SV with pathogenicity score computed.

    The pathogenicity score integrates:
    - Gene impact severity
    - Dosage sensitivity (HI/TS scores)
    - Population frequency rarity
    - SV size relative to gene content

    Parameters
    ----------
    annotated : AnnotatedSV
        Fully annotated SV.
    pathogenicity_score : Optional[float]
        Composite pathogenicity score (0.0-1.0). Higher means more
        likely pathogenic. None if scoring is impossible (intergenic
        SV with no frequency data).
    dosage_score : Optional[float]
        Dosage sensitivity contribution (0.0-1.0).
    size_score : float
        Size-based score component.
    frequency_score : float
        Rarity score (1.0 for absent, 0.0 for common).
    """

    annotated: AnnotatedSV
    pathogenicity_score: Optional[float] = None
    dosage_score: Optional[float] = None
    size_score: float = 0.0
    frequency_score: float = 1.0


@dataclass(frozen=True, slots=True)
class ClassifiedSV:
    """Scored SV with ClinGen-based classification.

    Parameters
    ----------
    scored : ScoredSV
        Scored SV with pathogenicity metrics.
    classification : SVClassification
        Final 5-tier classification.
    evidence_categories : frozenset[SVEvidenceCategory]
        Evidence categories satisfied during evaluation.
    evidence_score : float
        Accumulated point total from evidence categories.
    missing_data_sources : frozenset[str]
        Data sources unavailable during classification.
    syndrome_name : Optional[str]
        Name of the matching known syndrome (e.g., "22q11.2 deletion
        syndrome") when the SV overlaps a curated pathogenic region.
    """

    scored: ScoredSV
    classification: SVClassification = SVClassification.VUS
    evidence_categories: frozenset[SVEvidenceCategory] = field(
        default_factory=frozenset
    )
    evidence_score: float = 0.0
    missing_data_sources: frozenset[str] = field(default_factory=frozenset)
    syndrome_name: Optional[str] = None
