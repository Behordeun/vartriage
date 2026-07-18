"""Core data models for variant representation and classification.

Immutable dataclasses and enums used from VCF parsing through
ACMG/AMP classification.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional


class FunctionalConsequence(Enum):
    """Predicted biological effect of a variant on gene function.

    Members are ordered by severity (highest first). When a variant overlaps
    multiple transcripts with different consequences, the pipeline assigns the
    most severe value according to this ordering.

    Attributes
    ----------
    FRAMESHIFT : str
        Insertion or deletion that shifts the reading frame.
    NONSENSE : str
        Premature stop codon introduction.
    SPLICE_SITE : str
        Variant within 2 bases of an exon-intron junction.
    MISSENSE : str
        Single amino acid substitution.
    IN_FRAME_INSERTION : str
        Insertion that preserves the reading frame.
    IN_FRAME_DELETION : str
        Deletion that preserves the reading frame.
    SYNONYMOUS : str
        Codon change with no amino acid change.
    INTERGENIC : str
        Variant falls outside any known coding region.
    """

    FRAMESHIFT = "Frameshift"
    NONSENSE = "Nonsense"
    SPLICE_SITE = "Splice_Site"
    MISSENSE = "Missense"
    IN_FRAME_INSERTION = "In_Frame_Insertion"
    IN_FRAME_DELETION = "In_Frame_Deletion"
    SYNONYMOUS = "Synonymous"
    INTERGENIC = "Intergenic"


# Severity ordering: index 0 is most severe.
CONSEQUENCE_SEVERITY_ORDER: list[FunctionalConsequence] = [
    FunctionalConsequence.FRAMESHIFT,
    FunctionalConsequence.NONSENSE,
    FunctionalConsequence.SPLICE_SITE,
    FunctionalConsequence.MISSENSE,
    FunctionalConsequence.IN_FRAME_INSERTION,
    FunctionalConsequence.IN_FRAME_DELETION,
    FunctionalConsequence.SYNONYMOUS,
    FunctionalConsequence.INTERGENIC,
]


class ClinVarAssertion(Enum):
    """ClinVar clinical significance categories.

    Attributes
    ----------
    PATHOGENIC : str
        Variant is pathogenic.
    LIKELY_PATHOGENIC : str
        Variant is likely pathogenic.
    VUS : str
        Variant of uncertain significance.
    LIKELY_BENIGN : str
        Variant is likely benign.
    BENIGN : str
        Variant is benign.
    """

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_Pathogenic"
    VUS = "VUS"
    LIKELY_BENIGN = "Likely_Benign"
    BENIGN = "Benign"


class ACMGClassification(Enum):
    """Final ACMG/AMP 2015 classification for a variant.

    In v0.x, only PATHOGENIC, LIKELY_PATHOGENIC, and VUS are produced.
    LIKELY_BENIGN and BENIGN exist for forward compatibility but the
    current rules don't assign benign evidence tags yet.
    """

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_Pathogenic"
    VUS = "VUS"
    LIKELY_BENIGN = "Likely_Benign"
    BENIGN = "Benign"


class EvidenceTag(Enum):
    """ACMG/AMP evidence tags assigned during classification.

    Each tag corresponds to a specific criterion from the ACMG/AMP 2015
    guidelines, with an associated strength tier defined in
    ``EVIDENCE_STRENGTH_MAP``.

    Attributes
    ----------
    PVS1 : str
        Very Strong evidence: null variant (nonsense or frameshift).
    PM2 : str
        Moderate evidence: absent from population controls.
    PP3 : str
        Supporting evidence: computational pathogenicity prediction.
    PP5 : str
        Supporting evidence: reputable clinical source (ClinVar).
    """

    # Pathogenic evidence
    PVS1 = "PVS1"
    PS1 = "PS1"
    PM1 = "PM1"
    PM2 = "PM2"
    PM4 = "PM4"
    PM5 = "PM5"
    PP3 = "PP3"
    PP5 = "PP5"

    # Benign evidence
    BA1 = "BA1"
    BS1 = "BS1"
    BS2 = "BS2"
    BP4 = "BP4"
    BP7 = "BP7"


class EvidenceStrength(Enum):
    """ACMG evidence strength tiers used in combining rules.

    Attributes
    ----------
    VERY_STRONG : str
        Very strong level of evidence.
    STRONG : str
        Strong level of evidence.
    MODERATE : str
        Moderate level of evidence.
    SUPPORTING : str
        Supporting level of evidence.
    """

    STANDALONE = "Standalone"
    VERY_STRONG = "Very_Strong"
    STRONG = "Strong"
    MODERATE = "Moderate"
    SUPPORTING = "Supporting"


EVIDENCE_STRENGTH_MAP: dict[EvidenceTag, EvidenceStrength] = {
    # Pathogenic evidence
    EvidenceTag.PVS1: EvidenceStrength.VERY_STRONG,
    EvidenceTag.PS1: EvidenceStrength.STRONG,
    EvidenceTag.PM1: EvidenceStrength.MODERATE,
    EvidenceTag.PM2: EvidenceStrength.MODERATE,
    EvidenceTag.PM4: EvidenceStrength.MODERATE,
    EvidenceTag.PM5: EvidenceStrength.MODERATE,
    EvidenceTag.PP3: EvidenceStrength.SUPPORTING,
    EvidenceTag.PP5: EvidenceStrength.SUPPORTING,
    # Benign evidence
    EvidenceTag.BA1: EvidenceStrength.STANDALONE,
    EvidenceTag.BS1: EvidenceStrength.STRONG,
    EvidenceTag.BS2: EvidenceStrength.STRONG,
    EvidenceTag.BP4: EvidenceStrength.SUPPORTING,
    EvidenceTag.BP7: EvidenceStrength.SUPPORTING,
}
"""Mapping of evidence tags to their strength tiers.

Used by the ACMG combining rules to determine final classification
based on accumulated evidence.
"""


@dataclass(frozen=True, slots=True)
class Variant:
    """Raw variant record from VCF parsing.

    Parameters
    ----------
    chrom : str
        Chromosome name (e.g., "chr1", "1").
    pos : int
        1-based genomic position.
    id : Optional[str]
        Variant identifier or None if missing.
    ref : str
        Reference allele.
    alt : str
        Alternate allele (single ALT; multiallelic split upstream).
    qual : Optional[float]
        Phred-scaled quality score or None if missing.
    filter_status : str
        FILTER field value (e.g., "PASS", ".", "LowQual").
    info : dict[str, Any]
        INFO field key-value pairs.
    """

    chrom: str
    pos: int
    id: Optional[str]
    ref: str
    alt: str
    qual: Optional[float]
    filter_status: str
    info: dict[str, Any] = field(default_factory=dict)


class Zygosity(Enum):
    """Genotype zygosity state for a variant in a specific sample."""

    HETEROZYGOUS = "Heterozygous"
    HOMOZYGOUS_ALT = "Homozygous_Alt"
    HEMIZYGOUS = "Hemizygous"
    UNKNOWN = "Unknown"


@dataclass(frozen=True, slots=True)
class VariantQualityMetrics:
    """Per-variant sequencing quality metrics from VCF FORMAT fields."""

    depth: Optional[int] = None
    genotype_quality: Optional[int] = None
    allele_balance: Optional[float] = None
    is_low_confidence: bool = False


@dataclass(frozen=True, slots=True)
class PopulationFrequencies:
    """Per-population gnomAD allele frequencies."""

    global_af: Optional[float] = None
    afr: Optional[float] = None
    amr: Optional[float] = None
    asj: Optional[float] = None
    eas: Optional[float] = None
    fin: Optional[float] = None
    nfe: Optional[float] = None
    sas: Optional[float] = None

    @property
    def max_population_af(self) -> Optional[float]:
        """Highest frequency across all population subgroups."""
        values = [
            v
            for v in (
                self.afr,
                self.amr,
                self.asj,
                self.eas,
                self.fin,
                self.nfe,
                self.sas,
            )
            if v is not None
        ]
        if values:
            return max(values)
        return self.global_af

    def any_exceeds(self, threshold: float) -> bool:
        """True if any population-specific AF exceeds the threshold."""
        max_af = self.max_population_af
        return max_af is not None and max_af > threshold

    def all_below(self, threshold: float) -> bool:
        """True if ALL population AFs are strictly below threshold (or None).

        Uses >= to reject: a value AT the threshold is NOT below it.
        This matches the global AF path which uses `af < threshold`.
        """
        for af in (
            self.afr,
            self.amr,
            self.asj,
            self.eas,
            self.fin,
            self.nfe,
            self.sas,
            self.global_af,
        ):
            if af is not None and af >= threshold:
                return False
        return True


@dataclass(frozen=True, slots=True)
class AnnotatedVariant:
    """Variant enriched with functional and population annotations.

    Parameters
    ----------
    variant : Variant
        Original variant record.
    consequence : FunctionalConsequence
        Most severe predicted functional effect.
    allele_frequency : Optional[float]
        gnomAD population frequency (0.0-1.0) or None if unknown.
    clinvar_assertion : Optional[ClinVarAssertion]
        ClinVar clinical significance or None if not found.
    frequency_unknown : bool
        True if variant was absent from gnomAD.
    clinvar_unknown : bool
        True if variant was absent from ClinVar.
    gene_name : Optional[str]
        Gene symbol from consequence annotation, or None for intergenic variants.
    """

    variant: Variant
    consequence: FunctionalConsequence
    allele_frequency: Optional[float] = None
    clinvar_assertion: Optional[ClinVarAssertion] = None
    frequency_unknown: bool = False
    clinvar_unknown: bool = False
    gene_name: Optional[str] = None
    population_frequencies: Optional[PopulationFrequencies] = None
    zygosity: Zygosity = Zygosity.UNKNOWN
    quality_metrics: Optional[VariantQualityMetrics] = None


@dataclass(frozen=True, slots=True)
class ScoredVariant:
    """Annotated variant with pathogenicity scores attached.

    Parameters
    ----------
    annotated : AnnotatedVariant
        Fully annotated variant.
    cadd_phred : Optional[float]
        Raw CADD Phred score (0-99+).
    cadd_normalized : Optional[float]
        CADD normalized to 0.0-1.0 scale.
    revel_score : Optional[float]
        REVEL score (0.0-1.0).
    spliceai_score : Optional[float]
        SpliceAI delta score (0.0-1.0) predicting splice-disrupting
        effects. None when SpliceAI is not configured or the variant
        has no lookup match.
    composite_rank : Optional[float]
        Weighted composite rank derived from available scores using
        dynamic proportional weight redistribution.
    """

    annotated: AnnotatedVariant
    cadd_phred: Optional[float] = None
    cadd_normalized: Optional[float] = None
    revel_score: Optional[float] = None
    spliceai_score: Optional[float] = None
    composite_rank: Optional[float] = None
    prioritization_score: Optional[float] = None

    def __post_init__(self) -> None:
        # Sync: prioritization_score is the canonical field,
        # composite_rank kept for backward compatibility
        if self.prioritization_score is None and self.composite_rank is not None:
            object.__setattr__(self, "prioritization_score", self.composite_rank)
        elif self.composite_rank is None and self.prioritization_score is not None:
            object.__setattr__(self, "composite_rank", self.prioritization_score)


@dataclass(frozen=True, slots=True)
class ClassifiedVariant:
    """Scored variant with ACMG/AMP classification.

    Parameters
    ----------
    scored : ScoredVariant
        Scored variant with pathogenicity rank.
    evidence_tags : frozenset[EvidenceTag]
        Set of ACMG evidence tags assigned.
    classification : ACMGClassification
        Final ACMG/AMP classification result.
    missing_data_sources : frozenset[str]
        Names of data sources unavailable for classification.
    """

    scored: ScoredVariant
    evidence_tags: frozenset[EvidenceTag] = field(default_factory=frozenset)
    classification: ACMGClassification = ACMGClassification.VUS
    missing_data_sources: frozenset[str] = field(default_factory=frozenset)
