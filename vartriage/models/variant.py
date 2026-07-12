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

    PVS1 = "PVS1"
    PM2 = "PM2"
    PP3 = "PP3"
    PP5 = "PP5"


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

    VERY_STRONG = "Very_Strong"
    STRONG = "Strong"
    MODERATE = "Moderate"
    SUPPORTING = "Supporting"


EVIDENCE_STRENGTH_MAP: dict[EvidenceTag, EvidenceStrength] = {
    EvidenceTag.PVS1: EvidenceStrength.VERY_STRONG,
    EvidenceTag.PM2: EvidenceStrength.MODERATE,
    EvidenceTag.PP3: EvidenceStrength.SUPPORTING,
    EvidenceTag.PP5: EvidenceStrength.SUPPORTING,
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
