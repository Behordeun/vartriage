"""Data models for variant representation, configuration, and warnings."""

from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    ClinVarAssertion,
    CONSEQUENCE_SEVERITY_ORDER,
    EVIDENCE_STRENGTH_MAP,
    EvidenceStrength,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)

__all__ = [
    "ACMGClassification",
    "AnnotatedVariant",
    "ClassifiedVariant",
    "ClinVarAssertion",
    "CONSEQUENCE_SEVERITY_ORDER",
    "EVIDENCE_STRENGTH_MAP",
    "EvidenceStrength",
    "EvidenceTag",
    "FunctionalConsequence",
    "ScoredVariant",
    "Variant",
]
