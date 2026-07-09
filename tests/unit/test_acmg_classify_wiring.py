"""Unit tests for ACMGClassifier.classify combining rules integration.

Validates that the classify method correctly wires evidence tag assignment
with ACMG/AMP 2015 combining rules to produce the correct final classification.
"""

from __future__ import annotations

from vartriage.classification.acmg import ACMGClassifier
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)


def _make_scored_variant(
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    allele_frequency: float | None = 0.005,
    clinvar_assertion: ClinVarAssertion | None = None,
    frequency_unknown: bool = False,
    clinvar_unknown: bool = False,
    revel_score: float | None = 0.5,
    cadd_phred: float | None = 25.0,
) -> ScoredVariant:
    """Create a ScoredVariant with configurable fields for testing."""
    v = Variant(
        chrom="chr1",
        pos=100,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=v,
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
        frequency_unknown=frequency_unknown,
        clinvar_unknown=clinvar_unknown,
    )
    cadd_normalized = None
    if cadd_phred is not None:
        cadd_normalized = min(cadd_phred / 99.0, 1.0)
    return ScoredVariant(
        annotated=annotated,
        cadd_phred=cadd_phred,
        cadd_normalized=cadd_normalized,
        revel_score=revel_score,
        composite_rank=None,
    )


class TestClassifyCombiningWiring:
    """Verify classify method applies combining rules correctly."""

    def test_pvs1_pp3_pp5_yields_pathogenic(self) -> None:
        """PVS1 (Very Strong) + PP3 + PP5 (2 Supporting) -> Pathogenic.

        Combining rule: >=1 Very Strong AND >=2 Supporting.
        """
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.NONSENSE,
            revel_score=0.85,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))

        assert len(results) == 1
        result = results[0]
        assert EvidenceTag.PVS1 in result.evidence_tags
        assert EvidenceTag.PP3 in result.evidence_tags
        assert EvidenceTag.PP5 in result.evidence_tags
        assert result.classification == ACMGClassification.PATHOGENIC

    def test_pvs1_pm2_yields_likely_pathogenic(self) -> None:
        """PVS1 (Very Strong) + PM2 (Moderate) -> Likely_Pathogenic.

        Combining rule: 1 Very Strong AND 1 Moderate.
        """
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.FRAMESHIFT,
            allele_frequency=0.00005,
            revel_score=0.3,
            clinvar_assertion=None,
            clinvar_unknown=True,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))

        assert len(results) == 1
        result = results[0]
        assert EvidenceTag.PVS1 in result.evidence_tags
        assert EvidenceTag.PM2 in result.evidence_tags
        assert EvidenceTag.PP3 not in result.evidence_tags
        assert result.classification == ACMGClassification.LIKELY_PATHOGENIC

    def test_no_tags_yields_vus(self) -> None:
        """A variant with no evidence tags gets classified as VUS."""
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.05,
            revel_score=0.3,
            clinvar_assertion=ClinVarAssertion.VUS,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))

        assert len(results) == 1
        result = results[0]
        assert len(result.evidence_tags) == 0
        assert result.classification == ACMGClassification.VUS

    def test_single_pvs1_yields_vus(self) -> None:
        """PVS1 alone does not meet any combining rule threshold."""
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.NONSENSE,
            allele_frequency=0.01,
            revel_score=0.5,
            clinvar_assertion=ClinVarAssertion.VUS,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))

        assert len(results) == 1
        result = results[0]
        assert EvidenceTag.PVS1 in result.evidence_tags
        assert result.classification == ACMGClassification.VUS

    def test_missing_data_sources_populated(self) -> None:
        """Missing data sources are correctly tracked alongside classification."""
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.NONSENSE,
            allele_frequency=None,
            frequency_unknown=True,
            revel_score=None,
            clinvar_assertion=None,
            clinvar_unknown=True,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))

        result = results[0]
        assert EvidenceTag.PVS1 in result.evidence_tags
        assert "gnomAD" in result.missing_data_sources
        assert "REVEL" in result.missing_data_sources
        assert "ClinVar" in result.missing_data_sources
        assert result.classification == ACMGClassification.VUS

    def test_all_tags_assigned_yields_pathogenic(self) -> None:
        """All four tags: PVS1+PM2+PP3+PP5 -> Pathogenic."""
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.FRAMESHIFT,
            allele_frequency=0.00001,
            revel_score=0.9,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))

        result = results[0]
        assert EvidenceTag.PVS1 in result.evidence_tags
        assert EvidenceTag.PM2 in result.evidence_tags
        assert EvidenceTag.PP3 in result.evidence_tags
        assert EvidenceTag.PP5 in result.evidence_tags
        assert result.classification == ACMGClassification.PATHOGENIC
