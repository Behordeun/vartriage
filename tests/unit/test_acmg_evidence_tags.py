"""Unit tests for ACMG evidence tag assignment."""

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
    cadd_normalized: float | None = None,
    composite_rank: float | None = None,
) -> ScoredVariant:
    """Helper to create a ScoredVariant with configurable fields."""
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
    if cadd_normalized is None and cadd_phred is not None:
        cadd_normalized = min(cadd_phred / 99.0, 1.0)
    return ScoredVariant(
        annotated=annotated,
        cadd_phred=cadd_phred,
        cadd_normalized=cadd_normalized,
        revel_score=revel_score,
        composite_rank=composite_rank,
    )


class TestPVS1Assignment:
    """PVS1 is assigned when consequence is Nonsense or Frameshift."""

    def test_assigns_pvs1_for_nonsense(self) -> None:
        sv = _make_scored_variant(consequence=FunctionalConsequence.NONSENSE)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PVS1 in results[0].evidence_tags

    def test_assigns_pvs1_for_frameshift(self) -> None:
        sv = _make_scored_variant(consequence=FunctionalConsequence.FRAMESHIFT)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PVS1 in results[0].evidence_tags

    def test_does_not_assign_pvs1_for_missense(self) -> None:
        sv = _make_scored_variant(consequence=FunctionalConsequence.MISSENSE)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PVS1 not in results[0].evidence_tags

    def test_does_not_assign_pvs1_for_synonymous(self) -> None:
        sv = _make_scored_variant(consequence=FunctionalConsequence.SYNONYMOUS)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PVS1 not in results[0].evidence_tags

    def test_does_not_assign_pvs1_for_splice_site(self) -> None:
        sv = _make_scored_variant(consequence=FunctionalConsequence.SPLICE_SITE)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PVS1 not in results[0].evidence_tags


class TestPM2Assignment:
    """PM2 is assigned when allele frequency < 0.0001."""

    def test_assigns_pm2_for_very_rare_variant(self) -> None:
        sv = _make_scored_variant(allele_frequency=0.00005)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PM2 in results[0].evidence_tags

    def test_assigns_pm2_for_zero_frequency(self) -> None:
        sv = _make_scored_variant(allele_frequency=0.0)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PM2 in results[0].evidence_tags

    def test_does_not_assign_pm2_at_threshold(self) -> None:
        sv = _make_scored_variant(allele_frequency=0.0001)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PM2 not in results[0].evidence_tags

    def test_does_not_assign_pm2_above_threshold(self) -> None:
        sv = _make_scored_variant(allele_frequency=0.01)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PM2 not in results[0].evidence_tags

    def test_omits_pm2_when_frequency_unavailable(self) -> None:
        sv = _make_scored_variant(allele_frequency=None, frequency_unknown=True)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PM2 not in results[0].evidence_tags
        assert "gnomAD" in results[0].missing_data_sources


class TestPP3Assignment:
    """PP3 is assigned when REVEL > 0.7."""

    def test_assigns_pp3_for_high_revel(self) -> None:
        sv = _make_scored_variant(revel_score=0.85)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP3 in results[0].evidence_tags

    def test_does_not_assign_pp3_at_threshold(self) -> None:
        sv = _make_scored_variant(revel_score=0.7)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP3 not in results[0].evidence_tags

    def test_does_not_assign_pp3_below_threshold(self) -> None:
        sv = _make_scored_variant(revel_score=0.5)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP3 not in results[0].evidence_tags

    def test_omits_pp3_when_revel_unavailable(self) -> None:
        sv = _make_scored_variant(revel_score=None)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP3 not in results[0].evidence_tags
        assert "REVEL" in results[0].missing_data_sources


class TestPP5Assignment:
    """PP5 is assigned for ClinVar Pathogenic with no conflicting assertion."""

    def test_assigns_pp5_for_pathogenic(self) -> None:
        sv = _make_scored_variant(clinvar_assertion=ClinVarAssertion.PATHOGENIC)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP5 in results[0].evidence_tags

    def test_does_not_assign_pp5_for_benign(self) -> None:
        sv = _make_scored_variant(clinvar_assertion=ClinVarAssertion.BENIGN)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP5 not in results[0].evidence_tags

    def test_does_not_assign_pp5_for_likely_benign(self) -> None:
        sv = _make_scored_variant(clinvar_assertion=ClinVarAssertion.LIKELY_BENIGN)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP5 not in results[0].evidence_tags

    def test_does_not_assign_pp5_for_vus(self) -> None:
        sv = _make_scored_variant(clinvar_assertion=ClinVarAssertion.VUS)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP5 not in results[0].evidence_tags

    def test_does_not_assign_pp5_for_likely_pathogenic(self) -> None:
        sv = _make_scored_variant(
            clinvar_assertion=ClinVarAssertion.LIKELY_PATHOGENIC
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP5 not in results[0].evidence_tags

    def test_omits_pp5_when_clinvar_unavailable(self) -> None:
        sv = _make_scored_variant(clinvar_assertion=None, clinvar_unknown=True)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP5 not in results[0].evidence_tags
        assert "ClinVar" in results[0].missing_data_sources


class TestMissingDataSources:
    """Missing data sources are tracked when data is unavailable."""

    def test_records_all_missing_sources(self) -> None:
        sv = _make_scored_variant(
            allele_frequency=None,
            frequency_unknown=True,
            clinvar_assertion=None,
            clinvar_unknown=True,
            revel_score=None,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        missing = results[0].missing_data_sources
        assert "gnomAD" in missing
        assert "ClinVar" in missing
        assert "REVEL" in missing

    def test_no_missing_sources_when_all_available(self) -> None:
        sv = _make_scored_variant(
            allele_frequency=0.0005,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC,
            revel_score=0.85,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert len(results[0].missing_data_sources) == 0


class TestClassifyOutput:
    """Overall classify method output validation."""

    def test_classification_defaults_to_vus(self) -> None:
        sv = _make_scored_variant()
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert results[0].classification == ACMGClassification.VUS

    def test_processes_empty_iterator(self) -> None:
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([])))
        assert results == []

    def test_processes_multiple_variants(self) -> None:
        variants = [
            _make_scored_variant(consequence=FunctionalConsequence.NONSENSE),
            _make_scored_variant(consequence=FunctionalConsequence.MISSENSE),
            _make_scored_variant(
                consequence=FunctionalConsequence.FRAMESHIFT,
                allele_frequency=0.00001,
                revel_score=0.9,
                clinvar_assertion=ClinVarAssertion.PATHOGENIC,
            ),
        ]
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter(variants)))
        assert len(results) == 3

        # First variant: nonsense → PVS1
        assert EvidenceTag.PVS1 in results[0].evidence_tags

        # Second variant: missense, no special tags except based on AF/scores
        assert EvidenceTag.PVS1 not in results[1].evidence_tags

        # Third variant: frameshift + rare + high REVEL + ClinVar pathogenic
        tags = results[2].evidence_tags
        assert EvidenceTag.PVS1 in tags
        assert EvidenceTag.PM2 in tags
        assert EvidenceTag.PP3 in tags
        assert EvidenceTag.PP5 in tags

    def test_scored_variant_preserved_in_output(self) -> None:
        sv = _make_scored_variant(revel_score=0.9)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert results[0].scored is sv
