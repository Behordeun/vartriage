"""Unit tests for benign evidence criteria (BA1, BS1, BP4, BP7)."""

from __future__ import annotations

import pytest

from vartriage.classification.acmg import ACMGClassifier
from vartriage.models.variant import (ACMGClassification, AnnotatedVariant,
                                      ClinVarAssertion, EvidenceTag,
                                      FunctionalConsequence,
                                      PopulationFrequencies, ScoredVariant,
                                      Variant)


def _make_variant(
    consequence: FunctionalConsequence = FunctionalConsequence.SYNONYMOUS,
    allele_frequency: float | None = 0.001,
    population_frequencies: PopulationFrequencies | None = None,
    revel_score: float | None = None,
    spliceai_score: float | None = None,
    cadd_phred: float | None = None,
    clinvar_assertion: ClinVarAssertion | None = None,
) -> ScoredVariant:
    raw = Variant(
        chrom="chr1",
        pos=100,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=raw,
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
        population_frequencies=population_frequencies,
    )
    return ScoredVariant(
        annotated=annotated,
        revel_score=revel_score,
        spliceai_score=spliceai_score,
        cadd_phred=cadd_phred,
    )


def _classify(variant: ScoredVariant) -> set[EvidenceTag]:
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    return set(results[0].evidence_tags)


class TestBA1:
    """BA1: any population AF > 5%."""

    def test_fires_when_global_af_exceeds_five_percent(self) -> None:
        sv = _make_variant(allele_frequency=0.06)
        tags = _classify(sv)
        assert EvidenceTag.BA1 in tags

    def test_does_not_fire_at_exactly_five_percent(self) -> None:
        sv = _make_variant(allele_frequency=0.05)
        tags = _classify(sv)
        assert EvidenceTag.BA1 not in tags

    def test_fires_from_population_specific_frequency(self) -> None:
        pop = PopulationFrequencies(global_af=0.01, asj=0.08)
        sv = _make_variant(allele_frequency=0.01, population_frequencies=pop)
        tags = _classify(sv)
        assert EvidenceTag.BA1 in tags

    def test_does_not_fire_when_all_populations_below(self) -> None:
        pop = PopulationFrequencies(global_af=0.02, afr=0.03, nfe=0.04)
        sv = _make_variant(allele_frequency=0.02, population_frequencies=pop)
        tags = _classify(sv)
        assert EvidenceTag.BA1 not in tags

    def test_does_not_fire_when_af_is_none(self) -> None:
        sv = _make_variant(allele_frequency=None)
        tags = _classify(sv)
        assert EvidenceTag.BA1 not in tags


class TestBS1:
    """BS1: any population AF > 1%."""

    def test_fires_when_global_af_exceeds_one_percent(self) -> None:
        sv = _make_variant(allele_frequency=0.02)
        tags = _classify(sv)
        assert EvidenceTag.BS1 in tags

    def test_does_not_fire_below_one_percent(self) -> None:
        sv = _make_variant(allele_frequency=0.005)
        tags = _classify(sv)
        assert EvidenceTag.BS1 not in tags

    def test_does_not_fire_when_ba1_already_present(self) -> None:
        # AF=6% fires BA1, BS1 should not also fire
        sv = _make_variant(allele_frequency=0.06)
        tags = _classify(sv)
        assert EvidenceTag.BA1 in tags
        assert EvidenceTag.BS1 not in tags

    def test_fires_from_population_specific_frequency(self) -> None:
        pop = PopulationFrequencies(global_af=0.005, fin=0.015)
        sv = _make_variant(allele_frequency=0.005, population_frequencies=pop)
        tags = _classify(sv)
        assert EvidenceTag.BS1 in tags


class TestBP4:
    """BP4: computational benign evidence."""

    def test_fires_for_missense_with_low_revel(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            allele_frequency=0.005,
            revel_score=0.10,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP4 in tags

    def test_does_not_fire_for_missense_with_high_revel(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            allele_frequency=0.005,
            revel_score=0.50,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP4 not in tags

    def test_fires_for_non_missense_with_low_cadd(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.005,
            cadd_phred=5.0,
            spliceai_score=0.5,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP4 in tags

    def test_does_not_fire_for_non_missense_with_high_cadd(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.005,
            cadd_phred=25.0,
            spliceai_score=0.5,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP4 not in tags

    def test_does_not_fire_when_revel_is_none_for_missense(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            allele_frequency=0.005,
            revel_score=None,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP4 not in tags


class TestBP7:
    """BP7: synonymous + no splice impact."""

    def test_fires_for_synonymous_with_low_spliceai(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.005,
            spliceai_score=0.02,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP7 in tags

    def test_does_not_fire_for_synonymous_with_high_spliceai(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.005,
            spliceai_score=0.5,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP7 not in tags

    def test_does_not_fire_for_missense(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            allele_frequency=0.005,
            spliceai_score=0.02,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP7 not in tags

    def test_does_not_fire_when_spliceai_is_none(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.005,
            spliceai_score=None,
        )
        tags = _classify(sv)
        assert EvidenceTag.BP7 not in tags


class TestCombiningWithBenign:
    """Combining rules produce correct classification for benign evidence."""

    def test_ba1_alone_yields_benign(self) -> None:
        sv = _make_variant(allele_frequency=0.10)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert results[0].classification == ACMGClassification.BENIGN

    def test_bs1_plus_bp4_yields_likely_benign(self) -> None:
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            allele_frequency=0.02,
            cadd_phred=5.0,
            spliceai_score=0.5,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        # BS1 (AF > 1%) + BP4 (CADD < 10) = Likely Benign
        assert results[0].classification == ACMGClassification.LIKELY_BENIGN

    def test_conflicting_pathogenic_and_benign_yields_vus(self) -> None:
        # Frameshift (PVS1) + high AF (BA1) = conflicting = VUS
        sv = _make_variant(
            consequence=FunctionalConsequence.FRAMESHIFT,
            allele_frequency=0.10,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert results[0].classification == ACMGClassification.VUS
        assert EvidenceTag.PVS1 in results[0].evidence_tags
        assert EvidenceTag.BA1 in results[0].evidence_tags


class TestPopulationAwarePM2:
    """PM2 uses population-specific thresholds."""

    def test_fires_when_all_populations_below_threshold(self) -> None:
        pop = PopulationFrequencies(global_af=0.00005, afr=0.00008, nfe=0.00003)
        sv = _make_variant(
            allele_frequency=0.00005,
            population_frequencies=pop,
        )
        tags = _classify(sv)
        assert EvidenceTag.PM2 in tags

    def test_does_not_fire_when_one_population_exceeds(self) -> None:
        pop = PopulationFrequencies(global_af=0.00005, afr=0.0005, nfe=0.00003)
        sv = _make_variant(
            allele_frequency=0.00005,
            population_frequencies=pop,
        )
        tags = _classify(sv)
        assert EvidenceTag.PM2 not in tags

    def test_falls_back_to_global_af_without_population_data(self) -> None:
        sv = _make_variant(allele_frequency=0.00005)
        tags = _classify(sv)
        assert EvidenceTag.PM2 in tags
