"""Unit tests for the v0.19.0 clinical refinements in the ACMG classifier.

Covers PVS1 NMD-escape downgrade, disease-specific frequency thresholds,
PP5 review-status modulation, and the BS2 evaluator.
"""

from __future__ import annotations

from vartriage.annotation.transcript_index import TranscriptCDSIndex
from vartriage.classification.acmg import ACMGClassifier
from vartriage.knowledge.models import (
    DiseaseAssociation,
    GeneConstraint,
    GeneContext,
)
from vartriage.models.variant import (
    AnnotatedVariant,
    ClinVarAssertion,
    ClinVarReviewStatus,
    EvidenceTag,
    FunctionalConsequence,
    PopulationFrequencies,
    ScoredVariant,
    Variant,
)

_LOF_INTOLERANT = GeneConstraint(pli=0.99, loeuf=0.1, mis_z=1.0)
_AD = (DiseaseAssociation("Disorder A", "100000", "AD"),)
_AR = (DiseaseAssociation("Disorder B", "200000", "AR"),)


def _scored(
    consequence: FunctionalConsequence,
    *,
    gene: str | None = "GENE",
    pos: int = 250,
    chrom: str = "chr1",
    pop: PopulationFrequencies | None = None,
    clinvar: ClinVarAssertion | None = None,
    review: ClinVarReviewStatus | None = None,
    associations: tuple[DiseaseAssociation, ...] = (),
    constraint: GeneConstraint | None = None,
) -> ScoredVariant:
    gene_context = None
    if associations or constraint is not None:
        gene_context = GeneContext(
            disease_associations=associations, constraint=constraint
        )
    annotated = AnnotatedVariant(
        variant=Variant(chrom, pos, None, "A", "T", 99.0, "PASS", {}),
        consequence=consequence,
        gene_name=gene,
        population_frequencies=pop,
        clinvar_assertion=clinvar,
        clinvar_review_status=review,
        gene_context=gene_context,
    )
    return ScoredVariant(annotated=annotated)


def _plus_strand_two_exon_index() -> TranscriptCDSIndex:
    index = TranscriptCDSIndex()
    index.add_cds_exon("t1", "GENE", "chr1", 0, 100, "+", 0)
    index.add_cds_exon("t1", "GENE", "chr1", 200, 300, "+", 0)
    index.finalize()
    return index


def _single_exon_index() -> TranscriptCDSIndex:
    index = TranscriptCDSIndex()
    index.add_cds_exon("t1", "GENE", "chr1", 0, 300, "+", 0)
    index.finalize()
    return index


class TestPVS1NMDEscape:
    def test_last_exon_null_variant_downgrades_to_strong(self) -> None:
        classifier = ACMGClassifier(nmd_lookup=_plus_strand_two_exon_index())
        variant = _scored(
            FunctionalConsequence.NONSENSE, constraint=_LOF_INTOLERANT, pos=250
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags

    def test_single_exon_gene_downgrades_to_strong(self) -> None:
        classifier = ACMGClassifier(nmd_lookup=_single_exon_index())
        variant = _scored(
            FunctionalConsequence.FRAMESHIFT, constraint=_LOF_INTOLERANT, pos=150
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PVS1_STRONG in tags

    def test_variant_before_last_exon_stays_very_strong(self) -> None:
        classifier = ACMGClassifier(nmd_lookup=_plus_strand_two_exon_index())
        variant = _scored(
            FunctionalConsequence.NONSENSE, constraint=_LOF_INTOLERANT, pos=50
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PVS1 in tags

    def test_no_lookup_keeps_very_strong_and_records_source(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.NONSENSE, constraint=_LOF_INTOLERANT, pos=250
        )
        tags, missing = classifier._assign_tags(variant)
        assert EvidenceTag.PVS1 in tags
        assert "transcript_structure" in missing

    def test_downgrade_never_upgrades_a_strong_result(self) -> None:
        # No constraint data -> PVS1 resolves to Strong; NMD must not raise it.
        classifier = ACMGClassifier(nmd_lookup=_plus_strand_two_exon_index())
        variant = _scored(FunctionalConsequence.NONSENSE, pos=250)
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags


class TestDiseaseThresholds:
    def test_dominant_gene_uses_stricter_ba1(self) -> None:
        classifier = ACMGClassifier(use_disease_thresholds=True)
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.002),
            associations=_AD,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BA1 in tags

    def test_recessive_gene_uses_standard_ba1(self) -> None:
        classifier = ACMGClassifier(use_disease_thresholds=True)
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.002),
            associations=_AR,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BA1 not in tags

    def test_no_context_uses_default_thresholds(self) -> None:
        classifier = ACMGClassifier(use_disease_thresholds=True)
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.002),
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BA1 not in tags

    def test_flag_off_ignores_disease_context(self) -> None:
        classifier = ACMGClassifier(use_disease_thresholds=False)
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.002),
            associations=_AD,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BA1 not in tags


class TestPP5ReviewStatus:
    def test_expert_panel_upgrades_to_strong(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            clinvar=ClinVarAssertion.PATHOGENIC,
            review=ClinVarReviewStatus.EXPERT_PANEL,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PP5_STRONG in tags
        assert EvidenceTag.PP5 not in tags

    def test_single_submitter_stays_supporting(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            clinvar=ClinVarAssertion.PATHOGENIC,
            review=ClinVarReviewStatus.SINGLE_SUBMITTER,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PP5 in tags
        assert EvidenceTag.PP5_STRONG not in tags

    def test_no_criteria_does_not_fire_pp5(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            clinvar=ClinVarAssertion.PATHOGENIC,
            review=ClinVarReviewStatus.NO_CRITERIA,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PP5 not in tags
        assert EvidenceTag.PP5_STRONG not in tags

    def test_basic_format_without_review_status_keeps_supporting(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            clinvar=ClinVarAssertion.PATHOGENIC,
            review=None,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.PP5 in tags


class TestBS2:
    def test_dominant_gene_with_homozygotes_fires_bs2(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.001, hom_count=3),
            associations=_AD,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BS2 in tags

    def test_recessive_gene_does_not_fire_bs2(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.001, hom_count=3),
            associations=_AR,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BS2 not in tags

    def test_dominant_gene_missing_hom_count_records_source(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.001),
            associations=_AD,
        )
        tags, missing = classifier._assign_tags(variant)
        assert EvidenceTag.BS2 not in tags
        assert "gnomAD_homozygotes" in missing

    def test_no_context_does_not_record_homozygote_source(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(FunctionalConsequence.MISSENSE)
        _, missing = classifier._assign_tags(variant)
        assert "gnomAD_homozygotes" not in missing

    def test_zero_homozygotes_does_not_fire_bs2(self) -> None:
        classifier = ACMGClassifier()
        variant = _scored(
            FunctionalConsequence.MISSENSE,
            pop=PopulationFrequencies(global_af=0.001, hom_count=0),
            associations=_AD,
        )
        tags, _ = classifier._assign_tags(variant)
        assert EvidenceTag.BS2 not in tags
