"""Unit tests for ACMG evidence tag assignment."""

from __future__ import annotations

from vartriage.annotation.clinvar_protein_index import (
    ClinVarProteinIndex,
    PathogenicMissense,
)
from vartriage.classification.acmg import ACMGClassifier
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ProteinChange,
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
    protein_change: ProteinChange | None = None,
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
        protein_change=protein_change,
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
        # REVEL 0.7 is above supporting threshold (0.644) but below moderate (0.773)
        sv = _make_scored_variant(revel_score=0.7)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP3 in results[0].evidence_tags

    def test_assigns_pp3_moderate_for_very_high_revel(self) -> None:
        # REVEL 0.85 is above moderate threshold (0.773)
        sv = _make_scored_variant(revel_score=0.85)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert EvidenceTag.PP3_MODERATE in results[0].evidence_tags
        assert EvidenceTag.PP3 not in results[0].evidence_tags

    def test_does_not_assign_pp3_at_threshold(self) -> None:
        # REVEL 0.644 is at the boundary, not above — should not fire
        sv = _make_scored_variant(revel_score=0.644)
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
        sv = _make_scored_variant(clinvar_assertion=ClinVarAssertion.LIKELY_PATHOGENIC)
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
        # Uses NONSENSE consequence to avoid PS1/PM5 (which need protein_change data).
        # For missense with full coverage, see test_no_missing_sources_missense_with_protein_change.
        sv = _make_scored_variant(
            consequence=FunctionalConsequence.NONSENSE,
            allele_frequency=0.0005,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC,
            revel_score=0.85,
        )
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert len(results[0].missing_data_sources) == 0

    def test_no_missing_sources_missense_with_protein_change(self) -> None:
        """Missense with protein_change reports no missing data for PS1/PM5."""
        sv = _make_scored_variant(
            allele_frequency=0.0005,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC,
            revel_score=0.85,
            protein_change=ProteinChange(
                gene_name="BRCA1", position=100, reference_aa="R", altered_aa="H"
            ),
        )
        # Providing a protein_index=None still triggers missing for the index,
        # but codon_resolution is satisfied since protein_change is populated.
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert results[0].missing_data_sources == frozenset(
            {"ClinVar_protein_index", "functional_domain"}
        )


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
        # REVEL 0.9 > 0.773 fires PP3_Moderate (ClinGen-calibrated moderate threshold)
        assert EvidenceTag.PP3_MODERATE in tags
        assert EvidenceTag.PP5 in tags

    def test_scored_variant_preserved_in_output(self) -> None:
        sv = _make_scored_variant(revel_score=0.9)
        classifier = ACMGClassifier()
        results = list(classifier.classify(iter([sv])))
        assert results[0].scored is sv


def _make_scored_variant_at(
    chrom: str = "chr17",
    pos: int = 43091429,
    ref: str = "T",
    alt: str = "G",
    protein_change: ProteinChange | None = None,
    allele_frequency: float | None = 0.005,
    revel_score: float | None = 0.5,
) -> ScoredVariant:
    """Helper that lets us set specific genomic coordinates for PS1/PM5 tests."""
    v = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=v,
        consequence=FunctionalConsequence.MISSENSE,
        allele_frequency=allele_frequency,
        clinvar_assertion=None,
        frequency_unknown=False,
        clinvar_unknown=False,
        protein_change=protein_change,
    )
    cadd_phred = 25.0
    cadd_normalized = min(cadd_phred / 99.0, 1.0)
    return ScoredVariant(
        annotated=annotated,
        cadd_phred=cadd_phred,
        cadd_normalized=cadd_normalized,
        revel_score=revel_score,
        composite_rank=None,
    )


def _build_protein_index(entries: list[PathogenicMissense]) -> ClinVarProteinIndex:
    """Build an in-memory ClinVarProteinIndex from a list of PathogenicMissense."""
    return ClinVarProteinIndex.from_variants(entries)


class TestPS1Assignment:
    """PS1: same amino acid change as known pathogenic, different nucleotide."""

    def test_ps1_fires_for_same_aa_change_different_nucleotide(self) -> None:
        """Variant at different codon producing same AA change gets PS1."""
        # ClinVar pathogenic: chr17:43091429 T>C produces BRCA1 p.M1775R
        known = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="R",
            chrom="chr17",
            genomic_pos=43091429,
            ref_allele="T",
            alt_allele="C",
        )
        protein_index = _build_protein_index([known])

        # Query variant: chr17:43091430 A>G also produces p.M1775R (different codon)
        sv = _make_scored_variant_at(
            chrom="chr17",
            pos=43091430,
            ref="A",
            alt="G",
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="R"
            ),
        )

        classifier = ACMGClassifier(protein_index=protein_index)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        assert EvidenceTag.PS1 in tags
        assert EvidenceTag.PM5 not in tags

    def test_ps1_does_not_fire_for_same_nucleotide_change(self) -> None:
        """Exact same variant as the known pathogenic should not get PS1."""
        known = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="R",
            chrom="chr17",
            genomic_pos=43091429,
            ref_allele="T",
            alt_allele="C",
        )
        protein_index = _build_protein_index([known])

        # Query is identical to the known pathogenic entry
        sv = _make_scored_variant_at(
            chrom="chr17",
            pos=43091429,
            ref="T",
            alt="C",
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="R"
            ),
        )

        classifier = ACMGClassifier(protein_index=protein_index)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        assert EvidenceTag.PS1 not in tags

    def test_ps1_does_not_fire_for_different_aa_change(self) -> None:
        """Different amino acid substitution at same position should not get PS1."""
        known = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="R",
            chrom="chr17",
            genomic_pos=43091429,
            ref_allele="T",
            alt_allele="C",
        )
        protein_index = _build_protein_index([known])

        # Different alt AA: M1775K instead of M1775R
        sv = _make_scored_variant_at(
            chrom="chr17",
            pos=43091431,
            ref="G",
            alt="A",
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="K"
            ),
        )

        classifier = ACMGClassifier(protein_index=protein_index)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        assert EvidenceTag.PS1 not in tags

    def test_ps1_does_not_fire_without_protein_index(self) -> None:
        """Without a loaded protein index, PS1 cannot be evaluated."""
        sv = _make_scored_variant_at(
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="R"
            ),
        )

        classifier = ACMGClassifier(protein_index=None)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        assert EvidenceTag.PS1 not in tags
        assert "ClinVar_protein_index" in results[0].missing_data_sources


class TestPM5Assignment:
    """PM5: different amino acid change at same position as known pathogenic missense."""

    def test_pm5_fires_for_different_aa_at_same_position(self) -> None:
        """Novel missense at a position with a known pathogenic change gets PM5."""
        # Known pathogenic: BRCA1 p.M1775R
        known = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="R",
            chrom="chr17",
            genomic_pos=43091429,
            ref_allele="T",
            alt_allele="C",
        )
        protein_index = _build_protein_index([known])

        # Query: different substitution at same position — p.M1775K
        sv = _make_scored_variant_at(
            chrom="chr17",
            pos=43091431,
            ref="G",
            alt="A",
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="K"
            ),
        )

        classifier = ACMGClassifier(protein_index=protein_index)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        assert EvidenceTag.PM5 in tags
        assert EvidenceTag.PS1 not in tags

    def test_pm5_suppressed_when_ps1_fires(self) -> None:
        """When PS1 is assigned, PM5 should not fire (PS1 is stronger)."""
        # Known pathogenic: M1775R via T>C and M1775K via G>A
        known_r = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="R",
            chrom="chr17",
            genomic_pos=43091429,
            ref_allele="T",
            alt_allele="C",
        )
        known_k = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="K",
            chrom="chr17",
            genomic_pos=43091431,
            ref_allele="G",
            alt_allele="A",
        )
        protein_index = _build_protein_index([known_r, known_k])

        # Query variant: same AA change as known_r (M1775R) but different nucleotide
        sv = _make_scored_variant_at(
            chrom="chr17",
            pos=43091430,
            ref="A",
            alt="G",
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="R"
            ),
        )

        classifier = ACMGClassifier(protein_index=protein_index)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        # PS1 fires because same AA change via different nucleotide
        assert EvidenceTag.PS1 in tags
        # PM5 suppressed by the PS1 guard
        assert EvidenceTag.PM5 not in tags

    def test_pm5_does_not_fire_for_same_aa_change(self) -> None:
        """Same amino acid change (even same nucleotide) should not get PM5."""
        known = PathogenicMissense(
            gene="BRCA1",
            position=1775,
            ref_aa="M",
            alt_aa="R",
            chrom="chr17",
            genomic_pos=43091429,
            ref_allele="T",
            alt_allele="C",
        )
        protein_index = _build_protein_index([known])

        # Same exact variant — same AA, same nucleotide
        sv = _make_scored_variant_at(
            chrom="chr17",
            pos=43091429,
            ref="T",
            alt="C",
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="R"
            ),
        )

        classifier = ACMGClassifier(protein_index=protein_index)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        # Neither PS1 (same nucleotide) nor PM5 (same AA change) should fire
        assert EvidenceTag.PM5 not in tags
        assert EvidenceTag.PS1 not in tags

    def test_pm5_does_not_fire_without_protein_index(self) -> None:
        """Without a loaded protein index, PM5 cannot be evaluated."""
        sv = _make_scored_variant_at(
            protein_change=ProteinChange(
                gene_name="BRCA1", position=1775, reference_aa="M", altered_aa="K"
            ),
        )

        classifier = ACMGClassifier(protein_index=None)
        results = list(classifier.classify(iter([sv])))
        tags = results[0].evidence_tags

        assert EvidenceTag.PM5 not in tags
        assert "ClinVar_protein_index" in results[0].missing_data_sources
