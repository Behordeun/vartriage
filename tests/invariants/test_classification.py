"""Hypothesis tests for ACMG classification (evidence tagging and combining rules).

Verifies that evidence tags are assigned based on criteria satisfaction,
and that combining rules produce the correct final classification.
"""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from tests.generators.variants import evidence_tag_set
from vartriage.classification.acmg import ACMGClassifier
from vartriage.classification.combining import combine_evidence
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

_PVS1_CONSEQUENCES = [FunctionalConsequence.NONSENSE, FunctionalConsequence.FRAMESHIFT]
_NON_PVS1_CONSEQUENCES = [
    c for c in FunctionalConsequence if c not in _PVS1_CONSEQUENCES
]


@st.composite
def scored_variant_for_classification(draw: st.DrawFn) -> ScoredVariant:
    """Generate a ScoredVariant with controlled fields for classification testing.

    Produces variants with various combinations of:
    - Consequence (PVS1-triggering vs. non-triggering)
    - Allele frequency (below/above PM2 threshold, or None)
    - REVEL score (above/below PP3 threshold, or None)
    - ClinVar assertion (Pathogenic, Benign, VUS, None)
    """
    variant = Variant(
        chrom=draw(st.sampled_from([f"chr{i}" for i in range(1, 23)])),
        pos=draw(st.integers(min_value=1, max_value=250_000_000)),
        id=None,
        ref=draw(st.sampled_from(["A", "C", "G", "T"])),
        alt=draw(st.sampled_from(["A", "C", "G", "T"])),
        qual=30.0,
        filter_status="PASS",
        info={},
    )

    consequence = draw(st.sampled_from(list(FunctionalConsequence)))

    # Allele frequency: None (missing) or a value in [0, 1]
    af_available = draw(st.booleans())
    if af_available:
        allele_frequency = draw(
            st.floats(min_value=0.0, max_value=1.0, allow_nan=False)
        )
        frequency_unknown = False
    else:
        allele_frequency = None
        frequency_unknown = True

    # ClinVar: None (missing) or one of the assertion values
    clinvar_available = draw(st.booleans())
    if clinvar_available:
        clinvar_assertion = draw(st.sampled_from(list(ClinVarAssertion)))
        clinvar_unknown = False
    else:
        clinvar_assertion = None
        clinvar_unknown = True

    annotated = AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
        frequency_unknown=frequency_unknown,
        clinvar_unknown=clinvar_unknown,
    )

    # REVEL score: None (missing) or a value in [0, 1]
    revel_available = draw(st.booleans())
    if revel_available:
        revel_score = draw(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))
    else:
        revel_score = None

    cadd_phred = draw(
        st.one_of(
            st.none(),
            st.floats(min_value=0.0, max_value=60.0, allow_nan=False),
        )
    )
    cadd_normalized = min(cadd_phred / 99.0, 1.0) if cadd_phred is not None else None

    # SpliceAI score: None (missing) or a value in [0, 1]
    spliceai_available = draw(st.booleans())
    if spliceai_available:
        spliceai_score = draw(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))
    else:
        spliceai_score = None

    composite_rank = None
    if revel_score is not None and cadd_normalized is not None:
        composite_rank = (revel_score * 0.6) + (cadd_normalized * 0.4)
    elif revel_score is not None:
        composite_rank = revel_score
    elif cadd_normalized is not None:
        composite_rank = cadd_normalized

    return ScoredVariant(
        annotated=annotated,
        cadd_phred=cadd_phred,
        cadd_normalized=cadd_normalized,
        revel_score=revel_score,
        spliceai_score=spliceai_score,
        composite_rank=composite_rank,
    )


# ---------------------------------------------------------------------------


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_pvs1_assigned_iff_nonsense_or_frameshift(variant: ScoredVariant) -> None:
    """PVS1 is assigned for Nonsense/Frameshift or SPLICE_SITE with SpliceAI > 0.8."""
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    consequence = variant.annotated.consequence
    spliceai = variant.spliceai_score

    pvs1_expected = consequence in (
        FunctionalConsequence.NONSENSE,
        FunctionalConsequence.FRAMESHIFT,
    ) or (
        consequence == FunctionalConsequence.SPLICE_SITE
        and spliceai is not None
        and spliceai > 0.8
    )

    if pvs1_expected:
        assert EvidenceTag.PVS1 in classified.evidence_tags, (
            f"PVS1 should be assigned for {consequence.value} (spliceai={spliceai})"
        )
    else:
        assert EvidenceTag.PVS1 not in classified.evidence_tags, (
            f"PVS1 should NOT be assigned for {consequence.value} (spliceai={spliceai})"
        )


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_pm2_assigned_iff_af_below_threshold(variant: ScoredVariant) -> None:
    """PM2 is assigned when AF < 0.0001 or when AF is None (absent from controls)."""
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    af = variant.annotated.allele_frequency

    if af is None:
        # Absent from gnomAD = not observed in controls. PM2 fires.
        assert EvidenceTag.PM2 in classified.evidence_tags, (
            "PM2 should fire when allele frequency is None (absent from controls)"
        )
    elif af < 0.0001:
        assert EvidenceTag.PM2 in classified.evidence_tags, (
            f"PM2 should be assigned for AF={af} < 0.0001"
        )
    else:
        assert EvidenceTag.PM2 not in classified.evidence_tags, (
            f"PM2 should NOT be assigned for AF={af} >= 0.0001"
        )


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_pp3_assigned_iff_revel_or_spliceai_triggers(
    variant: ScoredVariant,
) -> None:
    """PP3 assigned when REVEL > 0.7, or SpliceAI > 0.5 on splice-adjacent."""
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    revel = variant.revel_score
    spliceai = variant.spliceai_score
    consequence = variant.annotated.consequence

    revel_available = revel is not None
    spliceai_available = spliceai is not None

    splice_adjacent = consequence in (
        FunctionalConsequence.SPLICE_SITE,
        FunctionalConsequence.MISSENSE,
    )

    # PP3 cannot fire if neither predictor is available
    if not revel_available and not spliceai_available:
        assert EvidenceTag.PP3 not in classified.evidence_tags
        assert EvidenceTag.PP3_MODERATE not in classified.evidence_tags
        assert EvidenceTag.PP3_STRONG not in classified.evidence_tags
        return

    # REVEL path: strong (> 0.932), moderate (> 0.773), or supporting (> 0.644)
    if revel_available and revel > 0.932:
        assert EvidenceTag.PP3_STRONG in classified.evidence_tags, (
            f"PP3_Strong should be assigned for REVEL={revel} > 0.932"
        )
        return

    if revel_available and revel > 0.773:
        assert EvidenceTag.PP3_MODERATE in classified.evidence_tags, (
            f"PP3_Moderate should be assigned for REVEL={revel} > 0.773"
        )
        return

    if revel_available and revel > 0.644:
        assert EvidenceTag.PP3 in classified.evidence_tags, (
            f"PP3 should be assigned for REVEL={revel} > 0.644"
        )
        return

    # SpliceAI path triggers PP3 on splice-adjacent
    if spliceai_available and spliceai > 0.5 and splice_adjacent:
        assert EvidenceTag.PP3 in classified.evidence_tags, (
            f"PP3 should be assigned for SpliceAI={spliceai} > 0.5 "
            f"on {consequence.value}"
        )
        return

    # Neither trigger met
    assert EvidenceTag.PP3 not in classified.evidence_tags, (
        f"PP3 should NOT be assigned: REVEL={revel}, "
        f"SpliceAI={spliceai}, consequence={consequence.value}"
    )
    assert EvidenceTag.PP3_MODERATE not in classified.evidence_tags, (
        f"PP3_Moderate should NOT be assigned: REVEL={revel}"
    )
    assert EvidenceTag.PP3_STRONG not in classified.evidence_tags, (
        f"PP3_Strong should NOT be assigned: REVEL={revel}"
    )


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_pp5_assigned_iff_clinvar_pathogenic(variant: ScoredVariant) -> None:
    """PP5 assigned when ClinVar is Pathogenic; omitted when None (missing)."""
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    assertion = variant.annotated.clinvar_assertion

    if assertion is None:
        assert EvidenceTag.PP5 not in classified.evidence_tags, (
            "PP5 should be omitted when ClinVar data is unavailable"
        )
        assert "ClinVar" in classified.missing_data_sources, (
            "ClinVar should be listed as missing when assertion is None"
        )
    elif assertion == ClinVarAssertion.PATHOGENIC:
        assert EvidenceTag.PP5 in classified.evidence_tags, (
            "PP5 should be assigned when ClinVar is Pathogenic"
        )
    else:
        assert EvidenceTag.PP5 not in classified.evidence_tags, (
            f"PP5 should NOT be assigned for ClinVar={assertion.value}"
        )


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_tag_set_is_exactly_satisfied_criteria(variant: ScoredVariant) -> None:
    """The full evidence tag set contains exactly the tags whose criteria are met.

    This is the completeness check: no extra tags, no missing tags.
    """
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    expected_tags: set[EvidenceTag] = set()

    consequence = variant.annotated.consequence
    spliceai = variant.spliceai_score

    # PVS1: consequence in {Nonsense, Frameshift} OR SPLICE_SITE + SpliceAI > 0.8
    if consequence in (
        FunctionalConsequence.NONSENSE,
        FunctionalConsequence.FRAMESHIFT,
    ) or (
        consequence == FunctionalConsequence.SPLICE_SITE
        and spliceai is not None
        and spliceai > 0.8
    ):
        expected_tags.add(EvidenceTag.PVS1)

    # PM4: in-frame indel or stop-loss
    if consequence in (
        FunctionalConsequence.IN_FRAME_INSERTION,
        FunctionalConsequence.IN_FRAME_DELETION,
        FunctionalConsequence.STOP_LOSS,
    ):
        expected_tags.add(EvidenceTag.PM4)

    # PM2: AF < 0.0001 or absent from gnomAD (AF=None)
    af = variant.annotated.allele_frequency
    pop_freq = variant.annotated.population_frequencies

    if pop_freq is not None:
        has_any_data = any(
            v is not None
            for v in (
                pop_freq.afr,
                pop_freq.amr,
                pop_freq.asj,
                pop_freq.eas,
                pop_freq.fin,
                pop_freq.nfe,
                pop_freq.sas,
                pop_freq.global_af,
            )
        )
        if not has_any_data:
            # All population fields None = absent from gnomAD → PM2
            expected_tags.add(EvidenceTag.PM2)
        elif pop_freq.all_below(0.0001):
            expected_tags.add(EvidenceTag.PM2)
    elif af is not None and af < 0.0001:
        expected_tags.add(EvidenceTag.PM2)
    elif af is None:
        # Absent from gnomAD = not observed in controls → PM2
        expected_tags.add(EvidenceTag.PM2)

    # PP3: ClinGen-calibrated REVEL thresholds (Pejaver et al. 2022)
    # Strong: REVEL > 0.932
    # Moderate: REVEL > 0.773
    # Supporting: REVEL > 0.644
    # OR SpliceAI > 0.5 on splice-adjacent (supporting only)
    revel = variant.revel_score
    revel_available = revel is not None
    spliceai_available = spliceai is not None

    if revel_available or spliceai_available:
        if revel_available and revel > 0.932:
            expected_tags.add(EvidenceTag.PP3_STRONG)
        elif revel_available and revel > 0.773:
            expected_tags.add(EvidenceTag.PP3_MODERATE)
        elif (
            revel_available
            and revel > 0.644
            or (
                spliceai_available
                and spliceai > 0.5
                and consequence
                in (
                    FunctionalConsequence.SPLICE_SITE,
                    FunctionalConsequence.MISSENSE,
                )
            )
        ):
            expected_tags.add(EvidenceTag.PP3)

    # PP5: ClinVar == Pathogenic (skip if assertion is None)
    assertion = variant.annotated.clinvar_assertion
    if assertion == ClinVarAssertion.PATHOGENIC:
        expected_tags.add(EvidenceTag.PP5)

    # BA1: any population AF > 5% (or global AF > 5% when no population data)
    pop_freq = variant.annotated.population_frequencies
    if pop_freq is not None:
        if pop_freq.any_exceeds(0.05):
            expected_tags.add(EvidenceTag.BA1)
    elif af is not None and af > 0.05:
        expected_tags.add(EvidenceTag.BA1)

    # BS1: any population AF > 1% (only if BA1 not already assigned)
    if EvidenceTag.BA1 not in expected_tags:
        if pop_freq is not None:
            if pop_freq.any_exceeds(0.01):
                expected_tags.add(EvidenceTag.BS1)
        elif af is not None and af > 0.01:
            expected_tags.add(EvidenceTag.BS1)

    # BP4: computational benign (ClinGen-calibrated thresholds)
    # Missense: REVEL < 0.183 (moderate) or REVEL < 0.290 (supporting)
    # Non-missense: CADD < 10 (supporting)
    # Does NOT fire for null/protein-altering variants
    if consequence not in (
        FunctionalConsequence.FRAMESHIFT,
        FunctionalConsequence.NONSENSE,
        FunctionalConsequence.STOP_LOSS,
        FunctionalConsequence.IN_FRAME_INSERTION,
        FunctionalConsequence.IN_FRAME_DELETION,
    ):
        if consequence == FunctionalConsequence.MISSENSE:
            if revel is not None:
                if revel < 0.183:
                    expected_tags.add(EvidenceTag.BP4_MODERATE)
                elif revel < 0.290:
                    expected_tags.add(EvidenceTag.BP4)
        else:
            cadd = variant.cadd_phred
            if cadd is not None and cadd < 10.0:
                expected_tags.add(EvidenceTag.BP4)

    # BP7: synonymous + SpliceAI < 0.1
    if (
        consequence == FunctionalConsequence.SYNONYMOUS
        and spliceai is not None
        and spliceai < 0.1
    ):
        expected_tags.add(EvidenceTag.BP7)

    assert classified.evidence_tags == frozenset(expected_tags), (
        f"Expected tags {expected_tags}, got {classified.evidence_tags}"
    )


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_missing_sources_reported_correctly(variant: ScoredVariant) -> None:
    """Missing data sources are reported when required data is unavailable."""
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    expected_missing: set[str] = set()

    # gnomAD: AF=None now fires PM2 (absent = rare), so it's no longer "missing"
    # gnomAD is only missing when pop_freq has partial data issues (not testable here)

    # ClinVar missing when assertion is None
    if variant.annotated.clinvar_assertion is None:
        expected_missing.add("ClinVar")

    # PP3 missing source tracking
    revel = variant.revel_score
    spliceai = variant.spliceai_score
    consequence = variant.annotated.consequence
    revel_available = revel is not None
    spliceai_available = spliceai is not None

    if not revel_available and not spliceai_available:
        # Both unavailable, both recorded
        expected_missing.add("REVEL")
        expected_missing.add("SpliceAI")
    elif revel_available and revel > 0.644:
        # REVEL triggered PP3 or PP3_Moderate, no missing sources from PP3
        pass
    elif (
        spliceai_available
        and spliceai > 0.5
        and consequence
        in (
            FunctionalConsequence.SPLICE_SITE,
            FunctionalConsequence.MISSENSE,
        )
    ):
        # SpliceAI triggered PP3, no missing sources from PP3
        pass
    else:
        # Neither triggered, record whichever is unavailable
        if not revel_available:
            expected_missing.add("REVEL")
        if not spliceai_available:
            expected_missing.add("SpliceAI")

    # PVS1 missing source tracking for SPLICE_SITE
    if (
        consequence == FunctionalConsequence.SPLICE_SITE
        and consequence
        not in (
            FunctionalConsequence.NONSENSE,
            FunctionalConsequence.FRAMESHIFT,
        )
        and spliceai is None
    ):
        expected_missing.add("SpliceAI")

    # PS1/PM5 missing source tracking for MISSENSE variants
    # Without protein_change (no codon resolution), classifier reports it as missing
    if consequence == FunctionalConsequence.MISSENSE:
        protein_change = variant.annotated.protein_change
        if protein_change is None:
            expected_missing.add("codon_resolution")
        else:
            # protein_change present but no protein_index → index is missing
            expected_missing.add("ClinVar_protein_index")

    # PM1 missing source tracking: missense without gene constraint data
    if consequence == FunctionalConsequence.MISSENSE:
        gene_context = variant.annotated.gene_context
        if gene_context is None or gene_context.constraint is None:
            expected_missing.add("gnomAD_constraint")

    assert classified.missing_data_sources == frozenset(expected_missing), (
        f"Expected missing sources {expected_missing}, "
        f"got {classified.missing_data_sources}"
    )


# ---------------------------------------------------------------------------


@given(data=st.data())
@settings(max_examples=100)
def test_empty_tags_produce_vus(data: st.DataObject) -> None:
    """Empty tag sets always produce VUS classification."""
    tags: frozenset[EvidenceTag] = frozenset()
    result = combine_evidence(tags)
    assert result == ACMGClassification.VUS, (
        f"Empty tag set should yield VUS, got {result.value}"
    )


@given(tags=evidence_tag_set())
@settings(max_examples=200)
def test_combining_matches_point_engine(
    tags: frozenset[EvidenceTag],
) -> None:
    """combine_evidence agrees with the SVI point engine on every tag set.

    combine_evidence is a thin delegation to classify_by_points, so the
    two must never disagree. This pins the delegation and guards against a
    future divergence.
    """
    from vartriage.classification.points import classify_by_points

    assert combine_evidence(tags) == classify_by_points(tags)


@given(tags=evidence_tag_set())
@settings(max_examples=200)
def test_ba1_forces_benign(tags: frozenset[EvidenceTag]) -> None:
    """BA1 is a stand-alone Benign override regardless of other evidence."""
    with_ba1 = tags | {EvidenceTag.BA1}
    assert combine_evidence(with_ba1) == ACMGClassification.BENIGN


@given(tags=evidence_tag_set())
@settings(max_examples=200)
def test_adding_pathogenic_support_never_lowers_the_tier(
    tags: frozenset[EvidenceTag],
) -> None:
    """Adding a supporting pathogenic criterion cannot move the tier down.

    Monotonicity of the point sum: PP3 adds a positive point, so the
    resulting classification is at least as pathogenic as before (unless a
    BA1 override is present, which pins Benign either way).
    """
    order = [
        ACMGClassification.BENIGN,
        ACMGClassification.LIKELY_BENIGN,
        ACMGClassification.VUS,
        ACMGClassification.LIKELY_PATHOGENIC,
        ACMGClassification.PATHOGENIC,
    ]
    if EvidenceTag.BA1 in tags:
        return
    before = combine_evidence(tags)
    after = combine_evidence(tags | {EvidenceTag.PP3})
    assert order.index(after) >= order.index(before)


@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_classifier_output_matches_combining_rules(
    variant: ScoredVariant,
) -> None:
    """The ACMGClassifier final classification matches combining rules applied
    to the assigned evidence tags.

    This is the end-to-end property: tag assignment feeds into combining,
    and the output classification must be consistent.
    """
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    # Re-derive classification from the tags assigned
    expected_classification = combine_evidence(classified.evidence_tags)

    assert classified.classification == expected_classification, (
        f"Classification {classified.classification.value} does not match "
        f"combine_evidence result {expected_classification.value} "
        f"for tags {[t.value for t in classified.evidence_tags]}"
    )


@given(data=st.data())
@settings(max_examples=100)
def test_pvs1_plus_pp3_pp5_yields_pathogenic(data: st.DataObject) -> None:
    """PVS1 (Very Strong) + PP3 + PP5 (2 Supporting) yields Pathogenic.

    This tests a specific combining path: >=1 VS + >=2 Supporting.
    """
    tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PP3, EvidenceTag.PP5})
    result = combine_evidence(tags)
    assert result == ACMGClassification.PATHOGENIC, (
        f"PVS1+PP3+PP5 should be Pathogenic, got {result.value}"
    )


@given(data=st.data())
@settings(max_examples=100)
def test_pvs1_plus_pm2_yields_pathogenic(data: st.DataObject) -> None:
    """PVS1 (8) + PM2 (2) sums to 10 points, which is Pathogenic."""
    tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PM2})
    result = combine_evidence(tags)
    assert result == ACMGClassification.PATHOGENIC, (
        f"PVS1+PM2 should be Pathogenic, got {result.value}"
    )


@given(data=st.data())
@settings(max_examples=100)
def test_single_supporting_tag_yields_vus(data: st.DataObject) -> None:
    """A single Supporting tag (PP3 or PP5) alone yields VUS.

    No combining rule is met with just one supporting tag.
    """
    tag = data.draw(st.sampled_from([EvidenceTag.PP3, EvidenceTag.PP5]))
    tags = frozenset({tag})
    result = combine_evidence(tags)
    assert result == ACMGClassification.VUS, (
        f"Single {tag.value} should yield VUS, got {result.value}"
    )
