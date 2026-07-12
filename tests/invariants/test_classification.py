"""Hypothesis tests for ACMG classification (evidence tagging and combining rules).

Verifies that evidence tags are assigned based on criteria satisfaction,
and that combining rules produce the correct final classification.
"""

from __future__ import annotations

from hypothesis import given, settings, assume
from hypothesis import strategies as st

from vartriage.classification.acmg import ACMGClassifier
from vartriage.classification.combining import combine_evidence
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClinVarAssertion,
    EvidenceStrength,
    EvidenceTag,
    EVIDENCE_STRENGTH_MAP,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)

from tests.generators.variants import scored_variant, evidence_tag_set

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
        revel_score = draw(
            st.floats(min_value=0.0, max_value=1.0, allow_nan=False)
        )
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
        spliceai_score = draw(
            st.floats(min_value=0.0, max_value=1.0, allow_nan=False)
        )
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
    """PVS1 is assigned for Nonsense/Frameshift or SPLICE_SITE with SpliceAI > 0.8.
    """
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    consequence = variant.annotated.consequence
    spliceai = variant.spliceai_score

    pvs1_expected = (
        consequence in (
            FunctionalConsequence.NONSENSE,
            FunctionalConsequence.FRAMESHIFT,
        )
        or (
            consequence == FunctionalConsequence.SPLICE_SITE
            and spliceai is not None
            and spliceai > 0.8
        )
    )

    if pvs1_expected:
        assert EvidenceTag.PVS1 in classified.evidence_tags, (
            f"PVS1 should be assigned for {consequence.value} "
            f"(spliceai={spliceai})"
        )
    else:
        assert EvidenceTag.PVS1 not in classified.evidence_tags, (
            f"PVS1 should NOT be assigned for {consequence.value} "
            f"(spliceai={spliceai})"
        )

@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_pm2_assigned_iff_af_below_threshold(variant: ScoredVariant) -> None:
    """PM2 is assigned when AF < 0.0001; omitted when AF is None (missing).
    """
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    af = variant.annotated.allele_frequency

    if af is None:
        # Data unavailable, PM2 should be omitted
        assert EvidenceTag.PM2 not in classified.evidence_tags, (
            "PM2 should be omitted when allele frequency is unavailable"
        )
        assert "gnomAD" in classified.missing_data_sources, (
            "gnomAD should be listed as missing when AF is None"
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
    """PP3 assigned when REVEL > 0.7, or SpliceAI > 0.5 on splice-adjacent.
    """
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
        return

    # REVEL path triggers PP3
    if revel_available and revel > 0.7:
        assert EvidenceTag.PP3 in classified.evidence_tags, (
            f"PP3 should be assigned for REVEL={revel} > 0.7"
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

@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_pp5_assigned_iff_clinvar_pathogenic(variant: ScoredVariant) -> None:
    """PP5 assigned when ClinVar is Pathogenic; omitted when None (missing).
    """
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
    ):
        expected_tags.add(EvidenceTag.PVS1)
    elif (
        consequence == FunctionalConsequence.SPLICE_SITE
        and spliceai is not None
        and spliceai > 0.8
    ):
        expected_tags.add(EvidenceTag.PVS1)

    # PM2: AF < 0.0001 (skip if AF is None)
    af = variant.annotated.allele_frequency
    if af is not None and af < 0.0001:
        expected_tags.add(EvidenceTag.PM2)

    # PP3: REVEL > 0.7 OR SpliceAI > 0.5 on splice-adjacent
    revel = variant.revel_score
    revel_available = revel is not None
    spliceai_available = spliceai is not None

    if revel_available or spliceai_available:
        if revel_available and revel > 0.7:
            expected_tags.add(EvidenceTag.PP3)
        elif (
            spliceai_available
            and spliceai > 0.5
            and consequence in (
                FunctionalConsequence.SPLICE_SITE,
                FunctionalConsequence.MISSENSE,
            )
        ):
            expected_tags.add(EvidenceTag.PP3)

    # PP5: ClinVar == Pathogenic (skip if assertion is None)
    assertion = variant.annotated.clinvar_assertion
    if assertion == ClinVarAssertion.PATHOGENIC:
        expected_tags.add(EvidenceTag.PP5)

    assert classified.evidence_tags == frozenset(expected_tags), (
        f"Expected tags {expected_tags}, got {classified.evidence_tags}"
    )

@given(variant=scored_variant_for_classification())
@settings(max_examples=200)
def test_missing_sources_reported_correctly(variant: ScoredVariant) -> None:
    """Missing data sources are reported when required data is unavailable.
    """
    classifier = ACMGClassifier()
    results = list(classifier.classify(iter([variant])))
    classified = results[0]

    expected_missing: set[str] = set()

    # gnomAD missing when AF is None
    if variant.annotated.allele_frequency is None:
        expected_missing.add("gnomAD")

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
    elif revel_available and revel > 0.7:
        # REVEL triggered PP3, no missing sources from PP3
        pass
    elif (
        spliceai_available
        and spliceai > 0.5
        and consequence in (
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
    if consequence == FunctionalConsequence.SPLICE_SITE:
        if consequence not in (
            FunctionalConsequence.NONSENSE,
            FunctionalConsequence.FRAMESHIFT,
        ):
            if spliceai is None:
                expected_missing.add("SpliceAI")

    assert classified.missing_data_sources == frozenset(expected_missing), (
        f"Expected missing sources {expected_missing}, "
        f"got {classified.missing_data_sources}"
    )

# ---------------------------------------------------------------------------

@given(data=st.data())
@settings(max_examples=100)
def test_empty_tags_produce_vus(data: st.DataObject) -> None:
    """Empty tag sets always produce VUS classification.
    """
    tags: frozenset[EvidenceTag] = frozenset()
    result = combine_evidence(tags)
    assert result == ACMGClassification.VUS, (
        f"Empty tag set should yield VUS, got {result.value}"
    )

@given(tags=evidence_tag_set())
@settings(max_examples=200)
def test_combining_rules_match_specification(
    tags: frozenset[EvidenceTag],
) -> None:
    """Combining rules produce the correct classification per ACMG/AMP 2015.

    Verifies:
    - Pathogenic: >=1 VS + >=1 S, OR >=2 S + >=1 Sup, OR >=1 VS + >=2 Sup
    - Likely Pathogenic: 1 VS + 1 M, OR 1 S + 1-2 M, OR 1 S + >=2 Sup
    - VUS: default
    """
    result = combine_evidence(tags)

    # Compute expected classification from the tag strengths
    counts: dict[EvidenceStrength, int] = {
        EvidenceStrength.VERY_STRONG: 0,
        EvidenceStrength.STRONG: 0,
        EvidenceStrength.MODERATE: 0,
        EvidenceStrength.SUPPORTING: 0,
    }
    for tag in tags:
        strength = EVIDENCE_STRENGTH_MAP[tag]
        counts[strength] += 1

    vs = counts[EvidenceStrength.VERY_STRONG]
    s = counts[EvidenceStrength.STRONG]
    m = counts[EvidenceStrength.MODERATE]
    sup = counts[EvidenceStrength.SUPPORTING]

    # Pathogenic rules
    is_pathogenic = (
        (vs >= 1 and s >= 1)
        or (s >= 2 and sup >= 1)
        or (vs >= 1 and sup >= 2)
    )

    # Likely Pathogenic rules
    is_likely_pathogenic = (
        (vs >= 1 and m >= 1)
        or (s >= 1 and 1 <= m <= 2)
        or (s >= 1 and sup >= 2)
    )

    if not tags:
        expected = ACMGClassification.VUS
    elif is_pathogenic:
        expected = ACMGClassification.PATHOGENIC
    elif is_likely_pathogenic:
        expected = ACMGClassification.LIKELY_PATHOGENIC
    else:
        expected = ACMGClassification.VUS

    assert result == expected, (
        f"For tags {[t.value for t in tags]} "
        f"(VS={vs}, S={s}, M={m}, Sup={sup}): "
        f"expected {expected.value}, got {result.value}"
    )

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
def test_pvs1_plus_pm2_yields_likely_pathogenic(data: st.DataObject) -> None:
    """PVS1 (Very Strong) + PM2 (Moderate) yields Likely Pathogenic.

    Tests the combining path: 1 VS + 1 Moderate.
    """
    tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PM2})
    result = combine_evidence(tags)
    assert result == ACMGClassification.LIKELY_PATHOGENIC, (
        f"PVS1+PM2 should be Likely_Pathogenic, got {result.value}"
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
