"""ACMG/AMP evidence tag assignment for variant classification.

This module implements the ACMGClassifier, which evaluates scored variants
against ACMG/AMP 2015 evidence criteria and assigns appropriate evidence tags.
The classifier handles missing data gracefully by omitting tags when required
data sources are unavailable and recording which sources were missing.
"""

from __future__ import annotations

from typing import Iterator

from vartriage.classification.combining import combine_evidence
from vartriage.models.variant import (ClassifiedVariant, ClinVarAssertion,
                                      EvidenceTag, FunctionalConsequence,
                                      ScoredVariant)

_PVS1_CONSEQUENCES: frozenset[FunctionalConsequence] = frozenset(
    {
        FunctionalConsequence.NONSENSE,
        FunctionalConsequence.FRAMESHIFT,
    }
)

_PM2_AF_THRESHOLD: float = 0.0001

_PP3_REVEL_THRESHOLD: float = 0.7

_PP3_SPLICEAI_THRESHOLD: float = 0.5

_PVS1_SPLICEAI_THRESHOLD: float = 0.8

_BA1_AF_THRESHOLD: float = 0.05

_BS1_AF_THRESHOLD: float = 0.01

_BP4_REVEL_THRESHOLD: float = 0.15

_BP4_CADD_THRESHOLD: float = 10.0

_BP7_SPLICEAI_THRESHOLD: float = 0.1

_PP3_SPLICE_ADJACENT: frozenset[FunctionalConsequence] = frozenset(
    {
        FunctionalConsequence.SPLICE_SITE,
        FunctionalConsequence.MISSENSE,
    }
)

_PP5_CONFLICTING_ASSERTIONS: frozenset[ClinVarAssertion] = frozenset(
    {
        ClinVarAssertion.BENIGN,
        ClinVarAssertion.LIKELY_BENIGN,
    }
)


class ACMGClassifier:
    """Assign ACMG/AMP evidence tags and final classification.

    The classifier evaluates each ScoredVariant against four evidence criteria:

    - PVS1: Nonsense or Frameshift consequence (null variant)
    - PM2: gnomAD allele frequency below 0.0001 (absent from controls)
    - PP3: REVEL score above 0.7 (computational evidence)
    - PP5: ClinVar Pathogenic with no conflicting Benign/Likely_Benign

    When a required data source is unavailable for a given criterion, that
    tag is omitted and the source name is recorded in the output.

    Notes
    -----
    This class handles evidence tag assignment only. The combining step
    (determining final Pathogenic/Likely_Pathogenic/VUS classification from
    accumulated tags) is handled separately by the combining module.
    """

    def classify(
        self, variants: Iterator[ScoredVariant]
    ) -> Iterator[ClassifiedVariant]:
        """Assign evidence tags and classify each scored variant.

        Evaluates ACMG/AMP 2015 evidence criteria for each variant, then
        applies combining rules to determine the final classification
        (Pathogenic, Likely_Pathogenic, or VUS).

        Parameters
        ----------
        variants : Iterator[ScoredVariant]
            Stream of scored variants to classify.

        Yields
        ------
        ClassifiedVariant
            Each variant with evidence tags assigned, classification
            determined by ACMG/AMP 2015 combining rules, and missing
            data sources recorded.
        """
        for variant in variants:
            tags, missing_sources = self._assign_tags(variant)
            evidence = frozenset(tags)
            classification = combine_evidence(evidence)
            yield ClassifiedVariant(
                scored=variant,
                evidence_tags=evidence,
                classification=classification,
                missing_data_sources=frozenset(missing_sources),
            )

    def _assign_tags(self, variant: ScoredVariant) -> tuple[set[EvidenceTag], set[str]]:
        """Evaluate all evidence criteria for a single variant.

        Parameters
        ----------
        variant : ScoredVariant
            The variant to evaluate.

        Returns
        -------
        tuple[set[EvidenceTag], set[str]]
            A tuple of (assigned tags, missing data source names).
        """
        tags: set[EvidenceTag] = set()
        missing_sources: set[str] = set()

        self._evaluate_pvs1(variant, tags, missing_sources)
        self._evaluate_pm2(variant, tags, missing_sources)
        self._evaluate_pp3(variant, tags, missing_sources)
        self._evaluate_pp5(variant, tags, missing_sources)

        # Benign criteria
        self._evaluate_ba1(variant, tags, missing_sources)
        self._evaluate_bs1(variant, tags, missing_sources)
        self._evaluate_bp4(variant, tags, missing_sources)
        self._evaluate_bp7(variant, tags, missing_sources)

        return tags, missing_sources

    def _evaluate_pvs1(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PVS1 for null variants or splice-site variants with high SpliceAI.

        PVS1 is assigned unconditionally for NONSENSE or FRAMESHIFT variants.
        For SPLICE_SITE variants, PVS1 is assigned only when SpliceAI > 0.8.
        If SpliceAI data is unavailable for a SPLICE_SITE variant, it is
        recorded as a missing source.

        Parameters
        ----------
        variant : ScoredVariant
            The variant to evaluate.
        tags : set[EvidenceTag]
            Accumulator for assigned tags (mutated in place).
        missing_sources : set[str]
            Accumulator for missing data sources (mutated in place).
        """
        consequence = variant.annotated.consequence

        if consequence in _PVS1_CONSEQUENCES:
            tags.add(EvidenceTag.PVS1)
            return

        if consequence == FunctionalConsequence.SPLICE_SITE:
            spliceai = variant.spliceai_score
            if spliceai is None:
                missing_sources.add("SpliceAI")
                return
            if spliceai > _PVS1_SPLICEAI_THRESHOLD:
                tags.add(EvidenceTag.PVS1)

    def _evaluate_pm2(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PM2 if allele frequency is below 0.0001 in ALL populations.

        Uses population-specific frequencies when available. If any
        population exceeds the threshold, PM2 does not fire (the variant
        is not truly rare). Falls back to global AF when per-population
        data is absent. Treats all-None population data as missing.
        """
        annotated = variant.annotated
        pop_freq = annotated.population_frequencies

        if pop_freq is not None:
            # Guard: if all population fields are None, treat as missing data
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
                missing_sources.add("gnomAD")
                return
            if pop_freq.all_below(_PM2_AF_THRESHOLD):
                tags.add(EvidenceTag.PM2)
            return

        # Fallback: global AF
        af = annotated.allele_frequency
        if af is None:
            missing_sources.add("gnomAD")
            return

        if af < _PM2_AF_THRESHOLD:
            tags.add(EvidenceTag.PM2)

    def _evaluate_pp3(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PP3 based on REVEL or SpliceAI computational evidence.

        PP3 is assigned when REVEL > 0.7, or when SpliceAI > 0.5 on a
        splice-adjacent variant (SPLICE_SITE or MISSENSE consequence).
        When neither predictor is available, both are recorded as missing.

        Parameters
        ----------
        variant : ScoredVariant
            The variant to evaluate.
        tags : set[EvidenceTag]
            Accumulator for assigned tags (mutated in place).
        missing_sources : set[str]
            Accumulator for missing data sources (mutated in place).
        """
        revel = variant.revel_score
        spliceai = variant.spliceai_score
        consequence = variant.annotated.consequence

        revel_available = revel is not None
        spliceai_available = spliceai is not None

        if not revel_available and not spliceai_available:
            missing_sources.add("REVEL")
            missing_sources.add("SpliceAI")
            return

        if revel is not None and revel > _PP3_REVEL_THRESHOLD:
            tags.add(EvidenceTag.PP3)
            return

        splice_adjacent = consequence in _PP3_SPLICE_ADJACENT
        if (
            spliceai is not None
            and spliceai > _PP3_SPLICEAI_THRESHOLD
            and splice_adjacent
        ):
            tags.add(EvidenceTag.PP3)
            return

        if not revel_available:
            missing_sources.add("REVEL")
        if not spliceai_available:
            missing_sources.add("SpliceAI")

    def _evaluate_pp5(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PP5 if ClinVar asserts Pathogenic without conflicts.

        PP5 is assigned when the ClinVar assertion is Pathogenic and
        there is no conflicting Benign or Likely_Benign assertion. If
        ClinVar data is unavailable (clinvar_unknown is True and
        assertion is None), PP5 is omitted and ClinVar is recorded as
        a missing data source.

        Parameters
        ----------
        variant : ScoredVariant
            The variant to evaluate.
        tags : set[EvidenceTag]
            Accumulator for assigned tags (mutated in place).
        missing_sources : set[str]
            Accumulator for missing data sources (mutated in place).
        """
        annotated = variant.annotated
        assertion = annotated.clinvar_assertion

        if assertion is None:
            missing_sources.add("ClinVar")
            return

        if assertion == ClinVarAssertion.PATHOGENIC:
            # In this simplified model, a single ClinVar assertion is stored.
            # If the assertion is Pathogenic, there's no conflicting
            # Benign/Likely_Benign assertion recorded, so PP5 applies.
            # A conflicting assertion would show up as one of the benign
            # categories in the assertion field itself.
            tags.add(EvidenceTag.PP5)
        elif assertion in _PP5_CONFLICTING_ASSERTIONS:
            # The assertion itself is Benign or Likely_Benign, so PP5
            # does not apply (this is the "conflicting" case).
            pass

    def _evaluate_ba1(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        _missing_sources: set[str],
    ) -> None:
        """Assign BA1 if any population AF exceeds 5%.

        BA1 is standalone benign evidence. If population-specific
        frequencies are available, checks each population. Falls back
        to global AF when per-population data is absent.
        """
        annotated = variant.annotated
        pop_freq = annotated.population_frequencies

        if pop_freq is not None:
            if pop_freq.any_exceeds(_BA1_AF_THRESHOLD):
                tags.add(EvidenceTag.BA1)

        else:
            # Fallback to global AF
            af = annotated.allele_frequency
            if af is not None and af > _BA1_AF_THRESHOLD:
                tags.add(EvidenceTag.BA1)

    def _evaluate_bs1(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        _missing_sources: set[str],
    ) -> None:
        """Assign BS1 if any population AF exceeds 1%.

        Only fires when BA1 has not already been assigned (BA1 is
        stronger and subsumes BS1 in the combining rules).
        """
        if EvidenceTag.BA1 in tags:
            return

        annotated = variant.annotated
        pop_freq = annotated.population_frequencies

        if pop_freq is not None:
            if pop_freq.any_exceeds(_BS1_AF_THRESHOLD):
                tags.add(EvidenceTag.BS1)

        else:
            af = annotated.allele_frequency
            if af is not None and af > _BS1_AF_THRESHOLD:
                tags.add(EvidenceTag.BS1)

    def _evaluate_bp4(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        _missing_sources: set[str],
    ) -> None:
        """Assign BP4 for computational benign evidence.

        BP4 fires when:
        - Missense with REVEL < 0.15 (computational evidence suggests no impact)
        - Non-missense with CADD Phred < 10

        Does NOT fire for null variants (frameshift/nonsense) where
        computational predictors are not appropriate for benign evidence.
        """
        consequence = variant.annotated.consequence

        # Null variants should not receive computational benign evidence
        if consequence in (
            FunctionalConsequence.FRAMESHIFT,
            FunctionalConsequence.NONSENSE,
        ):
            return

        if consequence == FunctionalConsequence.MISSENSE:
            revel = variant.revel_score
            if revel is not None and revel < _BP4_REVEL_THRESHOLD:
                tags.add(EvidenceTag.BP4)
        else:
            cadd = variant.cadd_phred
            if cadd is not None and cadd < _BP4_CADD_THRESHOLD:
                tags.add(EvidenceTag.BP4)

    def _evaluate_bp7(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        _missing_sources: set[str],
    ) -> None:
        """Assign BP7 for synonymous variants with no splice impact.

        Fires when the variant is synonymous AND SpliceAI < 0.1
        (no predicted splice disruption).
        """
        if variant.annotated.consequence != FunctionalConsequence.SYNONYMOUS:
            return

        spliceai = variant.spliceai_score
        if spliceai is not None and spliceai < _BP7_SPLICEAI_THRESHOLD:
            tags.add(EvidenceTag.BP7)
        elif spliceai is None:
            # Without SpliceAI, we can't confirm no splice impact
            # BP7 requires negative splice evidence, so don't fire
            pass
