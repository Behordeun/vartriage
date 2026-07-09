"""ACMG/AMP evidence tag assignment for variant classification.

This module implements the ACMGClassifier, which evaluates scored variants
against ACMG/AMP 2015 evidence criteria and assigns appropriate evidence tags.
The classifier handles missing data gracefully by omitting tags when required
data sources are unavailable and recording which sources were missing.
"""

from __future__ import annotations

from typing import Iterator

from vartriage.classification.combining import combine_evidence
from vartriage.models.variant import (
    ClassifiedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
)

_PVS1_CONSEQUENCES: frozenset[FunctionalConsequence] = frozenset({
    FunctionalConsequence.NONSENSE,
    FunctionalConsequence.FRAMESHIFT,
})

_PM2_AF_THRESHOLD: float = 0.0001

_PP3_REVEL_THRESHOLD: float = 0.7

_PP5_CONFLICTING_ASSERTIONS: frozenset[ClinVarAssertion] = frozenset({
    ClinVarAssertion.BENIGN,
    ClinVarAssertion.LIKELY_BENIGN,
})


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

    def _assign_tags(
        self, variant: ScoredVariant
    ) -> tuple[set[EvidenceTag], set[str]]:
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

        self._evaluate_pvs1(variant, tags)
        self._evaluate_pm2(variant, tags, missing_sources)
        self._evaluate_pp3(variant, tags, missing_sources)
        self._evaluate_pp5(variant, tags, missing_sources)

        return tags, missing_sources

    def _evaluate_pvs1(
        self, variant: ScoredVariant, tags: set[EvidenceTag]
    ) -> None:
        """Assign PVS1 if the variant consequence is Nonsense or Frameshift.

        Parameters
        ----------
        variant : ScoredVariant
            The variant to evaluate.
        tags : set[EvidenceTag]
            Accumulator for assigned tags (mutated in place).
        """
        consequence = variant.annotated.consequence
        if consequence in _PVS1_CONSEQUENCES:
            tags.add(EvidenceTag.PVS1)

    def _evaluate_pm2(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PM2 if allele frequency is below 0.0001.

        If gnomAD frequency data is unavailable (frequency_unknown is True and
        allele_frequency is None), PM2 is omitted and gnomAD is recorded as a
        missing data source.

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
        """Assign PP3 if REVEL score exceeds 0.7.

        If REVEL score is unavailable, PP3 is omitted and REVEL is recorded
        as a missing data source.

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

        if revel is None:
            missing_sources.add("REVEL")
            return

        if revel > _PP3_REVEL_THRESHOLD:
            tags.add(EvidenceTag.PP3)

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
            # Benign/Likely_Benign assertion recorded — PP5 applies.
            # A conflicting assertion would show up as one of the benign
            # categories in the assertion field itself.
            tags.add(EvidenceTag.PP5)
        elif assertion in _PP5_CONFLICTING_ASSERTIONS:
            # The assertion itself is Benign or Likely_Benign, so PP5
            # does not apply (this is the "conflicting" case).
            pass
