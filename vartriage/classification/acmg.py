"""ACMG/AMP evidence tag assignment for variant classification.

This module implements the ACMGClassifier, which evaluates scored variants
against ACMG/AMP 2015 evidence criteria and assigns appropriate evidence tags.
The classifier handles missing data gracefully by omitting tags when required
data sources are unavailable and recording which sources were missing.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING

from vartriage.classification.combining import (
    combine_evidence,
    has_conflicting_evidence,
)
from vartriage.models.variant import (
    ClassifiedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
)

if TYPE_CHECKING:
    from vartriage.annotation.clinvar_protein_index import ClinVarProteinIndex

_PVS1_CONSEQUENCES: frozenset[FunctionalConsequence] = frozenset(
    {
        FunctionalConsequence.NONSENSE,
        FunctionalConsequence.FRAMESHIFT,
    }
)

_PM2_AF_THRESHOLD: float = 0.0001

# ClinGen-calibrated PP3 thresholds (Pejaver et al., 2022)
# Supporting level: REVEL > 0.644
# Moderate level: REVEL > 0.773
_PP3_REVEL_THRESHOLD: float = 0.644

_PP3_REVEL_MODERATE_THRESHOLD: float = 0.773

_PP3_SPLICEAI_THRESHOLD: float = 0.5

_PVS1_SPLICEAI_THRESHOLD: float = 0.8

_BA1_AF_THRESHOLD: float = 0.05

_BS1_AF_THRESHOLD: float = 0.01

# ClinGen-calibrated BP4 thresholds (Pejaver et al., 2022)
# Supporting level: REVEL < 0.290
# Moderate level: REVEL < 0.183
_BP4_REVEL_THRESHOLD: float = 0.290

_BP4_REVEL_MODERATE_THRESHOLD: float = 0.183

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

_PM4_CONSEQUENCES: frozenset[FunctionalConsequence] = frozenset(
    {
        FunctionalConsequence.IN_FRAME_INSERTION,
        FunctionalConsequence.IN_FRAME_DELETION,
        FunctionalConsequence.STOP_LOSS,
    }
)

_MISSING_SOURCE_GNOMAD_CONSTRAINT = "gnomAD_constraint"


class ACMGClassifier:
    """Assign ACMG/AMP evidence tags and final classification.

    The classifier evaluates each ScoredVariant against ten evidence criteria:

    Pathogenic:
    - PVS1: Nonsense or Frameshift consequence (null variant)
    - PS1: Same amino acid change as established pathogenic (different nucleotide)
    - PM2: gnomAD allele frequency below 0.0001 (absent from controls)
    - PM5: Novel missense at amino acid position with known pathogenic change
    - PP3: REVEL score above threshold (computational evidence)
    - PP5: ClinVar Pathogenic with no conflicting Benign/Likely_Benign

    Benign:
    - BA1: Any population AF > 5% (standalone benign)
    - BS1: Any population AF > 1%
    - BP4: Low computational pathogenicity score
    - BP7: Synonymous with no splice impact

    When a required data source is unavailable for a given criterion, that
    tag is omitted and the source name is recorded in the output.
    """

    def __init__(
        self,
        protein_index: ClinVarProteinIndex | None = None,
        lof_gene_list: frozenset[str] | None = None,
    ) -> None:
        """Initialize the classifier with optional protein-level ClinVar index.

        Parameters
        ----------
        protein_index : Optional[ClinVarProteinIndex]
            Pre-loaded index of ClinVar pathogenic missense variants for
            PS1/PM5 evaluation. When None, PS1 and PM5 are omitted with
            the source recorded as missing.
        lof_gene_list : Optional[frozenset[str]]
            Explicit set of gene symbols where LoF is the established
            disease mechanism. When provided, PVS1 fires at Very Strong
            only for genes on this list. Genes not on the list get PVS1
            at Strong (downgraded). When None, pLI-based gating is used.
            Gene names are matched case-sensitively against
            ``variant.annotated.gene_name`` (HGNC symbols, e.g., "BRCA1").
        """
        self._protein_index = protein_index
        self._lof_gene_list = lof_gene_list

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
                has_conflicting_evidence=has_conflicting_evidence(evidence),
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
        self._evaluate_ps1(variant, tags, missing_sources)
        self._evaluate_pm1(variant, tags, missing_sources)
        self._evaluate_pm2(variant, tags, missing_sources)
        self._evaluate_pm4(variant, tags, missing_sources)
        self._evaluate_pm5(variant, tags, missing_sources)
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
        """Assign PVS1 for null variants with LoF constraint gating.

        For NONSENSE or FRAMESHIFT: fires PVS1 at Very Strong when LoF
        is the established mechanism (pLI > 0.9 or gene on lof_gene_list),
        at Strong otherwise. Novel genes without constraint data get
        benefit of the doubt (Very Strong).

        For SPLICE_SITE: fires PVS1 at Very Strong only when SpliceAI > 0.8.
        """
        consequence = variant.annotated.consequence

        if consequence in _PVS1_CONSEQUENCES:
            tag = self._resolve_pvs1_strength(variant)
            tags.add(tag)
            return

        if consequence == FunctionalConsequence.SPLICE_SITE:
            spliceai = variant.spliceai_score
            if spliceai is None:
                missing_sources.add("SpliceAI")
                return
            if spliceai > _PVS1_SPLICEAI_THRESHOLD:
                tags.add(EvidenceTag.PVS1)

    def _resolve_pvs1_strength(self, variant: ScoredVariant) -> EvidenceTag:
        """Determine PVS1 strength based on gene LoF mechanism evidence.

        Priority:
        1. Explicit lof_gene_list (gene on list → VS, not on list → Strong)
        2. gnomAD pLI constraint (pLI > 0.9 → VS, otherwise → Strong)
        3. No data available → VS (benefit of doubt for novel genes)
        """
        gene_name = variant.annotated.gene_name

        # Explicit gene list takes priority when provided
        if self._lof_gene_list is not None and gene_name is not None:
            if gene_name in self._lof_gene_list:
                return EvidenceTag.PVS1
            return EvidenceTag.PVS1_STRONG

        # Fall back to gnomAD pLI constraint
        gene_context = variant.annotated.gene_context
        if gene_context is None or gene_context.constraint is None:
            return EvidenceTag.PVS1

        if gene_context.constraint.is_lof_intolerant:
            return EvidenceTag.PVS1

        return EvidenceTag.PVS1_STRONG

    def _evaluate_ps1(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PS1 for same amino acid change as established pathogenic variant.

        PS1 fires when a different nucleotide change at the same codon
        produces the same amino acid substitution as a known ClinVar
        Pathogenic variant. Requires both the protein index and protein
        change annotation on the variant.
        """
        # PS1 only applies to missense variants
        if variant.annotated.consequence != FunctionalConsequence.MISSENSE:
            return

        protein_change = variant.annotated.protein_change
        if protein_change is None:
            # Missense but no codon resolution (no reference FASTA) — can't evaluate
            missing_sources.add("codon_resolution")
            return

        if self._protein_index is None or not self._protein_index.is_loaded:
            missing_sources.add("ClinVar_protein_index")
            return

        v = variant.annotated.variant
        if self._protein_index.check_ps1(
            gene=protein_change.gene_name,
            aa_position=protein_change.position,
            ref_aa=protein_change.reference_aa,
            alt_aa=protein_change.altered_aa,
            chrom=v.chrom,
            genomic_pos=v.pos,
            ref_allele=v.ref,
            alt_allele=v.alt,
        ):
            tags.add(EvidenceTag.PS1)

    def _evaluate_pm5(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PM5 for novel missense at position with known pathogenic change.

        PM5 fires when the variant introduces a different amino acid change
        at a position where another missense change is already classified
        as Pathogenic. Does not fire if PS1 already assigned (PS1 is stronger
        and the same-change case subsumes the different-change case).
        """
        # PM5 only applies to missense variants
        if variant.annotated.consequence != FunctionalConsequence.MISSENSE:
            return

        protein_change = variant.annotated.protein_change
        if protein_change is None:
            missing_sources.add("codon_resolution")
            return

        if self._protein_index is None or not self._protein_index.is_loaded:
            missing_sources.add("ClinVar_protein_index")
            return

        # Don't double-count: if PS1 already fired, PM5 is redundant
        if EvidenceTag.PS1 in tags:
            return

        if self._protein_index.check_pm5(
            gene=protein_change.gene_name,
            aa_position=protein_change.position,
            ref_aa=protein_change.reference_aa,
            alt_aa=protein_change.altered_aa,
        ):
            tags.add(EvidenceTag.PM5)

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
        """Assign PP3 based on ClinGen-calibrated REVEL or SpliceAI thresholds.

        Strength-modulated per Pejaver et al. (2022):
        - PP3_Moderate: REVEL > 0.773
        - PP3 (supporting): REVEL > 0.644
        - PP3 (supporting): SpliceAI > 0.5 on splice-adjacent variant

        Only the highest applicable strength fires. When neither predictor
        is available, both are recorded as missing.
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

        # Check REVEL at moderate threshold first (higher bar = stronger evidence)
        if revel is not None and revel > _PP3_REVEL_MODERATE_THRESHOLD:
            tags.add(EvidenceTag.PP3_MODERATE)
            return

        # Then supporting-level REVEL
        if revel is not None and revel > _PP3_REVEL_THRESHOLD:
            tags.add(EvidenceTag.PP3)
            return

        # SpliceAI-based PP3 (supporting only)
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

        Strength-modulated per Pejaver et al. (2022):
        - BP4_Moderate: REVEL < 0.183 (stronger benign evidence)
        - BP4 (supporting): REVEL < 0.290

        For non-missense variants, CADD Phred < 10 triggers supporting BP4.
        Does NOT fire for protein-altering variants (null variants or in-frame
        indels) where low computational scores are not appropriate evidence of
        benign impact. In-frame indels already fire PM4; awarding BP4
        simultaneously would create contradictory evidence for the same variant.
        """
        consequence = variant.annotated.consequence

        # Null and protein-length-altering variants should not receive
        # computational benign evidence. In-frame indels get PM4 (moderate
        # pathogenic); BP4 on the same variant is logically contradictory.
        if consequence in (
            FunctionalConsequence.FRAMESHIFT,
            FunctionalConsequence.NONSENSE,
            FunctionalConsequence.STOP_LOSS,
            FunctionalConsequence.IN_FRAME_INSERTION,
            FunctionalConsequence.IN_FRAME_DELETION,
        ):
            return

        if consequence == FunctionalConsequence.MISSENSE:
            revel = variant.revel_score
            if revel is not None:
                # Moderate level first (stricter threshold = more confident benign)
                if revel < _BP4_REVEL_MODERATE_THRESHOLD:
                    tags.add(EvidenceTag.BP4_MODERATE)
                elif revel < _BP4_REVEL_THRESHOLD:
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

    def _evaluate_pm1(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PM1 for missense in a functionally constrained gene region.

        Uses gnomAD missense constraint (mis_z > 3.09) as a proxy for
        functional domain intolerance. Only fires for missense variants.
        """
        if variant.annotated.consequence != FunctionalConsequence.MISSENSE:
            return

        gene_context = variant.annotated.gene_context
        if gene_context is None or gene_context.constraint is None:
            missing_sources.add(_MISSING_SOURCE_GNOMAD_CONSTRAINT)
            return

        if gene_context.constraint.is_missense_constrained:
            tags.add(EvidenceTag.PM1)

    def _evaluate_pm4(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        _missing_sources: set[str],
    ) -> None:
        """Assign PM4 for protein-length-changing variants.

        Fires for in-frame insertions, in-frame deletions, and
        stop-loss variants. These alter the protein without truncating
        the reading frame.
        """
        if variant.annotated.consequence in _PM4_CONSEQUENCES:
            tags.add(EvidenceTag.PM4)
