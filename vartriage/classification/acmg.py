"""ACMG/AMP evidence tag assignment for variant classification.

This module implements the ACMGClassifier, which evaluates scored variants
against ACMG/AMP 2015 evidence criteria and assigns appropriate evidence tags.
The classifier handles missing data gracefully by omitting tags when required
data sources are unavailable and recording which sources were missing.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING, Protocol

from vartriage.classification.combining import (
    combine_evidence,
    has_conflicting_evidence,
)
from vartriage.models.config import (
    DEFAULT_THRESHOLDS,
    DOMINANT_THRESHOLDS,
    RECESSIVE_THRESHOLDS,
    DiseaseThresholds,
)
from vartriage.models.variant import (
    ClassifiedVariant,
    ClinVarAssertion,
    ClinVarReviewStatus,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
)

if TYPE_CHECKING:
    from vartriage.annotation.clinvar_protein_index import ClinVarProteinIndex
    from vartriage.annotation.transcript_index import NMDEscapeZone

    class NMDLookup(Protocol):
        """Minimal interface the classifier needs for NMD-escape evaluation."""

        def escape_zone(
            self, gene_name: str, chrom: str, pos: int
        ) -> NMDEscapeZone | None: ...


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
# Strong level: REVEL > 0.932 (Pejaver et al. 2022 calibrated threshold)
_PP3_REVEL_STRONG_THRESHOLD: float = 0.932

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

_MISSING_SOURCE_GNOMAD = "gnomAD"
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
        nmd_lookup: NMDLookup | None = None,
        use_disease_thresholds: bool = False,
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
        nmd_lookup : Optional[NMDLookup]
            Transcript CDS index exposing ``escape_zone(gene, chrom, pos)``.
            When provided, PVS1 is downgraded from Very Strong to Strong for
            null variants that escape nonsense-mediated decay (last exon,
            within 50 nt of the final junction, or a single-exon gene). When
            None, ``transcript_structure`` is recorded as missing and PVS1
            keeps its current strength.
        use_disease_thresholds : bool
            When True, BA1/BS1/PM2 frequency gates are selected per variant
            from the gene's inheritance mode (via ``gene_context``). When
            False (default), fixed thresholds are used for backward
            compatibility.
        """
        self._protein_index = protein_index
        self._lof_gene_list = lof_gene_list
        self._nmd_lookup = nmd_lookup
        self._use_disease_thresholds = use_disease_thresholds

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
        self._evaluate_bs2(variant, tags, missing_sources)
        self._evaluate_bp4(variant, tags, missing_sources)
        self._evaluate_bp7(variant, tags, missing_sources)

        return tags, missing_sources

    def _resolve_thresholds(self, variant: ScoredVariant) -> DiseaseThresholds:
        """Select the BA1/BS1/PM2 frequency thresholds for a variant.

        When disease-aware thresholds are disabled, or the variant has no
        gene context with an inheritance mode, the fixed default thresholds
        are returned (backward compatible). Otherwise the gene's inheritance
        mode selects dominant (stricter) or recessive thresholds. The most
        clinically significant mode across the gene's disease associations
        wins: a dominant association tightens the gates even if a recessive
        one is also recorded.
        """
        if not self._use_disease_thresholds:
            return DEFAULT_THRESHOLDS

        gene_context = variant.annotated.gene_context
        if gene_context is None or not gene_context.disease_associations:
            return DEFAULT_THRESHOLDS

        modes = {
            a.inheritance_mode.upper()
            for a in gene_context.disease_associations
            if a.inheritance_mode
        }
        if not modes:
            return DEFAULT_THRESHOLDS
        if "AD" in modes:
            return DOMINANT_THRESHOLDS
        if modes & {"AR", "XL", "XLR", "XLD", "MT"}:
            return RECESSIVE_THRESHOLDS
        return DEFAULT_THRESHOLDS

    def _evaluate_pvs1(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PVS1 for null variants with LoF constraint gating.

        For NONSENSE or FRAMESHIFT: fires PVS1 at Very Strong when LoF
        is the established mechanism (pLI > 0.9 or gene on lof_gene_list),
        at Strong otherwise, including when the mechanism is unknown for
        want of constraint data.

        For SPLICE_SITE: fires PVS1 at Very Strong only when SpliceAI > 0.8.
        """
        consequence = variant.annotated.consequence

        if consequence in _PVS1_CONSEQUENCES:
            tag = self._resolve_pvs1_strength(variant)
            tags.add(tag)
            # When PVS1 stands at Very Strong but no transcript structure was
            # available to check NMD escape, record the source so the output
            # reflects that the downgrade check could not run.
            if tag is EvidenceTag.PVS1 and self._nmd_lookup is None:
                missing_sources.add("transcript_structure")
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

        PVS1 at Very Strong is warranted only when loss of function is an
        established disease mechanism for the gene. Priority:
        1. Explicit lof_gene_list (gene on list -> Very Strong, off -> Strong)
        2. gnomAD pLI constraint (pLI > 0.9 -> Very Strong, otherwise Strong)
        3. No mechanism evidence -> Strong (the mechanism is unknown, so the
           criterion does not reach Very Strong on absence of data).

        A Very Strong result is then downgraded to Strong when the variant
        escapes nonsense-mediated decay (NMD), because such a truncating
        variant may still produce a partially functional protein. The
        downgrade only ever lowers strength; it never raises it.
        """
        gene_name = variant.annotated.gene_name

        # Explicit gene list takes priority when provided
        if self._lof_gene_list is not None and gene_name is not None:
            if gene_name in self._lof_gene_list:
                return self._apply_nmd_downgrade(variant, EvidenceTag.PVS1)
            return EvidenceTag.PVS1_STRONG

        # Fall back to gnomAD pLI constraint
        gene_context = variant.annotated.gene_context
        if gene_context is None or gene_context.constraint is None:
            return EvidenceTag.PVS1_STRONG

        if gene_context.constraint.is_lof_intolerant:
            return self._apply_nmd_downgrade(variant, EvidenceTag.PVS1)

        return EvidenceTag.PVS1_STRONG

    def _apply_nmd_downgrade(
        self, variant: ScoredVariant, resolved: EvidenceTag
    ) -> EvidenceTag:
        """Downgrade a Very Strong PVS1 to Strong when the variant escapes NMD.

        Applies only to ``EvidenceTag.PVS1`` (Very Strong). A variant in the
        last exon, within 50 nt of the final exon-exon junction, or in a
        single-exon gene escapes nonsense-mediated decay and is downgraded to
        ``PVS1_STRONG``. When no NMD lookup is configured the strength is
        unchanged; the missing source is recorded by the caller via
        ``_evaluate_pvs1``.
        """
        if resolved is not EvidenceTag.PVS1 or self._nmd_lookup is None:
            return resolved

        gene_name = variant.annotated.gene_name
        if gene_name is None:
            return resolved

        v = variant.annotated.variant
        zone = self._nmd_lookup.escape_zone(gene_name, v.chrom, v.pos)
        if zone is None:
            return resolved

        if (
            zone.is_single_exon
            or zone.is_in_last_exon(v.pos)
            or zone.is_near_last_junction(v.pos, margin=50)
        ):
            return EvidenceTag.PVS1_STRONG
        return resolved

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
        """Assign PM2 when the variant is rare in consulted population data.

        PM2 encodes "absent from, or rare in, population controls", which
        requires the population database to have been consulted. Two inputs
        satisfy that: an observed allele frequency below the threshold in
        every available population, or ``frequency_unknown`` set, which the
        annotation stage sets when gnomAD was queried and returned no record
        for the variant.

        A variant with no frequency data and ``frequency_unknown`` not set
        was never consulted (annotation not run, or the lookup did not
        complete). That is the absence of evidence, not evidence of rarity,
        so PM2 does not fire and gnomAD is recorded as a missing source.
        """
        annotated = variant.annotated
        pop_freq = annotated.population_frequencies
        pm2_threshold = self._resolve_thresholds(variant).pm2_af

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
            if has_any_data:
                if pop_freq.all_below(pm2_threshold):
                    tags.add(EvidenceTag.PM2)
                return
            # Every per-population field is None. Only a confirmed gnomAD
            # miss (frequency_unknown) counts as absence; otherwise the data
            # was never obtained.
            if annotated.frequency_unknown:
                tags.add(EvidenceTag.PM2)
            else:
                missing_sources.add(_MISSING_SOURCE_GNOMAD)
            return

        af = annotated.allele_frequency
        if af is not None:
            if af < pm2_threshold:
                tags.add(EvidenceTag.PM2)
            return

        # No observed AF. A confirmed gnomAD miss is absence; a bare None
        # with no such confirmation is missing data.
        if annotated.frequency_unknown:
            tags.add(EvidenceTag.PM2)
        else:
            missing_sources.add(_MISSING_SOURCE_GNOMAD)

    def _evaluate_pp3(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign PP3 based on ClinGen-calibrated REVEL or SpliceAI thresholds.

        Strength-modulated per Pejaver et al. (2022):
        - PP3_Strong: REVEL > 0.932
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

        # Check REVEL at strong threshold first (highest bar = strongest evidence)
        if revel is not None and revel > _PP3_REVEL_STRONG_THRESHOLD:
            tags.add(EvidenceTag.PP3_STRONG)
            return

        # Then moderate threshold
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
            # Modulate PP5 strength by review status when the extended ClinVar
            # format supplied it. Without a review status (basic format), keep
            # the single-assertion supporting-level behavior.
            review_status = annotated.clinvar_review_status
            if review_status == ClinVarReviewStatus.EXPERT_PANEL:
                tags.add(EvidenceTag.PP5_STRONG)
            elif review_status in (
                ClinVarReviewStatus.NO_CRITERIA,
                ClinVarReviewStatus.CONFLICTING,
            ):
                # No assertion criteria (or conflicting) does not support PP5.
                return
            else:
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
        ba1_threshold = self._resolve_thresholds(variant).ba1_af

        if pop_freq is not None:
            if pop_freq.any_exceeds(ba1_threshold):
                tags.add(EvidenceTag.BA1)

        else:
            # Fallback to global AF
            af = annotated.allele_frequency
            if af is not None and af > ba1_threshold:
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
        bs1_threshold = self._resolve_thresholds(variant).bs1_af

        if pop_freq is not None:
            if pop_freq.any_exceeds(bs1_threshold):
                tags.add(EvidenceTag.BS1)

        else:
            af = annotated.allele_frequency
            if af is not None and af > bs1_threshold:
                tags.add(EvidenceTag.BS1)

    def _evaluate_bs2(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign BS2 for a variant observed homozygous in healthy controls.

        BS2 is strong benign evidence for a fully penetrant dominant
        disorder: a homozygous observation in gnomAD argues against
        pathogenicity. It fires only when the gnomAD homozygote count is
        available and positive AND the gene is associated with a dominant
        disorder. It does not fire for recessive disorders, where homozygous
        carriers are expected. When the homozygote count is unavailable, the
        source is recorded as missing rather than treated as absence.
        """
        annotated = variant.annotated
        gene_context = annotated.gene_context

        # BS2 only applies to dominant disorders. Without a dominant
        # association the criterion is not applicable, so the homozygote
        # count is not a required source and nothing is recorded missing.
        if gene_context is None or not gene_context.disease_associations:
            return

        is_dominant = any(
            (a.inheritance_mode or "").upper() == "AD"
            for a in gene_context.disease_associations
        )
        if not is_dominant:
            return

        # The gene is dominant, so BS2 is applicable and the homozygote count
        # is required. Record it as missing when unavailable.
        pop_freq = annotated.population_frequencies
        if pop_freq is None or pop_freq.hom_count is None:
            missing_sources.add("gnomAD_homozygotes")
            return

        if pop_freq.hom_count > 0:
            tags.add(EvidenceTag.BS2)

    def _evaluate_bp4(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        """Assign BP4 for computational benign evidence.

        Strength-modulated per Pejaver et al. (2022):
        - BP4_Moderate: REVEL < 0.183 (stronger benign evidence)
        - BP4 (supporting): REVEL < 0.290

        For missense variants without a REVEL score, a low CADD Phred
        (< 10) supports BP4, mirroring the CADD path used for other
        consequence classes and the fallback PP3 uses on the pathogenic
        side. When neither REVEL nor CADD is available for a missense
        variant, REVEL is recorded as a missing source.

        Does NOT fire for protein-altering variants (null variants or
        in-frame indels) where low computational scores are not
        appropriate evidence of benign impact. In-frame indels already
        fire PM4; awarding BP4 simultaneously would create contradictory
        evidence for the same variant.
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
            self._evaluate_bp4_missense(variant, tags, missing_sources)
        else:
            cadd = variant.cadd_phred
            if cadd is not None and cadd < _BP4_CADD_THRESHOLD:
                tags.add(EvidenceTag.BP4)

    def _evaluate_bp4_missense(
        self,
        variant: ScoredVariant,
        tags: set[EvidenceTag],
        missing_sources: set[str],
    ) -> None:
        revel = variant.revel_score
        if revel is not None:
            if revel < _BP4_REVEL_MODERATE_THRESHOLD:
                tags.add(EvidenceTag.BP4_MODERATE)
            elif revel < _BP4_REVEL_THRESHOLD:
                tags.add(EvidenceTag.BP4)
            return
        cadd = variant.cadd_phred
        if cadd is not None:
            if cadd < _BP4_CADD_THRESHOLD:
                tags.add(EvidenceTag.BP4)
            return
        missing_sources.add("REVEL")

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
