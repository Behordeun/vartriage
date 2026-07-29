"""Structural variant pathogenicity scoring.

Computes a composite pathogenicity score for each annotated SV based on:
- Gene impact severity (consequence type weight)
- Dosage sensitivity (ClinGen HI/TS scores of affected genes)
- Population frequency rarity (absent from gnomAD-SV = high score)
- Size relative to gene content (multi-gene SVs scored higher)

The composite is a weighted sum normalized to [0.0, 1.0]. SVs with no
gene overlap and no frequency data receive None (unscoreable).
"""

from __future__ import annotations

import logging
from typing import Iterator, Optional

from vartriage.structural.models import (
    AnnotatedSV,
    ScoredSV,
    SVConsequence,
    SVType,
)

logger = logging.getLogger(__name__)

# Weight distribution for the composite score components
_IMPACT_WEIGHT: float = 0.35
_DOSAGE_WEIGHT: float = 0.30
_FREQUENCY_WEIGHT: float = 0.20
_SIZE_WEIGHT: float = 0.15

# Base scores by consequence type (gene impact severity)
_CONSEQUENCE_BASE_SCORES: dict[SVConsequence, float] = {
    SVConsequence.WHOLE_GENE_DELETION: 1.0,
    SVConsequence.PARTIAL_GENE_DELETION: 0.7,
    SVConsequence.WHOLE_GENE_DUPLICATION: 0.8,
    SVConsequence.PARTIAL_GENE_DUPLICATION: 0.5,
    SVConsequence.GENE_DISRUPTION: 0.75,
    SVConsequence.INTRONIC: 0.2,
    SVConsequence.REGULATORY: 0.15,
    SVConsequence.INTERGENIC: 0.0,
}

# Size thresholds for the size score component (in bp)
_SIZE_LARGE: int = 1_000_000
_SIZE_MEDIUM: int = 100_000
_SIZE_SMALL: int = 10_000


class SVScorer:
    """Score annotated structural variants by pathogenicity.

    Processes a stream of AnnotatedSV records and attaches a composite
    pathogenicity score based on gene impact, dosage sensitivity,
    population frequency, and SV size.

    The scorer applies allele frequency filtering before scoring:
    SVs with population frequency above max_af are excluded from
    output (they are common and presumed benign).

    Parameters
    ----------
    max_allele_frequency : float
        Maximum gnomAD-SV frequency. SVs above this threshold are
        excluded. Default is 0.01 (1%).
    """

    def __init__(self, max_allele_frequency: float = 0.01) -> None:
        self._max_af = max_allele_frequency

    def score(self, variants: Iterator[AnnotatedSV]) -> Iterator[ScoredSV]:
        """Score and filter a stream of annotated SVs.

        Parameters
        ----------
        variants : Iterator[AnnotatedSV]
            Annotated SVs from the gene annotator.

        Yields
        ------
        ScoredSV
            SVs that pass frequency filtering, sorted within each
            batch by pathogenicity_score descending.
        """
        for annotated in variants:
            if not self._passes_frequency_filter(annotated):
                continue
            yield self._score_single(annotated)

    def _passes_frequency_filter(self, sv: AnnotatedSV) -> bool:
        """Exclude common SVs seen frequently in the population."""
        if sv.population_frequency is None:
            return True
        return sv.population_frequency <= self._max_af

    def _score_single(self, annotated: AnnotatedSV) -> ScoredSV:
        """Compute composite pathogenicity score for one SV."""
        impact = self._compute_impact_score(annotated)
        dosage = self._compute_dosage_score(annotated)
        frequency = self._compute_frequency_score(annotated)
        size = self._compute_size_score(annotated)

        # Intergenic SVs with no frequency data are unscoreable
        if (
            annotated.consequence == SVConsequence.INTERGENIC
            and annotated.frequency_unknown
        ):
            return ScoredSV(
                annotated=annotated,
                pathogenicity_score=None,
                dosage_score=dosage,
                size_score=size,
                frequency_score=frequency,
            )

        composite = (
            impact * _IMPACT_WEIGHT
            + dosage * _DOSAGE_WEIGHT
            + frequency * _FREQUENCY_WEIGHT
            + size * _SIZE_WEIGHT
        )

        # Clamp to [0.0, 1.0]
        composite = max(0.0, min(1.0, composite))

        return ScoredSV(
            annotated=annotated,
            pathogenicity_score=composite,
            dosage_score=dosage,
            size_score=size,
            frequency_score=frequency,
        )

    def _compute_impact_score(self, sv: AnnotatedSV) -> float:
        """Score based on the most severe gene-level consequence.

        Multi-gene deletions get a boost: deleting 3+ HI genes is
        worse than deleting 1.
        """
        base = _CONSEQUENCE_BASE_SCORES.get(sv.consequence, 0.0)

        # Boost for multi-gene events
        if sv.genes_affected > 1:
            gene_boost = min(0.2, sv.genes_affected * 0.05)
            base = min(1.0, base + gene_boost)

        return base

    def _compute_dosage_score(self, sv: AnnotatedSV) -> float:
        """Score based on dosage sensitivity of affected genes.

        Uses the highest HI or TS score among overlapped genes,
        scaled to [0.0, 1.0]. HI applies to losses (DEL/CNV<2),
        TS applies to gains (DUP/CNV>2).
        """
        if not sv.gene_overlaps:
            return 0.0

        is_loss = sv.sv.sv_type in (SVType.DEL,) or (
            sv.sv.sv_type == SVType.CNV
            and sv.sv.copy_number is not None
            and sv.sv.copy_number < 2
        )

        best_score = 0.0

        for overlap in sv.gene_overlaps:
            if is_loss:
                # HI score: 3 = sufficient evidence, scale 0-3 to 0-1
                if overlap.hi_score is not None:
                    normalized = min(1.0, overlap.hi_score / 3.0)
                    best_score = max(best_score, normalized)
            else:
                # TS score for gains
                if overlap.ts_score is not None:
                    normalized = min(1.0, overlap.ts_score / 3.0)
                    best_score = max(best_score, normalized)

        # If no dosage data available, use a modest default for
        # protein-coding gene overlap (some pathogenicity assumed)
        if best_score == 0.0 and sv.genes_affected > 0:
            best_score = 0.3

        return best_score

    def _compute_frequency_score(self, sv: AnnotatedSV) -> float:
        """Score based on population rarity.

        Absent from gnomAD-SV = 1.0 (rare, potentially pathogenic).
        Very common = 0.0 (likely benign, but these are already
        filtered by max_af so this mainly distinguishes among rare SVs).
        """
        if sv.frequency_unknown or sv.population_frequency is None:
            return 1.0

        af = sv.population_frequency

        # Linear decay: AF=0 → 1.0, AF=max_af → 0.0
        if self._max_af > 0:
            return max(0.0, 1.0 - (af / self._max_af))
        return 1.0

    def _compute_size_score(self, sv: AnnotatedSV) -> float:
        """Score based on SV size relative to clinical significance.

        Larger SVs affecting more genomic content are generally more
        likely to be pathogenic. The relationship is logarithmic:
        a 1Mb deletion isn't 10x worse than a 100kb deletion.
        """
        length = sv.sv.length

        if length >= _SIZE_LARGE:
            return 1.0
        if length >= _SIZE_MEDIUM:
            # Scale 100kb-1Mb to 0.6-1.0
            frac = (length - _SIZE_MEDIUM) / (_SIZE_LARGE - _SIZE_MEDIUM)
            return 0.6 + frac * 0.4
        if length >= _SIZE_SMALL:
            # Scale 10kb-100kb to 0.3-0.6
            frac = (length - _SIZE_SMALL) / (_SIZE_MEDIUM - _SIZE_SMALL)
            return 0.3 + frac * 0.3

        # Below 10kb: scale 0-10kb to 0.0-0.3
        frac = length / _SIZE_SMALL
        return frac * 0.3
