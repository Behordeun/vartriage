"""Mapping from Ensembl VEP Sequence Ontology terms to FunctionalConsequence.

VEP returns fine-grained SO terms (~50 possible values). This module
collapses them into vartriage's simplified 8-value enum using a
severity-ranked lookup. When multiple transcripts overlap a variant,
the most severe consequence wins.

Severity ranking matches Ensembl's own ordering with our enum granularity.
Terms not in the mapping default to INTERGENIC with a logged warning.
"""

from __future__ import annotations

import logging

from vartriage.models.variant import FunctionalConsequence

logger = logging.getLogger(__name__)

# Severity rank: lower number = more severe.
# Each SO term maps to (FunctionalConsequence, severity_rank).
_SO_TERM_MAP: dict[str, tuple[FunctionalConsequence, int]] = {
    # Rank 1: Splice-disrupting (loss of canonical splice signals)
    "transcript_ablation": (FunctionalConsequence.SPLICE_SITE, 1),
    "splice_acceptor_variant": (FunctionalConsequence.SPLICE_SITE, 1),
    "splice_donor_variant": (FunctionalConsequence.SPLICE_SITE, 1),
    "splice_donor_5th_base_variant": (FunctionalConsequence.SPLICE_SITE, 1),
    "splice_donor_region_variant": (FunctionalConsequence.SPLICE_SITE, 1),
    # Rank 2: Nonsense (premature termination)
    "stop_gained": (FunctionalConsequence.NONSENSE, 2),
    "start_lost": (FunctionalConsequence.NONSENSE, 2),
    # Rank 3: Frameshift
    "frameshift_variant": (FunctionalConsequence.FRAMESHIFT, 3),
    # Rank 4: Stop loss / start codon disruption (protein extension)
    "stop_lost": (FunctionalConsequence.MISSENSE, 4),
    "initiator_codon_variant": (FunctionalConsequence.MISSENSE, 4),
    # Rank 5: In-frame insertion
    "inframe_insertion": (FunctionalConsequence.IN_FRAME_INSERTION, 5),
    # Rank 6: In-frame deletion
    "inframe_deletion": (FunctionalConsequence.IN_FRAME_DELETION, 6),
    # Rank 7: Missense
    "missense_variant": (FunctionalConsequence.MISSENSE, 7),
    "protein_altering_variant": (FunctionalConsequence.MISSENSE, 7),
    # Rank 8: Splice region (weaker splice signal disruption)
    "splice_region_variant": (FunctionalConsequence.SPLICE_SITE, 8),
    "splice_polypyrimidine_tract_variant": (FunctionalConsequence.SPLICE_SITE, 8),
    # Rank 9: Incomplete terminal codon
    "incomplete_terminal_codon_variant": (FunctionalConsequence.MISSENSE, 9),
    # Rank 10: Synonymous
    "synonymous_variant": (FunctionalConsequence.SYNONYMOUS, 10),
    "stop_retained_variant": (FunctionalConsequence.SYNONYMOUS, 10),
    "start_retained_variant": (FunctionalConsequence.SYNONYMOUS, 10),
    # Rank 11: Coding sequence (ambiguous effect)
    "coding_sequence_variant": (FunctionalConsequence.SYNONYMOUS, 11),
    # Rank 12: UTR variants
    "5_prime_UTR_variant": (FunctionalConsequence.SYNONYMOUS, 12),
    "3_prime_UTR_variant": (FunctionalConsequence.SYNONYMOUS, 12),
    # Rank 13: Non-coding transcript / intronic
    "mature_miRNA_variant": (FunctionalConsequence.SYNONYMOUS, 13),
    "non_coding_transcript_exon_variant": (FunctionalConsequence.SYNONYMOUS, 13),
    "non_coding_transcript_variant": (FunctionalConsequence.SYNONYMOUS, 13),
    "NMD_transcript_variant": (FunctionalConsequence.SYNONYMOUS, 13),
    "intron_variant": (FunctionalConsequence.SYNONYMOUS, 13),
    # Rank 14: Intergenic / regulatory / distant
    "intergenic_variant": (FunctionalConsequence.INTERGENIC, 14),
    "upstream_gene_variant": (FunctionalConsequence.INTERGENIC, 14),
    "downstream_gene_variant": (FunctionalConsequence.INTERGENIC, 14),
    "regulatory_region_variant": (FunctionalConsequence.INTERGENIC, 14),
    "TF_binding_site_variant": (FunctionalConsequence.INTERGENIC, 14),
    "TFBS_ablation": (FunctionalConsequence.INTERGENIC, 14),
    "regulatory_region_ablation": (FunctionalConsequence.INTERGENIC, 14),
    "regulatory_region_amplification": (FunctionalConsequence.INTERGENIC, 14),
    "feature_elongation": (FunctionalConsequence.INTERGENIC, 14),
    "feature_truncation": (FunctionalConsequence.INTERGENIC, 14),
    "sequence_variant": (FunctionalConsequence.INTERGENIC, 14),
}

# Pre-compute the set for fast membership testing
_KNOWN_TERMS: frozenset[str] = frozenset(_SO_TERM_MAP.keys())


def map_so_term(term: str) -> tuple[FunctionalConsequence, int]:
    """Map a single SO consequence term to FunctionalConsequence with severity.

    Parameters
    ----------
    term
        Sequence Ontology term from VEP (e.g., "missense_variant").

    Returns
    -------
    tuple[FunctionalConsequence, int]
        The mapped consequence and its severity rank (lower = more severe).
        Unmapped terms return (INTERGENIC, 99) with a logged warning.
    """
    result = _SO_TERM_MAP.get(term)
    if result is not None:
        return result

    logger.warning("Unmapped VEP consequence term '%s', defaulting to INTERGENIC", term)
    return (FunctionalConsequence.INTERGENIC, 99)


def most_severe_consequence(terms: list[str]) -> FunctionalConsequence:
    """Select the most severe consequence from a list of SO terms.

    Used when a variant overlaps multiple transcripts and VEP returns
    multiple consequence terms. The term with the lowest severity rank wins.

    Parameters
    ----------
    terms
        List of SO consequence terms from VEP response.

    Returns
    -------
    FunctionalConsequence
        The most severe consequence across all terms.
        Returns INTERGENIC if the list is empty.
    """
    if not terms:
        return FunctionalConsequence.INTERGENIC

    best_consequence = FunctionalConsequence.INTERGENIC
    best_rank = 99

    for term in terms:
        consequence, rank = map_so_term(term)
        if rank < best_rank:
            best_rank = rank
            best_consequence = consequence

    return best_consequence


def map_vep_most_severe(
    most_severe_consequence_str: str | None,
) -> FunctionalConsequence:
    """Map VEP's top-level 'most_severe_consequence' field directly.

    VEP responses include a pre-computed most_severe_consequence at the
    variant level. This is faster than re-computing from transcript-level
    terms when we only need the top-level call.

    Parameters
    ----------
    most_severe_consequence_str
        The most_severe_consequence value from VEP response, or None.

    Returns
    -------
    FunctionalConsequence
        Mapped consequence enum value.
    """
    if most_severe_consequence_str is None:
        return FunctionalConsequence.INTERGENIC
    consequence, _ = map_so_term(most_severe_consequence_str)
    return consequence
