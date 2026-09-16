"""Refine a coding-SNV consequence from its resolved codon.

An interval backend can tell that a single-base substitution lands inside a
CDS, but not what it does to the protein: that needs the codon. A backend
that stops at "coding SNV" therefore calls everything missense, which hides
nonsense (stop-gain) variants that carry the strongest pathogenic evidence.

This refiner takes the base consequence and the resolved codon context and
returns the amino-acid-level consequence. It runs at the engine level, so
every interval backend produces the same call for the same variant.
"""

from __future__ import annotations

from typing import Protocol

from vartriage.models.variant import FunctionalConsequence

# Base calls that a codon context is allowed to refine. A splice-site or
# intergenic call carries information the codon does not, so it is never
# overridden here.
_REFINABLE_BASE = frozenset(
    {
        FunctionalConsequence.MISSENSE,
        FunctionalConsequence.SYNONYMOUS,
    }
)


class _CodonLike(Protocol):
    reference_aa: str
    altered_aa: str
    is_synonymous: bool
    is_nonsense: bool


def refine_coding_snv_consequence(
    base: FunctionalConsequence,
    ctx: _CodonLike | None,
) -> FunctionalConsequence:
    """Return the amino-acid-level consequence for a coding SNV.

    With no codon context the base call is returned unchanged. Otherwise a
    refinable coding base call becomes NONSENSE (stop gained), STOP_LOSS
    (stop removed), SYNONYMOUS (same amino acid), or stays MISSENSE.
    """
    if ctx is None or base not in _REFINABLE_BASE:
        return base

    if ctx.is_nonsense:
        return FunctionalConsequence.NONSENSE

    if ctx.reference_aa == "*" and ctx.altered_aa != "*":
        return FunctionalConsequence.STOP_LOSS

    if ctx.is_synonymous:
        return FunctionalConsequence.SYNONYMOUS

    return FunctionalConsequence.MISSENSE
