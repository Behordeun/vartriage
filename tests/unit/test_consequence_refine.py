"""A coding SNV's consequence is refined from its resolved codon.

A single-base substitution inside a CDS is not always missense: it may be
nonsense (introduces a stop), synonymous (same amino acid), or stop-loss
(removes an existing stop). The refiner turns a base coding call plus a
resolved codon context into the correct amino-acid-level consequence, so
the answer does not depend on which interval backend produced the base
call.
"""

from __future__ import annotations

from dataclasses import dataclass

from vartriage.annotation.consequence_refine import (
    refine_coding_snv_consequence,
)
from vartriage.models.variant import FunctionalConsequence


@dataclass(frozen=True)
class _Ctx:
    reference_aa: str
    altered_aa: str
    is_synonymous: bool
    is_nonsense: bool


def _ctx(ref_aa: str, alt_aa: str) -> _Ctx:
    return _Ctx(
        reference_aa=ref_aa,
        altered_aa=alt_aa,
        is_synonymous=(ref_aa == alt_aa),
        is_nonsense=(alt_aa == "*" and ref_aa != "*"),
    )


class TestRefineCodingSnvConsequence:
    def test_stop_gain_becomes_nonsense(self) -> None:
        result = refine_coding_snv_consequence(
            FunctionalConsequence.MISSENSE, _ctx("Q", "*")
        )
        assert result == FunctionalConsequence.NONSENSE

    def test_same_amino_acid_becomes_synonymous(self) -> None:
        result = refine_coding_snv_consequence(
            FunctionalConsequence.MISSENSE, _ctx("L", "L")
        )
        assert result == FunctionalConsequence.SYNONYMOUS

    def test_amino_acid_change_stays_missense(self) -> None:
        result = refine_coding_snv_consequence(
            FunctionalConsequence.MISSENSE, _ctx("A", "V")
        )
        assert result == FunctionalConsequence.MISSENSE

    def test_stop_loss_becomes_stop_loss(self) -> None:
        result = refine_coding_snv_consequence(
            FunctionalConsequence.MISSENSE, _ctx("*", "R")
        )
        assert result == FunctionalConsequence.STOP_LOSS

    def test_no_codon_context_leaves_base_unchanged(self) -> None:
        result = refine_coding_snv_consequence(FunctionalConsequence.MISSENSE, None)
        assert result == FunctionalConsequence.MISSENSE

    def test_non_coding_base_is_never_upgraded(self) -> None:
        # A splice-site base call must not be reinterpreted even if a codon
        # context is somehow present.
        result = refine_coding_snv_consequence(
            FunctionalConsequence.SPLICE_SITE, _ctx("Q", "*")
        )
        assert result == FunctionalConsequence.SPLICE_SITE

    def test_intergenic_base_is_never_upgraded(self) -> None:
        result = refine_coding_snv_consequence(
            FunctionalConsequence.INTERGENIC, _ctx("Q", "*")
        )
        assert result == FunctionalConsequence.INTERGENIC
