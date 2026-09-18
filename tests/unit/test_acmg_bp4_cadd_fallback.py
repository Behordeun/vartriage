"""BP4 has a CADD fallback for missense variants when REVEL is absent.

The pathogenic computational criterion (PP3) can fire from REVEL or, when
REVEL is absent, from SpliceAI. The benign criterion (BP4) had only REVEL
for missense, so a missense variant with no REVEL accrued no benign
computational evidence even when CADD was low. This restores the symmetry:
a low CADD supports BP4 for missense when REVEL is unavailable.
"""

from __future__ import annotations

from vartriage.classification.acmg import ACMGClassifier
from vartriage.models.variant import (
    AnnotatedVariant,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)


def _missense(
    *,
    revel_score: float | None,
    cadd_phred: float | None,
) -> ScoredVariant:
    variant = Variant(
        chrom="chr1",
        pos=3000,
        id=None,
        ref="A",
        alt="G",
        qual=None,
        filter_status="PASS",
        info={},
    )
    annotated = AnnotatedVariant(
        variant=variant,
        consequence=FunctionalConsequence.MISSENSE,
        gene_name="GENEX",
    )
    return ScoredVariant(
        annotated=annotated,
        revel_score=revel_score,
        cadd_phred=cadd_phred,
    )


def _bp4(variant: ScoredVariant) -> tuple[set[EvidenceTag], set[str]]:
    clf = ACMGClassifier()
    tags: set[EvidenceTag] = set()
    missing: set[str] = set()
    clf._evaluate_bp4(variant, tags, missing)
    return tags, missing


class TestBP4MissenseCaddFallback:
    def test_low_cadd_fires_bp4_when_revel_absent(self) -> None:
        tags, _ = _bp4(_missense(revel_score=None, cadd_phred=3.0))
        assert EvidenceTag.BP4 in tags

    def test_high_cadd_does_not_fire_bp4_when_revel_absent(self) -> None:
        tags, _ = _bp4(_missense(revel_score=None, cadd_phred=25.0))
        assert EvidenceTag.BP4 not in tags

    def test_revel_present_takes_precedence_over_cadd(self) -> None:
        # REVEL above the benign band: no BP4, even with a low CADD.
        tags, _ = _bp4(_missense(revel_score=0.8, cadd_phred=3.0))
        assert EvidenceTag.BP4 not in tags
        assert EvidenceTag.BP4_MODERATE not in tags

    def test_revel_benign_still_fires_bp4(self) -> None:
        tags, _ = _bp4(_missense(revel_score=0.05, cadd_phred=25.0))
        assert EvidenceTag.BP4_MODERATE in tags

    def test_no_revel_no_cadd_records_revel_missing(self) -> None:
        tags, missing = _bp4(_missense(revel_score=None, cadd_phred=None))
        assert EvidenceTag.BP4 not in tags
        assert "REVEL" in missing
