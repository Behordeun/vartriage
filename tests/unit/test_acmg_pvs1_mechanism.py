"""PVS1 strength is conservative when the LoF mechanism is unknown.

PVS1 at Very Strong is valid only when loss of function is an established
disease mechanism for the gene. With no gene-list membership and no
constraint data, the mechanism is unknown, so a truncating variant gets
Strong, not Very Strong. Established LoF intolerance still earns Very Strong.
"""

from __future__ import annotations

from vartriage.classification.acmg import ACMGClassifier
from vartriage.knowledge.models import GeneConstraint, GeneContext
from vartriage.models.variant import (
    AnnotatedVariant,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)


def _nonsense(
    *,
    gene_name: str | None,
    constraint: GeneConstraint | None,
) -> ScoredVariant:
    variant = Variant(
        chrom="chr1",
        pos=2000,
        id=None,
        ref="C",
        alt="T",
        qual=None,
        filter_status="PASS",
        info={},
    )
    gene_context = (
        None
        if constraint is None and gene_name is None
        else GeneContext(disease_associations=(), constraint=constraint)
    )
    annotated = AnnotatedVariant(
        variant=variant,
        consequence=FunctionalConsequence.NONSENSE,
        gene_name=gene_name,
        gene_context=gene_context,
    )
    return ScoredVariant(annotated=annotated)


def _pvs1_tag(variant: ScoredVariant, lof_gene_list=None) -> set[EvidenceTag]:
    clf = ACMGClassifier(lof_gene_list=lof_gene_list)
    tags: set[EvidenceTag] = set()
    missing: set[str] = set()
    clf._evaluate_pvs1(variant, tags, missing)
    return tags


class TestPVS1UnknownMechanism:
    def test_no_constraint_no_gene_list_is_strong_not_very_strong(self) -> None:
        variant = _nonsense(gene_name="NOVELGENE", constraint=None)
        tags = _pvs1_tag(variant)
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags

    def test_lof_intolerant_constraint_earns_very_strong(self) -> None:
        constraint = GeneConstraint(pli=0.99, loeuf=0.1, mis_z=1.0)
        variant = _nonsense(gene_name="LOFGENE", constraint=constraint)
        tags = _pvs1_tag(variant)
        assert EvidenceTag.PVS1 in tags

    def test_lof_tolerant_constraint_is_strong(self) -> None:
        constraint = GeneConstraint(pli=0.01, loeuf=1.5, mis_z=0.0)
        variant = _nonsense(gene_name="TOLGENE", constraint=constraint)
        tags = _pvs1_tag(variant)
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags

    def test_gene_on_lof_list_earns_very_strong(self) -> None:
        variant = _nonsense(gene_name="TP53", constraint=None)
        tags = _pvs1_tag(variant, lof_gene_list={"TP53"})
        assert EvidenceTag.PVS1 in tags

    def test_gene_not_on_lof_list_is_strong(self) -> None:
        variant = _nonsense(gene_name="OTHER", constraint=None)
        tags = _pvs1_tag(variant, lof_gene_list={"TP53"})
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags
