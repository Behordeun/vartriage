"""PVS1 Very Strong is gated on an established LoF disease mechanism, not population constraint.

ClinGen PVS1 (Abou Tayoun et al. 2018) assigns Very Strong to a null variant only when
loss of function is an established disease mechanism for the gene, supplied here as a
curated gene list. gnomAD pLI/LOEUF measure intolerance to HETEROZYGOUS loss of function
in the population, which misclassifies recessive and viable-carrier disease genes
(BRCA1/2, CFTR, the mismatch-repair genes, PAH, GAA), so constraint must never grant
Very Strong on its own. A gene whose mechanism is not established fails closed to Strong.
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


def _nonsense_in_gene(
    gene_name: str,
    pli: float = 0.0,
    loeuf: float = 1.0,
    mis_z: float = 0.0,
    with_constraint: bool = True,
) -> ScoredVariant:
    v = Variant(
        chrom="chr13",
        pos=32340000,
        id=None,
        ref="C",
        alt="T",
        qual=30.0,
        filter_status="PASS",
    )
    gene_context = None
    if with_constraint:
        constraint = GeneConstraint(pli=pli, loeuf=loeuf, mis_z=mis_z)
        gene_context = GeneContext(disease_associations=(), constraint=constraint)
    annotated = AnnotatedVariant(
        variant=v,
        consequence=FunctionalConsequence.NONSENSE,
        allele_frequency=0.00001,
        gene_name=gene_name,
        gene_context=gene_context,
    )
    return ScoredVariant(annotated=annotated, revel_score=None)


class TestPVS1LofMechanismGate:
    def test_listed_recessive_gene_with_low_pli_reaches_very_strong(self) -> None:
        # BRCA2-like: definitive LoF mechanism, pLI ~ 0 because heterozygous carriers
        # are viable. Constraint must not deny Very Strong for a listed gene.
        sv = _nonsense_in_gene("BRCA2", pli=0.0, loeuf=0.80)
        classifier = ACMGClassifier(lof_gene_list=frozenset({"BRCA2"}))
        tags = next(iter(classifier.classify(iter([sv])))).evidence_tags
        assert EvidenceTag.PVS1 in tags
        assert EvidenceTag.PVS1_STRONG not in tags

    def test_listed_gene_very_strong_even_without_constraint_data(self) -> None:
        sv = _nonsense_in_gene("CFTR", with_constraint=False)
        classifier = ACMGClassifier(lof_gene_list=frozenset({"CFTR"}))
        tags = next(iter(classifier.classify(iter([sv])))).evidence_tags
        assert EvidenceTag.PVS1 in tags
        assert EvidenceTag.PVS1_STRONG not in tags

    def test_unlisted_gene_fails_closed_to_strong(self) -> None:
        sv = _nonsense_in_gene("UNKNOWNGENE", pli=0.99, loeuf=0.10)
        classifier = ACMGClassifier(lof_gene_list=frozenset({"BRCA2"}))
        tags = next(iter(classifier.classify(iter([sv])))).evidence_tags
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags

    def test_high_pli_does_not_grant_very_strong_when_list_present(self) -> None:
        sv = _nonsense_in_gene("CONSTRAINEDBUTUNLISTED", pli=0.995, loeuf=0.05)
        classifier = ACMGClassifier(lof_gene_list=frozenset({"BRCA2"}))
        tags = next(iter(classifier.classify(iter([sv])))).evidence_tags
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags

    def test_pli_fallback_does_not_grant_very_strong_without_a_list(self) -> None:
        # The core fix: with no curated list supplied, a high-pLI gene must still
        # NOT reach Very Strong. Population constraint is not a mechanism signal.
        sv = _nonsense_in_gene("HIGHPLI_NO_LIST", pli=0.995, loeuf=0.05)
        classifier = ACMGClassifier()
        tags = next(iter(classifier.classify(iter([sv])))).evidence_tags
        assert EvidenceTag.PVS1_STRONG in tags
        assert EvidenceTag.PVS1 not in tags
