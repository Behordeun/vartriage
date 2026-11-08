"""Unit tests for ACMG benign criteria BP1 and BP3.

- BP1: Missense variant in a gene where loss of function is the established
  disease mechanism and missense variants are not known to cause disease.
- BP3: In-frame deletion or insertion in a repetitive sequence region without
  a known critical functional domain.
"""

from __future__ import annotations

from vartriage.annotation.clinvar_protein_index import (
    ClinVarProteinIndex,
    PathogenicMissense,
)
from vartriage.classification.acmg import ACMGClassifier
from vartriage.knowledge.models import GeneConstraint, GeneContext
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)


def _make_variant(
    *,
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    gene_name: str | None = "TESTGENE",
    info: dict[str, object] | None = None,
    allele_frequency: float | None = None,
    clinvar_assertion: ClinVarAssertion | None = None,
    mis_z: float = -0.5,
    revel_score: float | None = None,
    cadd_phred: float | None = None,
) -> ScoredVariant:
    raw = Variant(
        chrom="chr1",
        pos=1000,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
        info=info or {},
    )
    constraint = GeneConstraint(pli=0.1, loeuf=0.8, mis_z=mis_z)
    context = GeneContext(
        disease_associations=(),
        constraint=constraint,
    )
    annotated = AnnotatedVariant(
        variant=raw,
        consequence=consequence,
        gene_name=gene_name,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
        gene_context=context,
    )
    return ScoredVariant(
        annotated=annotated,
        revel_score=revel_score,
        cadd_phred=cadd_phred,
    )


class TestBP1:
    """BP1: Missense in a gene where LoF is the established mechanism."""

    def test_fires_for_missense_in_established_lof_gene(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
            mis_z=-0.8,
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 in tags

    def test_missing_protein_index_records_missing_source(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        classifier = ACMGClassifier(protein_index=None, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
            mis_z=-0.8,
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags
        assert "ClinVar_protein_index" in missing

    def test_does_not_fire_when_gene_not_in_lof_list(self) -> None:
        lof_genes = frozenset({"OTHERGENE"})
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags

    def test_missing_lof_list_records_missing_source(self) -> None:
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=None)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags
        assert "lof_mechanism_genes" in missing

    def test_does_not_fire_for_non_missense(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.SYNONYMOUS,
            gene_name="TESTGENE",
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags

    def test_does_not_fire_when_protein_index_has_pathogenic_missense(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        pm = PathogenicMissense(
            gene="TESTGENE",
            position=50,
            ref_aa="R",
            alt_aa="W",
            chrom="chr1",
            genomic_pos=1000,
            ref_allele="C",
            alt_allele="T",
        )
        pindex = ClinVarProteinIndex.from_variants([pm])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags

    def test_does_not_fire_when_gene_is_missense_constrained(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
            mis_z=3.5,  # Constrained
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags

    def test_does_not_fire_when_conflicting_pathogenic_tag_present(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
            mis_z=-0.5,
        )
        tags: set[EvidenceTag] = {EvidenceTag.PP3}
        missing: set[str] = set()

        classifier._evaluate_bp1(sv, tags, missing)

        assert EvidenceTag.BP1 not in tags


class TestBP3:
    """BP3: In-frame indels in repetitive regions without known function."""

    def test_fires_for_inframe_deletion_in_repeat_region(self) -> None:
        classifier = ACMGClassifier()
        sv = _make_variant(
            consequence=FunctionalConsequence.IN_FRAME_DELETION,
            info={"REPEAT": True},
            mis_z=-0.2,
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm4(sv, tags, missing)
        classifier._evaluate_bp3(sv, tags, missing)

        assert EvidenceTag.BP3 in tags
        assert EvidenceTag.PM4 not in tags

    def test_fires_for_inframe_insertion_in_repeat_region(self) -> None:
        classifier = ACMGClassifier()
        sv = _make_variant(
            consequence=FunctionalConsequence.IN_FRAME_INSERTION,
            info={"in_repeat": "1"},
            mis_z=-0.2,
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm4(sv, tags, missing)
        classifier._evaluate_bp3(sv, tags, missing)

        assert EvidenceTag.BP3 in tags
        assert EvidenceTag.PM4 not in tags

    def test_falsy_repeat_values_do_not_fire_bp3(self) -> None:
        classifier = ACMGClassifier()
        for falsy in ("0", False, "false", "no", 0):
            sv = _make_variant(
                consequence=FunctionalConsequence.IN_FRAME_DELETION,
                info={"REPEAT": falsy},
            )
            tags: set[EvidenceTag] = set()
            missing: set[str] = set()

            classifier._evaluate_pm4(sv, tags, missing)
            classifier._evaluate_bp3(sv, tags, missing)

            assert EvidenceTag.BP3 not in tags
            assert EvidenceTag.PM4 in tags

    def test_truthy_repeat_values_fire_bp3(self) -> None:
        classifier = ACMGClassifier()
        for key, truthy in [
            ("repeat", True),
            ("REPEAT", "1"),
            ("In_Repeat", "true"),
            ("str", 1),
            ("RMSK", "yes"),
        ]:
            sv = _make_variant(
                consequence=FunctionalConsequence.IN_FRAME_DELETION,
                info={key: truthy},
            )
            tags: set[EvidenceTag] = set()
            missing: set[str] = set()

            classifier._evaluate_pm4(sv, tags, missing)
            classifier._evaluate_bp3(sv, tags, missing)

            assert EvidenceTag.BP3 in tags
            assert EvidenceTag.PM4 not in tags

    def test_non_repeat_inframe_deletion_fires_pm4_not_bp3(self) -> None:
        classifier = ACMGClassifier()
        sv = _make_variant(
            consequence=FunctionalConsequence.IN_FRAME_DELETION,
            info={},  # Non-repeat
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm4(sv, tags, missing)
        classifier._evaluate_bp3(sv, tags, missing)

        assert EvidenceTag.PM4 in tags
        assert EvidenceTag.BP3 not in tags

    def test_does_not_fire_for_frameshift_in_repeat(self) -> None:
        classifier = ACMGClassifier()
        sv = _make_variant(
            consequence=FunctionalConsequence.FRAMESHIFT,
            info={"REPEAT": True},
        )
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_bp3(sv, tags, missing)

        assert EvidenceTag.BP3 not in tags


class TestBenignCombiningWithBP1BP3:
    """Combining rules integration with BP1 and BP3."""

    def test_bp1_plus_bp4_yields_likely_benign(self) -> None:
        lof_genes = frozenset({"TESTGENE"})
        pindex = ClinVarProteinIndex.from_variants([])
        classifier = ACMGClassifier(protein_index=pindex, lof_gene_list=lof_genes)
        # Low REVEL (0.25) triggers BP4; missense in LoF gene triggers BP1
        sv = _make_variant(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTGENE",
            mis_z=-0.5,
            revel_score=0.25,
        )
        results = list(classifier.classify(iter([sv])))
        res = results[0]

        assert EvidenceTag.BP1 in res.evidence_tags
        assert EvidenceTag.BP4 in res.evidence_tags
        assert res.classification == ACMGClassification.LIKELY_BENIGN

    def test_bp3_alone_yields_likely_benign(self) -> None:
        classifier = ACMGClassifier()
        sv = _make_variant(
            consequence=FunctionalConsequence.IN_FRAME_DELETION,
            info={"REPEAT": True},
        )
        results = list(classifier.classify(iter([sv])))
        res = results[0]

        assert EvidenceTag.BP3 in res.evidence_tags
        assert res.classification == ACMGClassification.LIKELY_BENIGN

    def test_bs1_plus_bp3_yields_likely_benign(self) -> None:
        classifier = ACMGClassifier()
        # BS1 (AF > 1%) contributes -4 points; BP3 contributes -1 point
        sv = _make_variant(
            consequence=FunctionalConsequence.IN_FRAME_DELETION,
            info={"REPEAT": True},
            allele_frequency=0.02,
        )
        results = list(classifier.classify(iter([sv])))
        res = results[0]

        assert EvidenceTag.BS1 in res.evidence_tags
        assert EvidenceTag.BP3 in res.evidence_tags
        assert res.classification == ACMGClassification.LIKELY_BENIGN
