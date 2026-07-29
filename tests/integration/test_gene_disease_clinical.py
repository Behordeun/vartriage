"""Clinical scenario integration test: affected pediatric epilepsy patient.

Simulates a real diagnostic workflow for a child presenting with:
- Febrile seizures progressing to afebrile (HP:0001250)
- Intellectual disability (HP:0001249)
- Developmental regression (HP:0002197)
- Speech delay (HP:0001263)

The synthetic VCF contains:
- A pathogenic SCN1A nonsense variant (Dravet syndrome candidate)
- A BRCA1 missense VUS (incidental, no phenotype overlap)
- A TP53 frameshift (incidental, partial phenotype overlap via cancer)
- An intergenic variant (noise)

This validates that the gene-disease linkage feature correctly:
1. Boosts SCN1A variants to the top via phenotype matching
2. Attaches disease associations (Dravet syndrome, GEFS+)
3. Reports constraint metrics (SCN1A is highly constrained)
4. Flags actionable genes where relevant
5. Degrades gracefully for intergenic variants
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from vartriage.knowledge.annotator import GeneKnowledgeAnnotator
from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.registry import GeneKnowledgeRegistry
from vartriage.models.variant import (
    AnnotatedVariant,
    ClinVarAssertion,
    FunctionalConsequence,
    Variant,
)


# --- Synthetic data builders ---


def _write_epilepsy_knowledge_dir(base: Path) -> Path:
    """Build a knowledge directory simulating real clinical gene data."""
    d = base / "knowledge"
    d.mkdir()

    (d / "omim_gene_disease.tsv").write_text(
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
        "SCN1A\tDravet syndrome\t607208\tAD\n"
        "SCN1A\tGeneralized epilepsy with febrile seizures plus, type 2\t604403\tAD\n"
        "BRCA1\tBreast-ovarian cancer, familial, 1\t604370\tAD\n"
        "TP53\tLi-Fraumeni syndrome\t151623\tAD\n"
    )

    # SCN1A overlaps strongly with our epilepsy patient's phenotype
    (d / "hpo_gene_annotations.tsv").write_text(
        "gene_symbol\thpo_terms\n"
        "SCN1A\tHP:0001250;HP:0001249;HP:0002197;HP:0001263\n"
        "BRCA1\tHP:0003002;HP:0002894\n"
        "TP53\tHP:0002664;HP:0003002\n"
    )

    (d / "clingen_validity.tsv").write_text(
        "gene_symbol\tvalidity_level\n"
        "SCN1A\tDefinitive\n"
        "BRCA1\tDefinitive\n"
        "TP53\tDefinitive\n"
    )

    # SCN1A is one of the most constrained genes in the genome
    (d / "gnomad_constraint.tsv").write_text(
        "gene_symbol\tpli\tloeuf\tmis_z\n"
        "SCN1A\t1.00\t0.07\t5.44\n"
        "BRCA1\t0.00\t1.17\t0.07\n"
        "TP53\t0.96\t0.20\t3.79\n"
    )

    (d / "clingen_actionability.tsv").write_text(
        "gene_symbol\tintervention_type\n"
        "BRCA1\tsurveillance\n"
    )

    return d


def _build_patient_variants() -> list[AnnotatedVariant]:
    """Construct annotated variants simulating a clinical VCF.

    Scenario:
    - SCN1A nonsense at chr2:166_848_884 (classic Dravet hotspot region)
    - BRCA1 missense (incidental finding)
    - TP53 frameshift (incidental)
    - Intergenic variant (background noise)
    """
    return [
        # SCN1A nonsense — primary diagnostic finding
        AnnotatedVariant(
            variant=Variant(
                chrom="chr2",
                pos=166848884,
                id="rs121918794",
                ref="C",
                alt="T",
                qual=99.0,
                filter_status="PASS",
                info={"DP": 45},
            ),
            consequence=FunctionalConsequence.NONSENSE,
            allele_frequency=0.000004,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC,
            frequency_unknown=False,
            clinvar_unknown=False,
            gene_name="SCN1A",
        ),
        # BRCA1 missense — incidental, no phenotype match for epilepsy patient
        AnnotatedVariant(
            variant=Variant(
                chrom="chr17",
                pos=43094000,
                id=None,
                ref="A",
                alt="G",
                qual=85.0,
                filter_status="PASS",
                info={"DP": 30},
            ),
            consequence=FunctionalConsequence.MISSENSE,
            allele_frequency=0.0003,
            clinvar_assertion=ClinVarAssertion.VUS,
            frequency_unknown=False,
            clinvar_unknown=False,
            gene_name="BRCA1",
        ),
        # TP53 frameshift — incidental
        AnnotatedVariant(
            variant=Variant(
                chrom="chr17",
                pos=7578000,
                id=None,
                ref="AG",
                alt="A",
                qual=92.0,
                filter_status="PASS",
                info={"DP": 38},
            ),
            consequence=FunctionalConsequence.FRAMESHIFT,
            allele_frequency=None,
            clinvar_assertion=None,
            frequency_unknown=True,
            clinvar_unknown=True,
            gene_name="TP53",
        ),
        # Intergenic noise
        AnnotatedVariant(
            variant=Variant(
                chrom="chr5",
                pos=50000000,
                id=None,
                ref="T",
                alt="C",
                qual=40.0,
                filter_status="PASS",
                info={},
            ),
            consequence=FunctionalConsequence.INTERGENIC,
            allele_frequency=0.35,
            clinvar_assertion=None,
            frequency_unknown=False,
            clinvar_unknown=True,
            gene_name=None,
        ),
    ]


# --- Tests ---


@pytest.fixture
def clinical_setup(tmp_path: Path) -> dict[str, object]:
    """Set up the clinical scenario: affected patient + knowledge base."""
    knowledge_dir = _write_epilepsy_knowledge_dir(tmp_path)

    # Patient HPO: seizures, intellectual disability, regression, speech delay
    patient_hpo = frozenset({
        "HP:0001250",  # Seizures
        "HP:0001249",  # Intellectual disability
        "HP:0002197",  # Seizure onset in first year of life / developmental regression
        "HP:0001263",  # Global developmental delay / speech delay
    })

    config = KnowledgeBaseConfig(
        data_dir=knowledge_dir,
        hpo_terms=patient_hpo,
    )

    annotator = GeneKnowledgeAnnotator(config)
    variants = _build_patient_variants()

    return {
        "annotator": annotator,
        "variants": variants,
        "config": config,
        "knowledge_dir": knowledge_dir,
    }


class TestAffectedPatientPhenotypeMatching:
    """Validate phenotype-driven gene prioritization for a real clinical case."""

    def test_scn1a_gets_perfect_phenotype_match(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """SCN1A has all 4 patient HPO terms -> overlap score = 1.0."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        scn1a_result = results[0]

        ctx = scn1a_result.gene_context
        assert ctx is not None
        assert ctx.phenotype_match_score == pytest.approx(1.0)

    def test_brca1_gets_zero_phenotype_match(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """BRCA1's HPO terms (breast cancer) don't overlap with epilepsy."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        brca1_result = results[1]

        ctx = brca1_result.gene_context
        assert ctx is not None
        assert ctx.phenotype_match_score == pytest.approx(0.0)

    def test_tp53_gets_zero_phenotype_match(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """TP53's HPO terms (neoplasm) don't overlap with epilepsy phenotype."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        tp53_result = results[2]

        ctx = tp53_result.gene_context
        assert ctx is not None
        assert ctx.phenotype_match_score == pytest.approx(0.0)

    def test_intergenic_gets_neutral_context(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """Intergenic variants have no gene -> neutral empty context."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        intergenic_result = results[3]

        ctx = intergenic_result.gene_context
        assert ctx is not None
        assert ctx.phenotype_match_score == pytest.approx(0.0)
        assert ctx.disease_associations == ()
        assert ctx.constraint is None


class TestDiseaseAssociationEnrichment:
    """Verify disease-gene associations are correctly attached."""

    def test_scn1a_has_two_disease_associations(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """SCN1A causes both Dravet and GEFS+."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[0].gene_context
        assert ctx is not None

        disease_names = [a.disease_name for a in ctx.disease_associations]
        assert "Dravet syndrome" in disease_names
        assert "Generalized epilepsy with febrile seizures plus, type 2" in disease_names

    def test_scn1a_inheritance_is_autosomal_dominant(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """Both SCN1A conditions are AD (de novo in Dravet, inherited in GEFS+)."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[0].gene_context
        assert ctx is not None

        modes = {a.inheritance_mode for a in ctx.disease_associations}
        assert modes == {"AD"}

    def test_brca1_has_breast_cancer_association(
        self, clinical_setup: dict[str, object]
    ) -> None:
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[1].gene_context
        assert ctx is not None

        assert len(ctx.disease_associations) == 1
        assert "Breast-ovarian cancer" in ctx.disease_associations[0].disease_name


class TestGeneConstraintMetrics:
    """Verify constraint data is correctly propagated for clinical interpretation."""

    def test_scn1a_is_highly_constrained(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """SCN1A has pLI=1.0, LOEUF=0.07 — extremely intolerant to LoF."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[0].gene_context
        assert ctx is not None
        assert ctx.constraint is not None
        assert ctx.constraint.pli == pytest.approx(1.0)
        assert ctx.constraint.loeuf == pytest.approx(0.07)
        assert ctx.constraint.mis_z == pytest.approx(5.44)
        assert ctx.constraint.is_lof_intolerant is True
        assert ctx.constraint.is_missense_constrained is True

    def test_brca1_is_not_lof_intolerant(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """BRCA1 has pLI=0.0 — tolerant to LoF (because AR conditions exist)."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[1].gene_context
        assert ctx is not None
        assert ctx.constraint is not None
        assert ctx.constraint.is_lof_intolerant is False


class TestActionabilityAnnotation:
    """Verify actionability flags for medically actionable genes."""

    def test_brca1_is_actionable(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """BRCA1 has ClinGen actionability (surveillance recommended)."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[1].gene_context
        assert ctx is not None
        assert ctx.is_actionable is True

    def test_scn1a_is_not_actionable(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """SCN1A doesn't have a ClinGen actionability curation in our data."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[0].gene_context
        assert ctx is not None
        assert ctx.is_actionable is False


class TestClinGenValidity:
    """Verify ClinGen gene-disease validity levels are reported."""

    def test_all_known_genes_have_definitive_validity(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """SCN1A, BRCA1, TP53 all have Definitive ClinGen validity."""
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        for i in range(3):  # First 3 variants are in known genes
            ctx = results[i].gene_context
            assert ctx is not None
            assert ctx.clingen_validity == "Definitive"

    def test_intergenic_has_no_validity(
        self, clinical_setup: dict[str, object]
    ) -> None:
        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        results = list(annotator.annotate(iter(variants)))
        ctx = results[3].gene_context
        assert ctx is not None
        assert ctx.clingen_validity is None


class TestPhenotypeBoostLogic:
    """Validate that the boost factor correctly separates signal from noise."""

    def test_boost_factor_for_perfect_match_is_two(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """score * (1 + 1.0) = score * 2.0 for perfect phenotype match."""
        from vartriage.knowledge.registry import apply_phenotype_boost

        boosted = apply_phenotype_boost(10.0, 1.0)
        assert boosted == pytest.approx(20.0)

    def test_boost_factor_for_no_match_is_one(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """score * (1 + 0.0) = score * 1.0 when no HPO overlap."""
        from vartriage.knowledge.registry import apply_phenotype_boost

        boosted = apply_phenotype_boost(10.0, 0.0)
        assert boosted == pytest.approx(10.0)

    def test_scn1a_would_outrank_brca1_after_boost(
        self, clinical_setup: dict[str, object]
    ) -> None:
        """Even if BRCA1 had a higher base score, phenotype boost lifts SCN1A.

        This models the tier-isolation requirement: a VUS BRCA1 with high
        base score should not outrank a Pathogenic SCN1A after phenotype
        boost within the same classification tier.
        """
        from vartriage.knowledge.registry import apply_phenotype_boost

        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        registry = annotator.registry

        # Hypothetical: BRCA1 has higher base score
        brca1_base = 15.0
        scn1a_base = 12.0

        brca1_overlap = registry.phenotype_overlap("BRCA1")
        scn1a_overlap = registry.phenotype_overlap("SCN1A")

        brca1_boosted = apply_phenotype_boost(brca1_base, brca1_overlap)
        scn1a_boosted = apply_phenotype_boost(scn1a_base, scn1a_overlap)

        # BRCA1: 15 * (1+0) = 15
        # SCN1A: 12 * (1+1) = 24
        assert scn1a_boosted is not None
        assert brca1_boosted is not None
        assert scn1a_boosted > brca1_boosted


class TestJsonOutputWithAffectedPatient:
    """Validate JSON serialization for gene-disease linkage fields."""

    def test_json_output_includes_gene_context(
        self, clinical_setup: dict[str, object], tmp_path: Path
    ) -> None:
        """Full JSON output includes disease associations and constraint."""
        from vartriage.models.variant import ScoredVariant, ClassifiedVariant
        from vartriage.models.variant import ACMGClassification, EvidenceTag
        from vartriage.reporting.json_writer import write_json

        annotator: GeneKnowledgeAnnotator = clinical_setup["annotator"]  # type: ignore[assignment]
        variants: list[AnnotatedVariant] = clinical_setup["variants"]  # type: ignore[assignment]

        enriched = list(annotator.annotate(iter(variants)))

        # Wrap in ScoredVariant -> ClassifiedVariant for JSON writer
        classified = []
        for av in enriched:
            scored = ScoredVariant(
                annotated=av,
                prioritization_score=10.0,
                composite_rank=1,
            )
            cv = ClassifiedVariant(
                scored=scored,
                classification=ACMGClassification.VUS,
                evidence_tags=frozenset(),
            )
            classified.append(cv)

        output = tmp_path / "results.json"
        write_json(classified, output)

        with open(output, encoding="utf-8") as f:
            data = json.load(f)

        # SCN1A variant should have full gene context in JSON
        scn1a_record = data[0]
        assert "disease_associations" in scn1a_record
        assert len(scn1a_record["disease_associations"]) == 2
        assert scn1a_record["disease_associations"][0]["disease_name"] == "Dravet syndrome"
        assert scn1a_record["disease_associations"][0]["inheritance_mode"] == "AD"
        assert scn1a_record["clingen_validity"] == "Definitive"
        assert scn1a_record["gene_constraint"]["pli"] == pytest.approx(1.0)
        assert scn1a_record["gene_constraint"]["loeuf"] == pytest.approx(0.07)
        assert scn1a_record["is_actionable"] is False
        assert scn1a_record["phenotype_match_score"] == pytest.approx(1.0)

        # Intergenic variant should have neutral context
        intergenic_record = data[3]
        assert intergenic_record["disease_associations"] == []
        assert intergenic_record["gene_constraint"] is None
        assert intergenic_record["phenotype_match_score"] == pytest.approx(0.0)
