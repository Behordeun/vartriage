"""Unit tests for GeneKnowledgeAnnotator pipeline integration."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.annotator import GeneKnowledgeAnnotator
from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.models.variant import AnnotatedVariant, FunctionalConsequence, Variant


def _make_variant(gene_name: str | None = "BRCA1") -> AnnotatedVariant:
    """Build a minimal AnnotatedVariant for testing."""
    raw = Variant(
        chrom="chr17",
        pos=43094000,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
    )
    return AnnotatedVariant(
        variant=raw,
        consequence=FunctionalConsequence.MISSENSE,
        gene_name=gene_name,
    )


@pytest.fixture
def knowledge_dir(tmp_path: Path) -> Path:
    """Minimal knowledge directory with BRCA1 data."""
    d = tmp_path / "knowledge"
    d.mkdir()

    (d / "omim_gene_disease.tsv").write_text(
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
        "BRCA1\tBreast cancer\t604370\tAD\n"
    )
    (d / "hpo_gene_annotations.tsv").write_text(
        "gene_symbol\thpo_terms\nBRCA1\tHP:0003002;HP:0002894;HP:0010619\n"
    )
    (d / "clingen_validity.tsv").write_text(
        "gene_symbol\tvalidity_level\nBRCA1\tDefinitive\n"
    )
    (d / "gnomad_constraint.tsv").write_text(
        "gene_symbol\tpli\tloeuf\tmis_z\nBRCA1\t0.00\t1.17\t0.07\n"
    )
    (d / "clingen_actionability.tsv").write_text(
        "gene_symbol\tintervention_type\nBRCA1\tsurveillance\n"
    )
    return d


def test_annotator_attaches_gene_context(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    annotator = GeneKnowledgeAnnotator(config)

    variant = _make_variant("BRCA1")
    results = list(annotator.annotate(iter([variant])))

    assert len(results) == 1
    ctx = results[0].gene_context
    assert ctx is not None
    assert ctx.is_actionable is True
    assert ctx.clingen_validity == "Definitive"
    assert len(ctx.disease_associations) == 1


def test_annotator_with_phenotype_matching(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(
        data_dir=knowledge_dir,
        hpo_terms=frozenset({"HP:0003002", "HP:0010619"}),
    )
    annotator = GeneKnowledgeAnnotator(config)

    variant = _make_variant("BRCA1")
    results = list(annotator.annotate(iter([variant])))

    ctx = results[0].gene_context
    assert ctx is not None
    # 2 of 2 patient terms match BRCA1's 3 HPO terms -> 2/2 = 1.0
    assert ctx.phenotype_match_score == pytest.approx(1.0)


def test_annotator_unknown_gene_gets_neutral_context(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    annotator = GeneKnowledgeAnnotator(config)

    variant = _make_variant("NOVEL_GENE")
    results = list(annotator.annotate(iter([variant])))

    ctx = results[0].gene_context
    assert ctx is not None
    assert ctx.disease_associations == ()
    assert ctx.is_actionable is False
    assert ctx.phenotype_match_score == pytest.approx(0.0)


def test_annotator_none_gene_gets_neutral_context(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    annotator = GeneKnowledgeAnnotator(config)

    variant = _make_variant(None)
    results = list(annotator.annotate(iter([variant])))

    ctx = results[0].gene_context
    assert ctx is not None
    assert ctx.disease_associations == ()


def test_annotator_processes_stream_in_order(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    annotator = GeneKnowledgeAnnotator(config)

    v1 = _make_variant("BRCA1")
    v2 = _make_variant("NOVEL_GENE")
    v3 = _make_variant("BRCA1")

    results = list(annotator.annotate(iter([v1, v2, v3])))
    assert len(results) == 3
    # First and third share gene context through flyweight
    assert results[0].gene_context is not None
    assert results[2].gene_context is not None
    assert results[0].gene_context.is_actionable is True
    assert results[1].gene_context is not None
    assert results[1].gene_context.is_actionable is False
