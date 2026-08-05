"""Unit tests for GeneKnowledgeRegistry composition and flyweight cache."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.models import EMPTY_GENE_ANNOTATION
from vartriage.knowledge.registry import GeneKnowledgeRegistry


@pytest.fixture
def knowledge_dir(tmp_path: Path) -> Path:
    """Create a minimal knowledge data directory."""
    d = tmp_path / "knowledge"
    d.mkdir()

    (d / "omim_gene_disease.tsv").write_text(
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
        "BRCA1\tBreast cancer\t604370\tAD\n"
    )
    (d / "hpo_gene_annotations.tsv").write_text(
        "gene_symbol\thpo_terms\nBRCA1\tHP:0003002;HP:0002894\n"
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


def test_annotate_known_gene(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    registry = GeneKnowledgeRegistry(config)

    ann = registry.annotate_gene("BRCA1")
    assert len(ann.disease_associations) == 1
    assert ann.disease_associations[0].disease_name == "Breast cancer"
    assert ann.clingen_validity == "Definitive"
    assert ann.constraint is not None
    assert ann.constraint.pli == pytest.approx(0.00)
    assert ann.is_actionable is True
    assert ann.hpo_terms == frozenset({"HP:0003002", "HP:0002894"})


def test_annotate_unknown_gene_returns_empty(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    registry = GeneKnowledgeRegistry(config)

    ann = registry.annotate_gene("TOTALLY_NOVEL")
    assert ann is EMPTY_GENE_ANNOTATION


def test_annotate_none_returns_empty(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    registry = GeneKnowledgeRegistry(config)

    ann = registry.annotate_gene(None)
    assert ann is EMPTY_GENE_ANNOTATION


def test_flyweight_cache_reuses_instances(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    registry = GeneKnowledgeRegistry(config)

    first = registry.annotate_gene("BRCA1")
    second = registry.annotate_gene("BRCA1")
    assert first is second
    assert registry.cached_gene_count == 1


def test_build_gene_context_with_phenotype_score(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(
        data_dir=knowledge_dir,
        hpo_terms=frozenset({"HP:0003002"}),
    )
    registry = GeneKnowledgeRegistry(config)

    ctx = registry.build_gene_context("BRCA1")
    # 1 of 1 patient HPO terms matches BRCA1's annotations -> 1.0
    assert ctx.phenotype_match_score == pytest.approx(1.0)
    assert ctx.is_actionable is True
    assert ctx.clingen_validity == "Definitive"


def test_build_gene_context_unknown_gene(knowledge_dir: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=knowledge_dir)
    registry = GeneKnowledgeRegistry(config)

    ctx = registry.build_gene_context("UNKNOWN")
    assert ctx.disease_associations == ()
    assert ctx.is_actionable is False
    assert ctx.phenotype_match_score == pytest.approx(0.0)
