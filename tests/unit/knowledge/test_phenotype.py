"""Unit tests for phenotype overlap computation and score boosting."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.registry import GeneKnowledgeRegistry, apply_phenotype_boost


@pytest.fixture
def knowledge_dir(tmp_path: Path) -> Path:
    """Create a minimal knowledge directory with HPO data."""
    d = tmp_path / "knowledge"
    d.mkdir()
    (d / "hpo_gene_annotations.tsv").write_text(
        "gene_symbol\thpo_terms\n"
        "SCN1A\tHP:0001250;HP:0001249;HP:0002197;HP:0001263\n"
        "CFTR\tHP:0002205;HP:0006528\n"
    )
    (d / "omim_gene_disease.tsv").write_text(
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
    )
    (d / "clingen_validity.tsv").write_text("gene_symbol\tvalidity_level\n")
    (d / "gnomad_constraint.tsv").write_text("gene_symbol\tpli\tloeuf\tmis_z\n")
    (d / "clingen_actionability.tsv").write_text("gene_symbol\tintervention_type\n")
    return d


class TestPhenotypeOverlap:
    def test_full_overlap_returns_one(self, knowledge_dir: Path) -> None:
        """All patient terms match the gene -> 1.0."""
        config = KnowledgeBaseConfig(
            data_dir=knowledge_dir,
            hpo_terms=frozenset({"HP:0001250", "HP:0001249"}),
        )
        registry = GeneKnowledgeRegistry(config)
        assert registry.phenotype_overlap("SCN1A") == pytest.approx(1.0)

    def test_partial_overlap(self, knowledge_dir: Path) -> None:
        """Some patient terms match -> fraction."""
        config = KnowledgeBaseConfig(
            data_dir=knowledge_dir,
            hpo_terms=frozenset({"HP:0001250", "HP:0001249", "HP:0099999"}),
        )
        registry = GeneKnowledgeRegistry(config)
        assert registry.phenotype_overlap("SCN1A") == pytest.approx(2.0 / 3.0)

    def test_no_overlap_returns_zero(self, knowledge_dir: Path) -> None:
        config = KnowledgeBaseConfig(
            data_dir=knowledge_dir,
            hpo_terms=frozenset({"HP:0099999"}),
        )
        registry = GeneKnowledgeRegistry(config)
        assert registry.phenotype_overlap("SCN1A") == pytest.approx(0.0)

    def test_unknown_gene_returns_zero(self, knowledge_dir: Path) -> None:
        config = KnowledgeBaseConfig(
            data_dir=knowledge_dir,
            hpo_terms=frozenset({"HP:0001250"}),
        )
        registry = GeneKnowledgeRegistry(config)
        assert registry.phenotype_overlap("NONEXISTENT") == pytest.approx(0.0)

    def test_none_gene_returns_zero(self, knowledge_dir: Path) -> None:
        config = KnowledgeBaseConfig(
            data_dir=knowledge_dir,
            hpo_terms=frozenset({"HP:0001250"}),
        )
        registry = GeneKnowledgeRegistry(config)
        assert registry.phenotype_overlap(None) == pytest.approx(0.0)

    def test_empty_patient_terms_returns_zero(self, knowledge_dir: Path) -> None:
        """No patient terms -> phenotype boosting is disabled."""
        config = KnowledgeBaseConfig(data_dir=knowledge_dir)
        registry = GeneKnowledgeRegistry(config)
        assert registry.phenotype_overlap("SCN1A") == pytest.approx(0.0)
        assert registry.phenotype_boost_active is False


class TestApplyPhenotypeBoost:
    def test_boost_with_full_overlap(self) -> None:
        """Perfect overlap -> score * 2.0 (maximum boost)."""
        assert apply_phenotype_boost(10.0, 1.0) == pytest.approx(20.0)

    def test_boost_with_no_overlap(self) -> None:
        """Zero overlap -> score * 1.0 (no change)."""
        assert apply_phenotype_boost(10.0, 0.0) == pytest.approx(10.0)

    def test_boost_with_partial_overlap(self) -> None:
        """50% overlap -> score * 1.5."""
        assert apply_phenotype_boost(10.0, 0.5) == pytest.approx(15.0)

    def test_boost_none_score_returns_none(self) -> None:
        """None base score passes through untouched."""
        assert apply_phenotype_boost(None, 1.0) is None
