"""Unit tests for PhenotypeRanker overlap computation and score boosting."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.hpo import HPODatabase
from vartriage.knowledge.phenotype import PhenotypeRanker


@pytest.fixture
def hpo_db(tmp_path: Path) -> HPODatabase:
    """Create a minimal HPO database."""
    content = (
        "gene_symbol\thpo_terms\n"
        "SCN1A\tHP:0001250;HP:0001249;HP:0002197;HP:0001263\n"
        "CFTR\tHP:0002205;HP:0006528\n"
    )
    tsv = tmp_path / "hpo_gene_annotations.tsv"
    tsv.write_text(content)
    return HPODatabase(tsv)


class TestComputeOverlap:
    def test_full_overlap_returns_one(self, hpo_db: HPODatabase) -> None:
        """All patient terms match the gene -> 1.0."""
        patient = frozenset({"HP:0001250", "HP:0001249"})
        ranker = PhenotypeRanker(patient, hpo_db)
        assert ranker.compute_overlap("SCN1A") == pytest.approx(1.0)

    def test_partial_overlap(self, hpo_db: HPODatabase) -> None:
        """Some patient terms match -> fraction."""
        patient = frozenset({"HP:0001250", "HP:0001249", "HP:0099999"})
        ranker = PhenotypeRanker(patient, hpo_db)
        # 2 out of 3 patient terms match SCN1A
        assert ranker.compute_overlap("SCN1A") == pytest.approx(2.0 / 3.0)

    def test_no_overlap_returns_zero(self, hpo_db: HPODatabase) -> None:
        patient = frozenset({"HP:0099999"})
        ranker = PhenotypeRanker(patient, hpo_db)
        assert ranker.compute_overlap("SCN1A") == pytest.approx(0.0)

    def test_unknown_gene_returns_zero(self, hpo_db: HPODatabase) -> None:
        patient = frozenset({"HP:0001250"})
        ranker = PhenotypeRanker(patient, hpo_db)
        assert ranker.compute_overlap("NONEXISTENT") == pytest.approx(0.0)

    def test_none_gene_returns_zero(self, hpo_db: HPODatabase) -> None:
        patient = frozenset({"HP:0001250"})
        ranker = PhenotypeRanker(patient, hpo_db)
        assert ranker.compute_overlap(None) == pytest.approx(0.0)

    def test_empty_patient_terms_returns_zero(self, hpo_db: HPODatabase) -> None:
        """No patient terms -> phenotype boosting is disabled."""
        ranker = PhenotypeRanker(frozenset(), hpo_db)
        assert ranker.compute_overlap("SCN1A") == pytest.approx(0.0)
        assert ranker.is_active is False


class TestBoostScore:
    def test_boost_with_full_overlap(self, hpo_db: HPODatabase) -> None:
        """Perfect overlap -> score * 2.0 (maximum boost)."""
        patient = frozenset({"HP:0001250", "HP:0001249"})
        ranker = PhenotypeRanker(patient, hpo_db)
        boosted = ranker.boost_score(10.0, 1.0)
        assert boosted == pytest.approx(20.0)

    def test_boost_with_no_overlap(self, hpo_db: HPODatabase) -> None:
        """Zero overlap -> score * 1.0 (no change)."""
        patient = frozenset({"HP:0001250"})
        ranker = PhenotypeRanker(patient, hpo_db)
        boosted = ranker.boost_score(10.0, 0.0)
        assert boosted == pytest.approx(10.0)

    def test_boost_with_partial_overlap(self, hpo_db: HPODatabase) -> None:
        """50% overlap -> score * 1.5."""
        patient = frozenset({"HP:0001250"})
        ranker = PhenotypeRanker(patient, hpo_db)
        boosted = ranker.boost_score(10.0, 0.5)
        assert boosted == pytest.approx(15.0)

    def test_boost_none_score_returns_none(self, hpo_db: HPODatabase) -> None:
        """None base score passes through untouched."""
        patient = frozenset({"HP:0001250"})
        ranker = PhenotypeRanker(patient, hpo_db)
        assert ranker.boost_score(None, 1.0) is None
