"""Unit tests for HPODatabase parsing and lookup."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.hpo import HPODatabase


@pytest.fixture
def hpo_tsv(tmp_path: Path) -> Path:
    """Write a minimal HPO annotations TSV."""
    content = (
        "gene_symbol\thpo_terms\n"
        "SCN1A\tHP:0001250;HP:0001249;HP:0002197\n"
        "MECP2\tHP:0001249;HP:0001250\n"
        "EMPTY_GENE\t\n"
    )
    tsv = tmp_path / "hpo_gene_annotations.tsv"
    tsv.write_text(content)
    return tsv


def test_loads_terms_as_frozenset(hpo_tsv: Path) -> None:
    db = HPODatabase(hpo_tsv)
    terms = db.get_terms("SCN1A")
    assert terms == frozenset({"HP:0001250", "HP:0001249", "HP:0002197"})


def test_shared_terms_across_genes(hpo_tsv: Path) -> None:
    db = HPODatabase(hpo_tsv)
    scn1a = db.get_terms("SCN1A")
    mecp2 = db.get_terms("MECP2")
    overlap = scn1a & mecp2
    assert overlap == frozenset({"HP:0001250", "HP:0001249"})


def test_unknown_gene_returns_empty_frozenset(hpo_tsv: Path) -> None:
    db = HPODatabase(hpo_tsv)
    assert db.get_terms("NOPE") == frozenset()


def test_empty_terms_row_skipped(hpo_tsv: Path) -> None:
    db = HPODatabase(hpo_tsv)
    assert db.get_terms("EMPTY_GENE") == frozenset()
    assert db.gene_count == 2


def test_missing_file_loads_empty(tmp_path: Path) -> None:
    db = HPODatabase(tmp_path / "missing.tsv")
    assert db.gene_count == 0
