"""Unit tests for OMIMDatabase parsing and lookup."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.omim import OMIMDatabase


@pytest.fixture
def omim_tsv(tmp_path: Path) -> Path:
    """Write a minimal OMIM TSV for testing."""
    content = (
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
        "BRCA1\tBreast-ovarian cancer, familial, 1\t604370\tAD\n"
        "BRCA1\tFanconi anemia, complementation group S\t617883\tAR\n"
        "CFTR\tCystic fibrosis\t219700\tAR\n"
    )
    tsv = tmp_path / "omim_gene_disease.tsv"
    tsv.write_text(content)
    return tsv


def test_loads_multiple_associations_per_gene(omim_tsv: Path) -> None:
    db = OMIMDatabase(omim_tsv)
    assocs = db.lookup("BRCA1")
    assert len(assocs) == 2
    assert assocs[0].disease_name == "Breast-ovarian cancer, familial, 1"
    assert assocs[0].mim_number == "604370"
    assert assocs[0].inheritance_mode == "AD"
    assert assocs[1].inheritance_mode == "AR"


def test_single_disease_gene(omim_tsv: Path) -> None:
    db = OMIMDatabase(omim_tsv)
    assocs = db.lookup("CFTR")
    assert len(assocs) == 1
    assert assocs[0].disease_name == "Cystic fibrosis"


def test_unknown_gene_returns_empty_tuple(omim_tsv: Path) -> None:
    db = OMIMDatabase(omim_tsv)
    assert db.lookup("NONEXISTENT") == ()


def test_get_inheritance_modes(omim_tsv: Path) -> None:
    db = OMIMDatabase(omim_tsv)
    modes = db.get_inheritance_modes("BRCA1")
    assert modes == frozenset({"AD", "AR"})


def test_missing_file_loads_empty(tmp_path: Path) -> None:
    db = OMIMDatabase(tmp_path / "does_not_exist.tsv")
    assert db.gene_count == 0
    assert db.lookup("BRCA1") == ()


def test_skips_rows_with_empty_gene_symbol(tmp_path: Path) -> None:
    content = (
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
        "\tSome disease\t123456\tAD\n"
        "VALID\tReal disease\t654321\tAR\n"
    )
    tsv = tmp_path / "test.tsv"
    tsv.write_text(content)
    db = OMIMDatabase(tsv)
    assert db.gene_count == 1
    assert db.lookup("VALID")[0].disease_name == "Real disease"


def test_skips_rows_with_empty_disease_name(tmp_path: Path) -> None:
    content = (
        "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
        "GENE1\t\t123456\tAD\n"
        "GENE2\tActual disease\t654321\tAR\n"
    )
    tsv = tmp_path / "test.tsv"
    tsv.write_text(content)
    db = OMIMDatabase(tsv)
    assert db.gene_count == 1
