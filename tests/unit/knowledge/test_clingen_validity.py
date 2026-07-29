"""Unit tests for ClinGenValidityDB parsing and lookup."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.clingen_validity import ClinGenValidityDB


@pytest.fixture
def validity_tsv(tmp_path: Path) -> Path:
    """Write a minimal ClinGen validity TSV."""
    content = (
        "gene_symbol\tvalidity_level\n"
        "BRCA1\tDefinitive\n"
        "LMNA\tStrong\n"
        "GENE_X\tLimited\n"
    )
    tsv = tmp_path / "clingen_validity.tsv"
    tsv.write_text(content)
    return tsv


def test_lookup_definitive(validity_tsv: Path) -> None:
    db = ClinGenValidityDB(validity_tsv)
    assert db.lookup("BRCA1") == "Definitive"


def test_lookup_strong(validity_tsv: Path) -> None:
    db = ClinGenValidityDB(validity_tsv)
    assert db.lookup("LMNA") == "Strong"


def test_unknown_gene_returns_none(validity_tsv: Path) -> None:
    db = ClinGenValidityDB(validity_tsv)
    assert db.lookup("UNKNOWN") is None


def test_keeps_strongest_when_duplicated(tmp_path: Path) -> None:
    """When a gene appears twice, the more definitive level wins."""
    content = (
        "gene_symbol\tvalidity_level\n"
        "MLH1\tModerate\n"
        "MLH1\tDefinitive\n"
    )
    tsv = tmp_path / "test.tsv"
    tsv.write_text(content)
    db = ClinGenValidityDB(tsv)
    assert db.lookup("MLH1") == "Definitive"


def test_skips_unrecognized_levels(tmp_path: Path) -> None:
    content = (
        "gene_symbol\tvalidity_level\n"
        "GENE1\tInventedLevel\n"
        "GENE2\tStrong\n"
    )
    tsv = tmp_path / "test.tsv"
    tsv.write_text(content)
    db = ClinGenValidityDB(tsv)
    assert db.lookup("GENE1") is None
    assert db.lookup("GENE2") == "Strong"


def test_missing_file_loads_empty(tmp_path: Path) -> None:
    db = ClinGenValidityDB(tmp_path / "nope.tsv")
    assert db.gene_count == 0
