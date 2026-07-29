"""Unit tests for ConstraintDB parsing and lookup."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.constraint import ConstraintDB


@pytest.fixture
def constraint_tsv(tmp_path: Path) -> Path:
    """Write a minimal gnomAD constraint TSV."""
    content = (
        "gene_symbol\tpli\tloeuf\tmis_z\n"
        "TP53\t0.96\t0.20\t3.79\n"
        "BRCA1\t0.00\t1.17\t0.07\n"
        "DOT_GENE\t.\t.\t.\n"
    )
    tsv = tmp_path / "gnomad_constraint.tsv"
    tsv.write_text(content)
    return tsv


def test_parses_numeric_constraint(constraint_tsv: Path) -> None:
    db = ConstraintDB(constraint_tsv)
    c = db.lookup("TP53")
    assert c is not None
    assert c.pli == pytest.approx(0.96)
    assert c.loeuf == pytest.approx(0.20)
    assert c.mis_z == pytest.approx(3.79)


def test_is_lof_intolerant_property(constraint_tsv: Path) -> None:
    db = ConstraintDB(constraint_tsv)
    tp53 = db.lookup("TP53")
    brca1 = db.lookup("BRCA1")
    assert tp53 is not None and tp53.is_lof_intolerant is True
    assert brca1 is not None and brca1.is_lof_intolerant is False


def test_is_missense_constrained_property(constraint_tsv: Path) -> None:
    db = ConstraintDB(constraint_tsv)
    tp53 = db.lookup("TP53")
    brca1 = db.lookup("BRCA1")
    assert tp53 is not None and tp53.is_missense_constrained is True
    assert brca1 is not None and brca1.is_missense_constrained is False


def test_dot_values_skipped(constraint_tsv: Path) -> None:
    db = ConstraintDB(constraint_tsv)
    assert db.lookup("DOT_GENE") is None


def test_unknown_gene_returns_none(constraint_tsv: Path) -> None:
    db = ConstraintDB(constraint_tsv)
    assert db.lookup("NOPE") is None


def test_missing_file_loads_empty(tmp_path: Path) -> None:
    db = ConstraintDB(tmp_path / "missing.tsv")
    assert db.gene_count == 0
