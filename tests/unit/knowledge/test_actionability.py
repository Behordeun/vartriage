"""Unit tests for ActionabilityDB parsing and lookup."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.actionability import ActionabilityDB


@pytest.fixture
def actionability_tsv(tmp_path: Path) -> Path:
    """Write a minimal ClinGen actionability TSV."""
    content = (
        "gene_symbol\tintervention_type\n"
        "BRCA1\tsurveillance\n"
        "LDLR\ttherapeutic\n"
        "NO_TYPE\t\n"
    )
    tsv = tmp_path / "clingen_actionability.tsv"
    tsv.write_text(content)
    return tsv


def test_is_actionable_true(actionability_tsv: Path) -> None:
    db = ActionabilityDB(actionability_tsv)
    assert db.is_actionable("BRCA1") is True
    assert db.is_actionable("LDLR") is True


def test_is_actionable_false_for_unknown(actionability_tsv: Path) -> None:
    db = ActionabilityDB(actionability_tsv)
    assert db.is_actionable("RANDOM_GENE") is False


def test_get_intervention_type(actionability_tsv: Path) -> None:
    db = ActionabilityDB(actionability_tsv)
    assert db.get_intervention_type("BRCA1") == "surveillance"
    assert db.get_intervention_type("LDLR") == "therapeutic"


def test_empty_intervention_defaults_to_unspecified(actionability_tsv: Path) -> None:
    db = ActionabilityDB(actionability_tsv)
    assert db.is_actionable("NO_TYPE") is True
    assert db.get_intervention_type("NO_TYPE") == "unspecified"


def test_missing_file_loads_empty(tmp_path: Path) -> None:
    db = ActionabilityDB(tmp_path / "nope.tsv")
    assert db.gene_count == 0
    assert db.is_actionable("BRCA1") is False
