"""Knowledge loaders validate their TSV header instead of loading nothing.

A one-character header rename (pli -> pLI) used to load zero rows with only
an info log, so a downstream gene looked unconstrained or unassociated for a
reason no one could see. Each loader now checks its required columns up front
and raises a clear error naming the file and the missing columns.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.constraint import ConstraintDB
from vartriage.knowledge.omim import OMIMDatabase


def _write(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


class TestConstraintHeaderValidation:
    def test_renamed_column_raises_naming_the_missing_column(
        self, tmp_path: Path
    ) -> None:
        # "pLI" instead of "pli"
        tsv = _write(
            tmp_path,
            "constraint.tsv",
            "gene_symbol\tpLI\tloeuf\tmis_z\nBRCA1\t0.99\t0.1\t3.5\n",
        )
        with pytest.raises(ValueError, match="pli"):
            ConstraintDB(tsv)

    def test_correct_header_loads(self, tmp_path: Path) -> None:
        tsv = _write(
            tmp_path,
            "constraint.tsv",
            "gene_symbol\tpli\tloeuf\tmis_z\nBRCA1\t0.99\t0.1\t3.5\n",
        )
        db = ConstraintDB(tsv)
        assert db.gene_count == 1


class TestOMIMHeaderValidation:
    def test_missing_column_raises(self, tmp_path: Path) -> None:
        tsv = _write(
            tmp_path,
            "omim.tsv",
            "gene_symbol\tdisease_name\tmim_number\nBRCA1\tcancer\t113705\n",
        )
        with pytest.raises(ValueError, match="inheritance_mode"):
            OMIMDatabase(tsv)

    def test_correct_header_loads(self, tmp_path: Path) -> None:
        tsv = _write(
            tmp_path,
            "omim.tsv",
            "gene_symbol\tdisease_name\tmim_number\tinheritance_mode\n"
            "BRCA1\tcancer\t113705\tAD\n",
        )
        db = OMIMDatabase(tsv)
        assert db is not None
