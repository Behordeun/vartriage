"""Unit tests for the SpliceAI SQLite backend loader."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from vartriage.models.config import PrioritizationConfig
from vartriage.prioritization.spliceai_db import SpliceAISQLiteLoader


@pytest.fixture
def spliceai_db(tmp_path: Path) -> Path:
    """Create a minimal SpliceAI SQLite test database."""
    db_path = tmp_path / "spliceai_test.sqlite"
    conn = sqlite3.connect(str(db_path))

    conn.execute(
        "CREATE TABLE chr17 ("
        "pos int, ref text, alt text, "
        "ds_ag real, ds_al real, ds_dg real, ds_dl real, "
        "dp_ag int, dp_al int, dp_dg int, dp_dl int)"
    )
    conn.execute(
        "CREATE TABLE chr22 ("
        "pos int, ref text, alt text, "
        "ds_ag real, ds_al real, ds_dg real, ds_dl real, "
        "dp_ag int, dp_al int, dp_dg int, dp_dl int)"
    )

    # chr17 test data: splice-disrupting variant
    conn.execute(
        "INSERT INTO chr17 VALUES (43091429, 'T', 'G', 0.02, 0.01, 0.92, 0.03, "
        "-5, 10, -2, 8)"
    )
    # chr17 test data: benign variant (low scores)
    conn.execute(
        "INSERT INTO chr17 VALUES (7674220, 'C', 'T', 0.01, 0.0, 0.02, 0.0, "
        "3, -1, 5, -2)"
    )
    # chr22 test data: moderate score
    conn.execute(
        "INSERT INTO chr22 VALUES (28695868, 'G', 'A', 0.55, 0.02, 0.03, 0.01, "
        "4, -3, 7, -1)"
    )

    conn.commit()
    conn.close()
    return db_path


class TestSpliceAISQLiteLoader:
    """Tests for SpliceAISQLiteLoader."""

    def test_lookup_returns_max_delta_score(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("chr17", 43091429, "T", "G")
        assert score == pytest.approx(0.92)

    def test_lookup_low_score_variant(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("chr17", 7674220, "C", "T")
        assert score == pytest.approx(0.02)

    def test_lookup_missing_variant_returns_none(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("chr17", 99999999, "A", "G")
        assert score is None

    def test_lookup_missing_table_returns_none(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("chr1", 100, "A", "G")
        assert score is None

    def test_chromosome_normalization_bare_number(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("17", 43091429, "T", "G")
        assert score == pytest.approx(0.92)

    def test_chromosome_normalization_uppercase(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("CHR17", 43091429, "T", "G")
        assert score == pytest.approx(0.92)

    def test_chromosome_normalization_mixed_case(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            score = loader.lookup("Chr22", 28695868, "G", "A")
        assert score == pytest.approx(0.55)

    def test_lookup_batch_returns_correct_order(self, spliceai_db: Path) -> None:
        variants = [
            ("chr22", 28695868, "G", "A"),
            ("chr17", 43091429, "T", "G"),
            ("chr17", 99999999, "A", "G"),
        ]
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            results = loader.lookup_batch(variants)

        assert len(results) == 3
        assert results[0] == pytest.approx(0.55)
        assert results[1] == pytest.approx(0.92)
        assert results[2] is None

    def test_lookup_batch_empty_list(self, spliceai_db: Path) -> None:
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            results = loader.lookup_batch([])
        assert results == []

    def test_lookup_batch_missing_table_skips_gracefully(
        self, spliceai_db: Path
    ) -> None:
        variants = [
            ("chrX", 100, "A", "G"),
            ("chr17", 43091429, "T", "G"),
        ]
        with SpliceAISQLiteLoader(spliceai_db) as loader:
            results = loader.lookup_batch(variants)

        assert results[0] is None
        assert results[1] == pytest.approx(0.92)

    def test_nonexistent_db_raises_valueerror(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="not found"):
            SpliceAISQLiteLoader(tmp_path / "nonexistent.sqlite")

    def test_context_manager_closes_connection(self, spliceai_db: Path) -> None:
        loader = SpliceAISQLiteLoader(spliceai_db)
        loader.close()
        # After close, operations should fail
        with pytest.raises(sqlite3.ProgrammingError):
            loader.lookup("chr17", 43091429, "T", "G")


class TestPrioritizationConfigMutualExclusion:
    """Tests for PrioritizationConfig spliceai backend mutual exclusion."""

    def test_both_backends_raises_valueerror(self) -> None:
        with pytest.raises(ValueError, match="Cannot configure both"):
            PrioritizationConfig(
                spliceai_scores_path=Path("scores.tsv"),
                spliceai_db_path=Path("scores.sqlite"),
            )

    def test_only_tsv_is_valid(self) -> None:
        config = PrioritizationConfig(spliceai_scores_path=Path("scores.tsv"))
        assert config.spliceai_scores_path == Path("scores.tsv")
        assert config.spliceai_db_path is None

    def test_only_sqlite_is_valid(self) -> None:
        config = PrioritizationConfig(spliceai_db_path=Path("scores.sqlite"))
        assert config.spliceai_db_path == Path("scores.sqlite")
        assert config.spliceai_scores_path is None

    def test_neither_backend_is_valid(self) -> None:
        config = PrioritizationConfig()
        assert config.spliceai_scores_path is None
        assert config.spliceai_db_path is None
