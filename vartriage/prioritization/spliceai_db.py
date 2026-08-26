"""SpliceAI SQLite backend for precomputed delta score lookups.

Queries the OpenCRAVAT SpliceAI SQLite database directly, returning
the maximum delta score (across acceptor gain/loss, donor gain/loss)
for each variant. Eliminates the need to pre-filter scores into TSV
files for per-analysis use.

Database schema (one table per chromosome: chr1..chr22, chrX, chrY, chrM):
    pos int, ref text, alt text,
    ds_ag real, ds_al real, ds_dg real, ds_dl real,
    dp_ag int, dp_al int, dp_dg int, dp_dl int
"""

from __future__ import annotations

import logging
import sqlite3
from pathlib import Path

logger = logging.getLogger(__name__)


class SpliceAISQLiteLoader:
    """Query SpliceAI precomputed scores from an OpenCRAVAT SQLite database.

    Returns max(ds_ag, ds_al, ds_dg, ds_dl) per variant — the maximum
    predicted splice-altering probability across all four mechanisms.

    Parameters
    ----------
    db_path : Path
        Path to the OpenCRAVAT SpliceAI SQLite file.

    Raises
    ------
    ValueError
        If the database file does not exist.
    """

    def __init__(self, db_path: Path) -> None:
        if not db_path.exists():
            raise ValueError(f"SpliceAI database not found: {db_path}")
        uri = f"file:{db_path}?mode=ro"
        self._conn = sqlite3.connect(uri, uri=True)
        self._conn.row_factory = None
        self._available_tables: set[str] | None = None

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> float | None:
        """Look up the max SpliceAI delta score for a single variant.

        Parameters
        ----------
        chrom : str
            Chromosome identifier (accepts "chr22", "22", "Chr22").
        pos : int
            1-based genomic position.
        ref : str
            Reference allele.
        alt : str
            Alternate allele.

        Returns
        -------
        float or None
            Max delta score, or None if variant not in database.
        """
        table = self._normalize_chrom(chrom)
        if not self._table_exists(table):
            return None
        cursor = self._conn.execute(
            f"SELECT max(ds_ag, ds_al, ds_dg, ds_dl) FROM [{table}] "
            "WHERE pos = ? AND ref = ? AND alt = ?",
            (pos, ref, alt),
        )
        row = cursor.fetchone()
        if row is None or row[0] is None:
            return None
        return float(row[0])

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[float | None]:
        """Look up scores for a batch of variants, grouped by chromosome.

        Parameters
        ----------
        variants : list[tuple[str, int, str, str]]
            List of (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[float | None]
            Scores in the same order as input; None where not found.
        """
        results: list[float | None] = [None] * len(variants)

        by_chrom: dict[str, list[tuple[int, int, str, str]]] = {}
        for idx, (chrom, pos, ref, alt) in enumerate(variants):
            table = self._normalize_chrom(chrom)
            by_chrom.setdefault(table, []).append((idx, pos, ref, alt))

        for table, entries in by_chrom.items():
            if not self._table_exists(table):
                logger.warning(
                    "SpliceAI database has no table '%s', skipping %d variants",
                    table,
                    len(entries),
                )
                continue
            for idx, pos, ref, alt in entries:
                cursor = self._conn.execute(
                    f"SELECT max(ds_ag, ds_al, ds_dg, ds_dl) FROM [{table}] "
                    "WHERE pos = ? AND ref = ? AND alt = ?",
                    (pos, ref, alt),
                )
                row = cursor.fetchone()
                if row is not None and row[0] is not None:
                    results[idx] = float(row[0])

        return results

    def close(self) -> None:
        """Close the database connection."""
        self._conn.close()

    def __enter__(self) -> SpliceAISQLiteLoader:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def _table_exists(self, table: str) -> bool:
        """Check if a chromosome table exists in the database."""
        if self._available_tables is None:
            cursor = self._conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
            self._available_tables = {row[0] for row in cursor.fetchall()}
        return table in self._available_tables

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        """Normalize chromosome identifier to table name format.

        Accepts "chr22", "22", "Chr22", "CHR22" and returns "chr22".
        """
        c = chrom.lower().removeprefix("chr")
        return f"chr{c}"
