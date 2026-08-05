"""ClinGen gene-disease validity classification database.

Parses a pre-processed TSV mapping gene symbols to their ClinGen
validity level (Definitive, Strong, Moderate, Limited, Disputed, Refuted).

Expected TSV columns: gene_symbol, validity_level
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from vartriage._internal.path_safety import resolve_path

logger = logging.getLogger(__name__)

VALID_LEVELS = frozenset(
    {
        "Definitive",
        "Strong",
        "Moderate",
        "Limited",
        "Disputed",
        "Refuted",
        "No Known Disease Relationship",
    }
)


class ClinGenValidityDB:
    """ClinGen gene-disease validity lookup.

    Parameters
    ----------
    tsv_path : Path
        Path to the pre-processed clingen_validity.tsv file.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._index: dict[str, str] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        """Parse the TSV and build the gene->validity index."""
        if not tsv_path.exists():
            logger.warning("ClinGen validity file not found: %s", tsv_path)
            return

        tsv_path = resolve_path(tsv_path)
        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("gene_symbol", "").strip()
                level = row.get("validity_level", "").strip()

                if not gene or not level:
                    continue

                if level not in VALID_LEVELS:
                    logger.debug(
                        "Skipping unrecognized validity level '%s' for gene %s",
                        level,
                        gene,
                    )
                    continue

                # A gene may appear multiple times if curated for different
                # conditions. Keep the strongest (most definitive) level.
                existing = self._index.get(gene)
                if existing is None or _level_rank(level) < _level_rank(existing):
                    self._index[gene] = level

        logger.info(
            "ClinGen validity loaded: %d genes curated",
            len(self._index),
        )

    def lookup(self, gene_symbol: str) -> str | None:
        """Return the ClinGen validity level for a gene.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (case-sensitive).

        Returns
        -------
        str | None
            Validity level string or None if gene is not curated.
        """
        return self._index.get(gene_symbol)

    @property
    def gene_count(self) -> int:
        """Number of genes with ClinGen validity curations."""
        return len(self._index)


# Lower rank = stronger evidence
_LEVEL_RANKING = {
    "Definitive": 0,
    "Strong": 1,
    "Moderate": 2,
    "Limited": 3,
    "Disputed": 4,
    "Refuted": 5,
    "No Known Disease Relationship": 6,
}


def _level_rank(level: str) -> int:
    """Numeric rank for ordering validity levels (lower = stronger)."""
    return _LEVEL_RANKING.get(level, 99)
