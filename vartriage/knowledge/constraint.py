"""gnomAD gene constraint metrics database.

Parses a pre-processed TSV containing per-gene constraint scores
(pLI, LOEUF, mis_z) from gnomAD.

Expected TSV columns: gene_symbol, pli, loeuf, mis_z
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional

from vartriage._internal.path_safety import safe_read_path
from vartriage.knowledge.models import GeneConstraint

logger = logging.getLogger(__name__)


class ConstraintDB:
    """gnomAD gene constraint metric lookup.

    Parameters
    ----------
    tsv_path : Path
        Path to the pre-processed gnomad_constraint.tsv file.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._index: dict[str, GeneConstraint] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        """Parse the TSV and build the gene->constraint index."""
        if not tsv_path.exists():
            logger.warning("gnomAD constraint file not found: %s", tsv_path)
            return

        tsv_path = safe_read_path(tsv_path, "Constraint data")
        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("gene_symbol", "").strip()
                if not gene:
                    continue

                pli_raw = row.get("pli", "").strip()
                loeuf_raw = row.get("loeuf", "").strip()
                mis_z_raw = row.get("mis_z", "").strip()

                # Skip rows with missing numeric data (dot or empty)
                if not pli_raw or pli_raw == "." or \
                   not loeuf_raw or loeuf_raw == "." or \
                   not mis_z_raw or mis_z_raw == ".":
                    continue

                try:
                    constraint = GeneConstraint(
                        pli=float(pli_raw),
                        loeuf=float(loeuf_raw),
                        mis_z=float(mis_z_raw),
                    )
                except (ValueError, TypeError):
                    logger.debug(
                        "Skipping unparseable constraint row for gene %s", gene
                    )
                    continue

                self._index[gene] = constraint

        logger.info(
            "gnomAD constraint loaded: %d genes with metrics",
            len(self._index),
        )

    def lookup(self, gene_symbol: str) -> Optional[GeneConstraint]:
        """Return constraint metrics for a gene.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (case-sensitive).

        Returns
        -------
        GeneConstraint | None
            Constraint metrics or None if gene is not in gnomAD.
        """
        return self._index.get(gene_symbol)

    @property
    def gene_count(self) -> int:
        """Number of genes with constraint data."""
        return len(self._index)
