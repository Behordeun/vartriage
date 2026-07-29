"""ClinGen actionability curations database.

Parses a pre-processed TSV mapping gene symbols to their actionability
classification (whether established medical interventions exist).

Expected TSV columns: gene_symbol, intervention_type
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional

from vartriage._internal.path_safety import resolve_path

logger = logging.getLogger(__name__)


class ActionabilityDB:
    """ClinGen actionability lookup for medically actionable genes.

    Parameters
    ----------
    tsv_path : Path
        Path to the pre-processed clingen_actionability.tsv file.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._index: dict[str, str] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        """Parse the TSV and build the gene->intervention index."""
        if not tsv_path.exists():
            logger.warning("ClinGen actionability file not found: %s", tsv_path)
            return

        tsv_path = resolve_path(tsv_path)
        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("gene_symbol", "").strip()
                intervention = row.get("intervention_type", "").strip()

                if not gene:
                    continue

                # A gene is actionable if it has any intervention type
                self._index[gene] = intervention if intervention else "unspecified"

        logger.info(
            "ClinGen actionability loaded: %d actionable genes",
            len(self._index),
        )

    def is_actionable(self, gene_symbol: str) -> bool:
        """Check whether a gene has established medical interventions."""
        return gene_symbol in self._index

    def get_intervention_type(self, gene_symbol: str) -> Optional[str]:
        """Return the intervention type for an actionable gene.

        Returns
        -------
        str | None
            Intervention category (surveillance, therapeutic, etc.)
            or None if the gene is not actionable.
        """
        return self._index.get(gene_symbol)

    @property
    def gene_count(self) -> int:
        """Number of actionable genes loaded."""
        return len(self._index)
