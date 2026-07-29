"""HPO gene-phenotype annotation database.

Parses a pre-processed TSV mapping gene symbols to their associated
HPO term IDs. Used for phenotype-driven prioritization when the
clinician provides patient HPO terms.

Expected TSV columns: gene_symbol, hpo_terms (semicolon-separated HP:NNNNNNN IDs)
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from vartriage._internal.path_safety import safe_read_path

logger = logging.getLogger(__name__)


class HPODatabase:
    """Gene-to-HPO term mapping for phenotype overlap scoring.

    Parameters
    ----------
    tsv_path : Path
        Path to the pre-processed hpo_gene_annotations.tsv file.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._index: dict[str, frozenset[str]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        """Parse the TSV and build the gene->HPO terms index."""
        if not tsv_path.exists():
            logger.warning("HPO annotations file not found: %s", tsv_path)
            return

        tsv_path = safe_read_path(tsv_path, "HPO data")
        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("gene_symbol", "").strip()
                if not gene:
                    continue

                raw_terms = row.get("hpo_terms", "").strip()
                if not raw_terms:
                    continue

                terms = frozenset(
                    t.strip() for t in raw_terms.split(";") if t.strip()
                )
                if terms:
                    self._index[gene] = terms

        logger.info(
            "HPO database loaded: %d genes with phenotype annotations",
            len(self._index),
        )

    def get_terms(self, gene_symbol: str) -> frozenset[str]:
        """Return HPO terms associated with a gene.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (case-sensitive).

        Returns
        -------
        frozenset[str]
            Set of HPO term IDs. Empty frozenset if gene is not annotated.
        """
        return self._index.get(gene_symbol, frozenset())

    @property
    def gene_count(self) -> int:
        """Number of genes with HPO annotations."""
        return len(self._index)
