"""OMIM gene-disease association database.

Parses a pre-processed TSV mapping gene symbols to their associated
diseases, MIM numbers, and inheritance modes. Multiple diseases per
gene are supported (one row per association in the TSV).

Expected TSV columns: gene_symbol, disease_name, mim_number, inheritance_mode
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional

from vartriage.knowledge.models import DiseaseAssociation

logger = logging.getLogger(__name__)


class OMIMDatabase:
    """Gene-disease association lookup from OMIM data.

    Loads a TSV file at construction and builds an internal dict
    for O(1) gene symbol lookups.

    Parameters
    ----------
    tsv_path : Path
        Path to the pre-processed omim_gene_disease.tsv file.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._index: dict[str, tuple[DiseaseAssociation, ...]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        """Parse the TSV and build the gene->diseases index."""
        if not tsv_path.exists():
            logger.warning("OMIM data file not found: %s", tsv_path)
            return

        accumulator: dict[str, list[DiseaseAssociation]] = {}

        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("gene_symbol", "").strip()
                if not gene:
                    continue

                disease_name = row.get("disease_name", "").strip()
                mim_number = row.get("mim_number", "").strip()
                inheritance = row.get("inheritance_mode", "").strip()

                if not disease_name:
                    continue

                assoc = DiseaseAssociation(
                    disease_name=disease_name,
                    mim_number=mim_number,
                    inheritance_mode=inheritance,
                )

                if gene not in accumulator:
                    accumulator[gene] = []
                accumulator[gene].append(assoc)

        # Freeze lists into tuples for immutability
        self._index = {
            gene: tuple(assocs) for gene, assocs in accumulator.items()
        }

        logger.info(
            "OMIM database loaded: %d genes, %d associations",
            len(self._index),
            sum(len(v) for v in self._index.values()),
        )

    def lookup(self, gene_symbol: str) -> tuple[DiseaseAssociation, ...]:
        """Return disease associations for a gene symbol.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (case-sensitive).

        Returns
        -------
        tuple[DiseaseAssociation, ...]
            Associated diseases. Empty tuple if the gene is not in OMIM.
        """
        return self._index.get(gene_symbol, ())

    def get_inheritance_modes(self, gene_symbol: str) -> frozenset[str]:
        """Return all inheritance modes reported for a gene.

        Useful for filtering variants by expected zygosity.
        """
        assocs = self._index.get(gene_symbol, ())
        return frozenset(a.inheritance_mode for a in assocs if a.inheritance_mode)

    @property
    def gene_count(self) -> int:
        """Number of genes loaded from the OMIM data."""
        return len(self._index)
