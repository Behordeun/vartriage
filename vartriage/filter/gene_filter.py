"""Gene-list-based variant filtering.

Restricts the annotated variant stream to only those variants whose gene
symbol appears in a user-supplied plain text file. Matching is
case-insensitive; gene symbols are stored uppercase-normalized in a
frozenset.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator

from vartriage.models.config import GeneFilterConfig
from vartriage.models.variant import AnnotatedVariant

logger = logging.getLogger(__name__)


class GeneFilter:
    """Filter annotated variants by gene symbol membership.

    Loads a gene list file at construction time, builds a frozenset of
    uppercase gene symbols. Only variants whose gene_name matches a member
    of this set will pass through when ``apply`` is called.

    Parameters
    ----------
    config : GeneFilterConfig
        Must contain the path to the gene list file.

    Raises
    ------
    FileNotFoundError
        If the gene list file does not exist at the specified path.
    ValueError
        If the gene list file contains zero valid gene symbols after
        blank and comment lines are removed.
    """

    def __init__(self, config: GeneFilterConfig) -> None:
        self._genes: frozenset[str] = self._load_gene_list(config.gene_list_path)
        self._matched_genes: set[str] = set()

    @property
    def genes(self) -> frozenset[str]:
        """The loaded gene set (uppercase-normalized)."""
        return self._genes

    @property
    def unmatched_genes(self) -> frozenset[str]:
        """Genes from the list that received zero matching variants."""
        return frozenset(self._genes - self._matched_genes)

    def apply(self, variants: Iterator[AnnotatedVariant]) -> Iterator[AnnotatedVariant]:
        """Yield variants whose gene_name is in the loaded gene set.

        After the iterator is fully consumed, logs a WARNING for any
        genes in the list that had zero matching variants.

        Parameters
        ----------
        variants : Iterator[AnnotatedVariant]
            Input annotated variant stream.

        Yields
        ------
        AnnotatedVariant
            Variants matching the gene list.
        """
        self._matched_genes = set()

        for variant in variants:
            gene = self._extract_gene(variant)
            if gene is not None and gene in self._genes:
                self._matched_genes.add(gene)
                yield variant

        self._log_unmatched()

    def _extract_gene(self, variant: AnnotatedVariant) -> str | None:
        """Extract and normalize gene symbol from a variant.

        Parameters
        ----------
        variant : AnnotatedVariant
            The variant to extract a gene name from.

        Returns
        -------
        str | None
            Uppercased gene name, or None for intergenic variants.
        """
        gene_name = variant.gene_name
        if gene_name is None:
            return None
        return gene_name.upper()

    def _log_unmatched(self) -> None:
        """Log WARNING for genes with zero matching variants."""
        unmatched = self._genes - self._matched_genes
        if unmatched:
            sorted_genes = sorted(unmatched)
            logger.warning(
                "Gene list genes with no matching variants: %s",
                ", ".join(sorted_genes),
            )

    @staticmethod
    def _load_gene_list(path: Path) -> frozenset[str]:
        """Parse a gene list file into a frozenset of uppercase symbols.

        Each non-blank, non-comment line is treated as a single gene
        symbol. Leading and trailing whitespace is stripped, and the
        symbol is uppercased before storage.

        Parameters
        ----------
        path : Path
            Path to the gene list file.

        Returns
        -------
        frozenset[str]
            Unique uppercase gene symbols found in the file.

        Raises
        ------
        FileNotFoundError
            File does not exist.
        ValueError
            File contains zero valid gene symbols.
        """
        if not path.exists():
            raise FileNotFoundError(f"Gene list file not found: {path}")

        genes: set[str] = set()
        with open(path, "r") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#"):
                    continue
                genes.add(stripped.upper())

        if not genes:
            raise ValueError(
                f"Gene list file contains no valid gene " f"symbols: {path}"
            )

        return frozenset(genes)
