"""ACMG Secondary Findings (SF v3.2) gene filter.

Flags variants in medically actionable genes regardless of the primary
gene panel filter. When enabled, variants in ACMG SF genes bypass
gene-list filtering and appear in a dedicated "Secondary Findings"
section of the clinical report.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from vartriage.models.variant import AnnotatedVariant

logger = logging.getLogger(__name__)

_DEFAULT_SF_PATH = Path(__file__).parent.parent / "data" / "acmg_sf_v3.2.txt"


class SecondaryFindingsFilter:
    """Identify variants in ACMG Secondary Findings genes.

    Loads the SF gene list (one symbol per line, # comments skipped)
    and provides membership testing and stream splitting.

    Parameters
    ----------
    gene_list_path
        Path to the SF gene list file. Defaults to the shipped
        ACMG SF v3.2 list.
    """

    def __init__(self, gene_list_path: Optional[Path] = None) -> None:
        path = gene_list_path or _DEFAULT_SF_PATH
        self._genes = self._load_genes(path)
        logger.info("SecondaryFindingsFilter loaded %d SF genes", len(self._genes))

    @property
    def gene_count(self) -> int:
        """Number of genes in the SF list."""
        return len(self._genes)

    def is_secondary_finding(self, gene_name: Optional[str]) -> bool:
        """Check if a gene is in the ACMG SF list.

        Parameters
        ----------
        gene_name
            Gene symbol to check. None returns False.

        Returns
        -------
        bool
            True if the gene is in the secondary findings list.
        """
        if gene_name is None:
            return False
        normalized = gene_name.strip().upper()
        if not normalized:
            return False
        return normalized in self._genes

    def split_stream(
        self, variants: list[AnnotatedVariant]
    ) -> tuple[list[AnnotatedVariant], list[AnnotatedVariant]]:
        """Split a materialized variant list into primary and secondary findings.

        The clinical report pipeline materializes variants before this
        point (for sorting by tier). Primary-panel filtering happens
        earlier in the pipeline via GeneFilter; all variants reaching
        this method have already passed panel filtering. Secondary
        findings are variants that additionally match an SF gene.

        Parameters
        ----------
        variants
            Materialized variant list (post-panel-filter).

        Returns
        -------
        tuple[list, list]
            (all_variants, secondary_findings). The first list contains
            all input variants unchanged. The second contains the subset
            matching SF genes.
        """
        primary: list[AnnotatedVariant] = []
        secondary: list[AnnotatedVariant] = []

        for variant in variants:
            primary.append(variant)
            if self.is_secondary_finding(variant.gene_name):
                secondary.append(variant)

        return primary, secondary

    def _load_genes(self, path: Path) -> frozenset[str]:
        """Load gene symbols from file, uppercased for case-insensitive matching.

        Raises
        ------
        FileNotFoundError
            If the gene list file does not exist. A missing SF list is
            a configuration error that could cause clinically relevant
            variants to be silently excluded.
        """
        if not path.is_file():
            raise FileNotFoundError(
                f"ACMG Secondary Findings gene list not found or is not a regular file: {path}. "
                f"This file is required when --secondary-findings is enabled."
            )

        genes: set[str] = set()
        with open(path, encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                genes.add(stripped.upper())

        return frozenset(genes)
