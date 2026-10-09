"""ClinGen-curated loss-of-function disease-mechanism gene set.

Parses the pre-processed lof_mechanism_genes.tsv (produced by
scripts/derive_lof_mechanism_genes.py from the ClinGen gene curation list).
A gene is present when ClinGen curated it with a Haploinsufficiency Score of
3 (dominant haploinsufficiency) or 30 (autosomal recessive phenotype), either
of which establishes loss of function as a disease mechanism per the ClinGen
PVS1 decision tree (Abou Tayoun et al. 2018).

Expected TSV columns: gene_symbol, hi_score (lines beginning with # are comments).

The set gates PVS1 Very Strong: a null variant in a gene on this list may reach
Very Strong, a gene off the list fails closed to Strong. Population constraint
(pLI/LOEUF) is deliberately not used for this decision, because it measures
intolerance to heterozygous loss of function and misclassifies recessive and
viable-carrier disease genes.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from vartriage._internal.path_safety import resolve_path
from vartriage.knowledge._validation import require_columns

logger = logging.getLogger(__name__)


def default_lof_mechanism_path() -> Path:
    """Path to the bundled LoF-mechanism gene set shipped with the package."""
    return Path(__file__).resolve().parent.parent / "data" / "lof_mechanism_genes.tsv"


class LofMechanismGeneSet:
    """Lookup for genes with an established LoF disease mechanism.

    Parameters
    ----------
    tsv_path : Path
        Path to the pre-processed lof_mechanism_genes.tsv file.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._genes: frozenset[str] = frozenset()
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            logger.warning("LoF-mechanism gene set not found: %s", tsv_path)
            return

        tsv_path = resolve_path(tsv_path)
        genes: set[str] = set()
        with open(tsv_path, newline="", encoding="utf-8") as fh:
            rows = (line for line in fh if not line.startswith("#"))
            reader = csv.DictReader(rows, delimiter="\t")
            require_columns(reader, ["gene_symbol", "hi_score"], tsv_path)
            for row in reader:
                gene = row.get("gene_symbol", "").strip()
                if gene:
                    genes.add(gene)
        self._genes = frozenset(genes)
        logger.info("Loaded %d LoF-mechanism genes from %s", len(self._genes), tsv_path)

    def __contains__(self, gene_symbol: str) -> bool:
        return gene_symbol in self._genes

    @property
    def genes(self) -> frozenset[str]:
        """The curated LoF-mechanism gene symbols."""
        return self._genes

    @property
    def gene_count(self) -> int:
        return len(self._genes)
