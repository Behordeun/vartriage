"""Phenotype-driven gene prioritization via HPO term overlap.

Computes how well a gene's associated phenotypes match the patient's
presenting symptoms (provided as HPO terms). The overlap score drives
a multiplicative boost on variant prioritization scores.
"""

from __future__ import annotations

import logging
from typing import Optional

from vartriage.knowledge.hpo import HPODatabase

logger = logging.getLogger(__name__)


class PhenotypeRanker:
    """Boost variant prioritization based on HPO phenotype overlap.

    Computes exact term overlap between patient HPO terms and each
    gene's annotated phenotype terms. Future versions may incorporate
    semantic similarity via the HPO ontology graph.

    Parameters
    ----------
    patient_hpo_terms : frozenset[str]
        Patient's HPO term IDs (e.g., {"HP:0001250", "HP:0001249"}).
    hpo_db : HPODatabase
        Loaded gene-to-HPO annotation database.
    """

    def __init__(
        self, patient_hpo_terms: frozenset[str], hpo_db: HPODatabase
    ) -> None:
        self._patient_terms = patient_hpo_terms
        self._hpo_db = hpo_db

    @property
    def is_active(self) -> bool:
        """True if patient HPO terms were provided (phenotype boosting enabled)."""
        return len(self._patient_terms) > 0

    def compute_overlap(self, gene_symbol: Optional[str]) -> float:
        """Fraction of patient HPO terms matching this gene's phenotype.

        Returns 0.0 when no patient terms are configured, gene is None,
        or the gene has no HPO annotations.

        Parameters
        ----------
        gene_symbol : str | None
            HGNC gene symbol.

        Returns
        -------
        float
            Overlap fraction in range [0.0, 1.0].
        """
        if not self._patient_terms or gene_symbol is None:
            return 0.0

        gene_terms = self._hpo_db.get_terms(gene_symbol)
        if not gene_terms:
            return 0.0

        overlap = self._patient_terms & gene_terms
        return len(overlap) / len(self._patient_terms)

    def boost_score(
        self, base_score: Optional[float], overlap: float
    ) -> Optional[float]:
        """Apply phenotype boost: score * (1 + overlap).

        The boost factor ranges from 1.0 (no overlap) to 2.0 (perfect
        overlap). This keeps phenotype boosting bounded and prevents
        runaway inflation of prioritization scores.

        Parameters
        ----------
        base_score : float | None
            Original prioritization score. None passes through unchanged.
        overlap : float
            Phenotype overlap fraction from compute_overlap().

        Returns
        -------
        float | None
            Boosted score, or None if base_score was None.
        """
        if base_score is None:
            return None
        return base_score * (1.0 + overlap)
