"""Central gene knowledge registry composing all knowledge databases.

Loads all data sources once at pipeline start, provides O(1) lookups
per gene symbol, and maintains a flyweight cache so identical gene
contexts are never allocated twice. Phenotype overlap scoring is an
internal concern of this registry.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from vartriage.knowledge.actionability import ActionabilityDB
from vartriage.knowledge.clingen_validity import ClinGenValidityDB
from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.constraint import ConstraintDB
from vartriage.knowledge.hpo import HPODatabase
from vartriage.knowledge.models import (
    EMPTY_GENE_ANNOTATION,
    GeneAnnotation,
    GeneContext,
)
from vartriage.knowledge.omim import OMIMDatabase

logger = logging.getLogger(__name__)


def apply_phenotype_boost(
    base_score: Optional[float], overlap: float
) -> Optional[float]:
    """Apply phenotype boost: score * (1 + overlap).

    The boost factor ranges from 1.0 (no overlap) to 2.0 (perfect
    overlap). None base scores pass through unchanged.
    """
    if base_score is None:
        return None
    return base_score * (1.0 + overlap)


class GeneKnowledgeRegistry:
    """Central lookup for all gene-level annotations.

    Loaded once at pipeline start. All lookups are O(1) dict access.
    Maintains a flyweight cache to reuse GeneAnnotation instances
    across variants hitting the same gene. Owns phenotype overlap
    scoring internally.

    Parameters
    ----------
    config : KnowledgeBaseConfig
        Data directory and patient HPO terms.
    """

    def __init__(self, config: KnowledgeBaseConfig) -> None:
        data_dir = config.resolved_data_dir

        self._omim = OMIMDatabase(data_dir / "omim_gene_disease.tsv")
        self._hpo = HPODatabase(data_dir / "hpo_gene_annotations.tsv")
        self._clingen = ClinGenValidityDB(data_dir / "clingen_validity.tsv")
        self._constraint = ConstraintDB(data_dir / "gnomad_constraint.tsv")
        self._actionability = ActionabilityDB(data_dir / "clingen_actionability.tsv")

        self._patient_hpo_terms = config.hpo_terms

        # Flyweight cache: gene_symbol -> GeneAnnotation
        self._cache: dict[str, GeneAnnotation] = {}

        logger.info(
            "GeneKnowledgeRegistry initialized: data_dir=%s, patient_hpo_terms=%d",
            data_dir,
            len(self._patient_hpo_terms),
        )

    @property
    def omim(self) -> OMIMDatabase:
        """Direct access to the OMIM database (for inheritance mode queries)."""
        return self._omim

    @property
    def hpo(self) -> HPODatabase:
        """Direct access to the HPO database."""
        return self._hpo

    @property
    def patient_hpo_terms(self) -> frozenset[str]:
        """Patient HPO terms configured for this run."""
        return self._patient_hpo_terms

    @property
    def phenotype_boost_active(self) -> bool:
        """True if patient HPO terms were provided (boosting is enabled)."""
        return len(self._patient_hpo_terms) > 0

    def phenotype_overlap(self, gene_symbol: Optional[str]) -> float:
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
        if not self._patient_hpo_terms or gene_symbol is None:
            return 0.0

        gene_terms = self._hpo.get_terms(gene_symbol)
        if not gene_terms:
            return 0.0

        overlap = self._patient_hpo_terms & gene_terms
        return len(overlap) / len(self._patient_hpo_terms)

    def annotate_gene(self, gene_symbol: Optional[str]) -> GeneAnnotation:
        """Return all gene-level annotations for a gene symbol.

        Uses a flyweight cache so repeated lookups for the same gene
        reuse the same immutable GeneAnnotation object.

        Parameters
        ----------
        gene_symbol : str | None
            HGNC gene symbol. Returns EMPTY_GENE_ANNOTATION for None
            or unknown genes that have no data in any source.

        Returns
        -------
        GeneAnnotation
            Aggregated annotations from OMIM, ClinGen, gnomAD, HPO.
        """
        if gene_symbol is None:
            return EMPTY_GENE_ANNOTATION

        cached = self._cache.get(gene_symbol)
        if cached is not None:
            return cached

        disease_associations = self._omim.lookup(gene_symbol)
        clingen_validity = self._clingen.lookup(gene_symbol)
        constraint = self._constraint.lookup(gene_symbol)
        is_actionable = self._actionability.is_actionable(gene_symbol)
        actionability_type = self._actionability.get_intervention_type(gene_symbol)
        hpo_terms = self._hpo.get_terms(gene_symbol)

        has_any_data = (
            disease_associations
            or clingen_validity is not None
            or constraint is not None
            or is_actionable
            or hpo_terms
        )

        if not has_any_data:
            self._cache[gene_symbol] = EMPTY_GENE_ANNOTATION
            return EMPTY_GENE_ANNOTATION

        annotation = GeneAnnotation(
            disease_associations=disease_associations,
            clingen_validity=clingen_validity,
            constraint=constraint,
            is_actionable=is_actionable,
            actionability_type=actionability_type,
            hpo_terms=hpo_terms,
        )

        self._cache[gene_symbol] = annotation
        return annotation

    def build_gene_context(self, gene_symbol: Optional[str]) -> GeneContext:
        """Build a GeneContext for attachment to a variant.

        Computes phenotype overlap internally and combines it with
        the cached gene annotation.

        Parameters
        ----------
        gene_symbol : str | None
            Gene symbol from consequence annotation.

        Returns
        -------
        GeneContext
            Variant-facing gene context data including phenotype score.
        """
        annotation = self.annotate_gene(gene_symbol)
        overlap = self.phenotype_overlap(gene_symbol)

        return GeneContext(
            disease_associations=annotation.disease_associations,
            clingen_validity=annotation.clingen_validity,
            constraint=annotation.constraint,
            is_actionable=annotation.is_actionable,
            phenotype_match_score=overlap,
        )

    @property
    def cached_gene_count(self) -> int:
        """Number of genes currently in the flyweight cache."""
        return len(self._cache)
