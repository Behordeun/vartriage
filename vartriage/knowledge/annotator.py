"""Gene knowledge annotator stage for the pipeline.

Sits between the AnnotationEngine and PrioritizationEngine. Enriches
each AnnotatedVariant with gene-level context (disease associations,
constraint, actionability, phenotype match score).
"""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Iterator, Optional

from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.models import GeneContext
from vartriage.knowledge.phenotype import PhenotypeRanker
from vartriage.knowledge.registry import GeneKnowledgeRegistry
from vartriage.models.variant import AnnotatedVariant

logger = logging.getLogger(__name__)


class GeneKnowledgeAnnotator:
    """Enrich AnnotatedVariants with gene-disease linkage context.

    Processes a stream of AnnotatedVariants (which already have
    gene_name resolved by the AnnotationEngine) and attaches a
    GeneContext to each one.

    Parameters
    ----------
    config : KnowledgeBaseConfig
        Knowledge base configuration with data dir and HPO terms.
    """

    def __init__(self, config: KnowledgeBaseConfig) -> None:
        self._registry = GeneKnowledgeRegistry(config)
        self._ranker = PhenotypeRanker(
            patient_hpo_terms=config.hpo_terms,
            hpo_db=self._registry.hpo,
        )

    @property
    def registry(self) -> GeneKnowledgeRegistry:
        """Access the underlying gene knowledge registry."""
        return self._registry

    @property
    def phenotype_ranker(self) -> PhenotypeRanker:
        """Access the phenotype ranker for score boosting."""
        return self._ranker

    def annotate(
        self, variants: Iterator[AnnotatedVariant]
    ) -> Iterator[AnnotatedVariant]:
        """Attach gene context to each variant in the stream.

        Parameters
        ----------
        variants : Iterator[AnnotatedVariant]
            Input stream of annotated variants (gene_name already resolved).

        Yields
        ------
        AnnotatedVariant
            Variants with gene_context field populated.
        """
        for variant in variants:
            gene_symbol = variant.gene_name

            # Compute phenotype overlap for this gene
            overlap = self._ranker.compute_overlap(gene_symbol)

            # Build the context object (uses flyweight cache internally)
            gene_context = self._registry.build_gene_context(
                gene_symbol=gene_symbol,
                phenotype_match_score=overlap,
            )

            # Attach context via dataclass replace (AnnotatedVariant is frozen)
            yield replace(variant, gene_context=gene_context)
