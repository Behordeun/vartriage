"""Gene knowledge annotator stage for the pipeline.

Sits between the AnnotationEngine and PrioritizationEngine. Enriches
each AnnotatedVariant with gene-level context (disease associations,
constraint, actionability, phenotype match score). Optionally filters
by inheritance mode and actionability.

Also provides score boosting via boost_scores() for the prioritization
stage.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import replace

from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.registry import GeneKnowledgeRegistry, apply_phenotype_boost
from vartriage.models.variant import AnnotatedVariant, ScoredVariant

logger = logging.getLogger(__name__)


class GeneKnowledgeAnnotator:
    """Enrich AnnotatedVariants with gene-disease linkage context.

    Processes a stream of AnnotatedVariants and attaches a GeneContext
    to each. Filters by inheritance mode and actionability when configured.
    Provides phenotype score boosting for the prioritization stage.

    Parameters
    ----------
    config : KnowledgeBaseConfig
        Knowledge base configuration with data dir, HPO terms,
        inheritance mode filter, and actionable flag.
    """

    def __init__(self, config: KnowledgeBaseConfig) -> None:
        self._registry = GeneKnowledgeRegistry(config)
        self._inheritance_mode: str | None = config.inheritance_mode
        self._flag_actionable: bool = config.flag_actionable

    @property
    def registry(self) -> GeneKnowledgeRegistry:
        """Access the underlying gene knowledge registry."""
        return self._registry

    def annotate(
        self, variants: Iterator[AnnotatedVariant]
    ) -> Iterator[AnnotatedVariant]:
        """Attach gene context to each variant in the stream.

        Filtering order:
        1. Inheritance mode check (cheap OMIM lookup, avoids full context build)
        2. Build gene context (registry lookup + phenotype overlap)
        3. Actionability filter (when --flag-actionable is set)

        Intergenic variants (gene_name=None) pass through all filters.

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

            # Inheritance mode filter runs first to skip context build
            if self._inheritance_mode is not None and gene_symbol is not None:
                gene_modes = self._registry.omim.get_inheritance_modes(gene_symbol)
                if gene_modes and self._inheritance_mode not in gene_modes:
                    continue

            gene_context = self._registry.build_gene_context(gene_symbol)

            # Actionability filter: when --flag-actionable is set, only
            # yield variants in actionable genes (intergenic passes through)
            if (
                self._flag_actionable
                and gene_symbol is not None
                and not gene_context.is_actionable
            ):
                continue

            yield replace(variant, gene_context=gene_context)

    def boost_scores(self, scored: Iterator[ScoredVariant]) -> Iterator[ScoredVariant]:
        """Apply phenotype-based boost to prioritization scores.

        Multiplies prioritization_score by (1 + phenotype_match_score).
        Variants without gene_context or with zero overlap pass through
        unchanged.

        Parameters
        ----------
        scored : Iterator[ScoredVariant]
            Scored variants from the prioritization engine.

        Yields
        ------
        ScoredVariant
            Variants with boosted prioritization_score where applicable.
        """
        for variant in scored:
            ctx = variant.annotated.gene_context
            if ctx is None or ctx.phenotype_match_score == 0.0:
                yield variant
                continue

            boosted = apply_phenotype_boost(
                variant.prioritization_score, ctx.phenotype_match_score
            )
            yield replace(variant, prioritization_score=boosted)
