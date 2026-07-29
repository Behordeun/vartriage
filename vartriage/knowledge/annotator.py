"""Gene knowledge annotator stage for the pipeline.

Sits between the AnnotationEngine and PrioritizationEngine. Enriches
each AnnotatedVariant with gene-level context (disease associations,
constraint, actionability, phenotype match score). Optionally filters
by inheritance mode.
"""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Iterator, Optional

from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.registry import GeneKnowledgeRegistry
from vartriage.models.variant import AnnotatedVariant

logger = logging.getLogger(__name__)


class GeneKnowledgeAnnotator:
    """Enrich AnnotatedVariants with gene-disease linkage context.

    Thin wrapper around GeneKnowledgeRegistry that processes a stream
    of AnnotatedVariants and attaches a GeneContext to each one.
    Optionally filters variants by gene inheritance mode.

    Parameters
    ----------
    config : KnowledgeBaseConfig
        Knowledge base configuration with data dir, HPO terms,
        inheritance mode filter, and actionable flag.
    """

    def __init__(self, config: KnowledgeBaseConfig) -> None:
        self._registry = GeneKnowledgeRegistry(config)
        self._inheritance_mode: Optional[str] = config.inheritance_mode
        self._flag_actionable: bool = config.flag_actionable

    @property
    def registry(self) -> GeneKnowledgeRegistry:
        """Access the underlying gene knowledge registry."""
        return self._registry

    def annotate(
        self, variants: Iterator[AnnotatedVariant]
    ) -> Iterator[AnnotatedVariant]:
        """Attach gene context to each variant in the stream.

        When inheritance_mode is set, variants in genes that don't
        match the specified mode are filtered out. Intergenic variants
        (gene_name=None) pass through regardless.

        Parameters
        ----------
        variants : Iterator[AnnotatedVariant]
            Input stream of annotated variants (gene_name already resolved).

        Yields
        ------
        AnnotatedVariant
            Variants with gene_context field populated, filtered by
            inheritance mode when configured.
        """
        for variant in variants:
            gene_context = self._registry.build_gene_context(variant.gene_name)

            # Inheritance mode filtering: drop variants in genes that
            # don't match the requested mode. Intergenic variants and
            # genes with no OMIM data pass through unfiltered.
            if self._inheritance_mode is not None and variant.gene_name is not None:
                gene_modes = self._registry.omim.get_inheritance_modes(
                    variant.gene_name
                )
                if gene_modes and self._inheritance_mode not in gene_modes:
                    continue

            yield replace(variant, gene_context=gene_context)
