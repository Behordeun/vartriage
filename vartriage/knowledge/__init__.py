"""Gene-disease linkage knowledge base.

Provides gene-level annotations from OMIM, ClinGen, HPO, and gnomAD
constraint data. All lookups are O(1) dict-based after initial load.
"""

from vartriage.knowledge.config import KnowledgeBaseConfig
from vartriage.knowledge.models import (
    DiseaseAssociation,
    GeneAnnotation,
    GeneConstraint,
    GeneContext,
)
from vartriage.knowledge.registry import GeneKnowledgeRegistry

__all__ = [
    "DiseaseAssociation",
    "GeneAnnotation",
    "GeneConstraint",
    "GeneContext",
    "GeneKnowledgeRegistry",
    "KnowledgeBaseConfig",
]
