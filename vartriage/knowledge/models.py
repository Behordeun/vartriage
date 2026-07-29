"""Data models for gene-disease linkage knowledge base.

Immutable dataclasses representing gene-level annotations from OMIM,
ClinGen, HPO, and gnomAD constraint databases.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True, slots=True)
class DiseaseAssociation:
    """A single gene-disease relationship from OMIM.

    Parameters
    ----------
    disease_name : str
        Human-readable disease name (e.g., "Breast-ovarian cancer, familial, 1").
    mim_number : str
        OMIM MIM number for the phenotype entry.
    inheritance_mode : str
        Inheritance pattern code: AD, AR, XL, XLD, XLR, or MT.
    """

    disease_name: str
    mim_number: str
    inheritance_mode: str


@dataclass(frozen=True, slots=True)
class GeneConstraint:
    """gnomAD gene constraint metrics indicating intolerance to variation.

    Parameters
    ----------
    pli : float
        Probability of loss-of-function intolerance (0.0-1.0).
        Values >0.9 suggest haploinsufficiency.
    loeuf : float
        Loss-of-function observed/expected upper bound fraction.
        Lower values indicate stronger constraint.
    mis_z : float
        Missense Z-score. Higher values indicate intolerance to
        missense variation.
    """

    pli: float
    loeuf: float
    mis_z: float

    @property
    def is_lof_intolerant(self) -> bool:
        """True if gene is highly intolerant to loss-of-function (pLI > 0.9)."""
        return self.pli > 0.9

    @property
    def is_missense_constrained(self) -> bool:
        """True if gene is intolerant to missense variation (mis_z > 3.09)."""
        return self.mis_z > 3.09


@dataclass(frozen=True, slots=True)
class GeneAnnotation:
    """Complete gene-level annotation aggregated from all knowledge sources.

    Returned by GeneKnowledgeRegistry.annotate_gene(). This is the
    internal representation used by the registry cache.

    Parameters
    ----------
    disease_associations : tuple[DiseaseAssociation, ...]
        All known disease relationships from OMIM.
    clingen_validity : str | None
        ClinGen gene-disease validity level (Definitive, Strong,
        Moderate, Limited, Disputed, Refuted) or None if not curated.
    constraint : GeneConstraint | None
        gnomAD constraint metrics or None if unavailable.
    is_actionable : bool
        True if the gene has ClinGen actionability curations.
    actionability_type : str | None
        Type of intervention (surveillance, therapeutic, etc.) or None.
    hpo_terms : frozenset[str]
        HPO term IDs associated with this gene's phenotypes.
    """

    disease_associations: tuple[DiseaseAssociation, ...]
    clingen_validity: Optional[str] = None
    constraint: Optional[GeneConstraint] = None
    is_actionable: bool = False
    actionability_type: Optional[str] = None
    hpo_terms: frozenset[str] = frozenset()


# Neutral default for genes missing from all knowledge sources
EMPTY_GENE_ANNOTATION = GeneAnnotation(
    disease_associations=(),
    clingen_validity=None,
    constraint=None,
    is_actionable=False,
    actionability_type=None,
    hpo_terms=frozenset(),
)


@dataclass(frozen=True, slots=True)
class GeneContext:
    """Gene-level context attached to each variant in the output.

    This is the variant-facing view of gene knowledge, including a
    phenotype match score computed per-run based on patient HPO terms.

    Parameters
    ----------
    disease_associations : tuple[DiseaseAssociation, ...]
        Known disease relationships for the gene.
    clingen_validity : str | None
        ClinGen validity level or None.
    constraint : GeneConstraint | None
        gnomAD constraint metrics or None.
    is_actionable : bool
        Whether the gene has established medical interventions.
    phenotype_match_score : float
        Fraction of patient HPO terms matching this gene (0.0-1.0).
        Defaults to 0.0 when no patient HPO terms are provided.
    """

    disease_associations: tuple[DiseaseAssociation, ...]
    clingen_validity: Optional[str] = None
    constraint: Optional[GeneConstraint] = None
    is_actionable: bool = False
    phenotype_match_score: float = 0.0
