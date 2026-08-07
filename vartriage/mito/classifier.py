"""Mitochondrial variant classifier using mtDNA-specific criteria.

Applies classification rules distinct from the nuclear ACMG/AMP 2015
framework. mtDNA variants are classified based on:
- MITOMAP confirmation status
- Heteroplasmy level
- HelixMTdb population frequency
- Gene context (protein-coding vs tRNA vs rRNA vs control region)

Classification hierarchy:
1. Pathogenic: confirmed MITOMAP + high/homoplasmic heteroplasmy + rare
2. Benign: common haplogroup marker (AF > 5%)
3. Likely Benign: moderate frequency (AF > 0.1%)
4. Likely Pathogenic: reported MITOMAP + functional region + moderate-high level
5. VUS: default for novel variants without strong evidence
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import Enum

from vartriage.mito.frequency import HelixMTdbDatabase
from vartriage.mito.gene_map import MtGeneContext, MtGeneMap
from vartriage.mito.heteroplasmy import HeteroplasmyLevel, extract_heteroplasmy
from vartriage.mito.mitomap import MitomapDatabase, MitomapEntry
from vartriage.models.variant import Variant

logger = logging.getLogger(__name__)

# Frequency thresholds for classification
_BENIGN_AF_THRESHOLD: float = 0.05
_LIKELY_BENIGN_AF_THRESHOLD: float = 0.001
_RARE_AF_THRESHOLD: float = 0.0001


class MitoClassification(Enum):
    """Final classification for a mitochondrial variant."""

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_Pathogenic"
    VUS = "VUS"
    LIKELY_BENIGN = "Likely_Benign"
    BENIGN = "Benign"


@dataclass(frozen=True, slots=True)
class MitoClassifiedVariant:
    """A mitochondrial variant with full annotation and classification.

    Parameters
    ----------
    variant
        Original VCF variant record.
    classification
        mtDNA-specific classification result.
    heteroplasmy
        Heteroplasmy level, or None if extraction failed.
    gene_context
        MT gene map annotation for this position.
    mitomap_entry
        MITOMAP disease association, or None if not found.
    helix_af
        HelixMTdb population allele frequency, or None if novel.
    classification_reason
        Human-readable explanation of why this classification was assigned.
    """

    variant: Variant
    classification: MitoClassification
    heteroplasmy: HeteroplasmyLevel | None
    gene_context: MtGeneContext
    mitomap_entry: MitomapEntry | None
    helix_af: float | None
    classification_reason: str


class MitochondrialClassifier:
    """Classify mtDNA variants using mitochondrial-specific criteria.

    Composes gene map, MITOMAP, and HelixMTdb lookups with heteroplasmy
    extraction to produce a final classification per variant.

    Parameters
    ----------
    gene_map
        Pre-loaded mitochondrial gene interval index.
    mitomap_db
        Pre-loaded MITOMAP pathogenic mutation database.
    helix_db
        Pre-loaded HelixMTdb population frequency database.
    """

    def __init__(
        self,
        gene_map: MtGeneMap,
        mitomap_db: MitomapDatabase,
        helix_db: HelixMTdbDatabase,
    ) -> None:
        self._gene_map = gene_map
        self._mitomap_db = mitomap_db
        self._helix_db = helix_db

    def classify(self, variant: Variant) -> MitoClassifiedVariant:
        """Classify a single mitochondrial variant.

        Parameters
        ----------
        variant
            Raw variant record (must be chrM/MT).

        Returns
        -------
        MitoClassifiedVariant
            Fully annotated and classified mitochondrial variant.
        """
        # Extract heteroplasmy from the variant info dict
        heteroplasmy = extract_heteroplasmy(variant.info)

        # Gene context lookup
        gene_context = self._gene_map.query(variant.pos)

        # MITOMAP lookup
        mitomap_entry = self._mitomap_db.lookup(variant.pos, variant.ref, variant.alt)

        # HelixMTdb frequency lookup
        helix_af = self._helix_db.get_af(variant.pos, variant.ref, variant.alt)

        # Apply classification logic
        classification, reason = self._apply_rules(
            heteroplasmy=heteroplasmy,
            gene_context=gene_context,
            mitomap_entry=mitomap_entry,
            helix_af=helix_af,
        )

        return MitoClassifiedVariant(
            variant=variant,
            classification=classification,
            heteroplasmy=heteroplasmy,
            gene_context=gene_context,
            mitomap_entry=mitomap_entry,
            helix_af=helix_af,
            classification_reason=reason,
        )

    def _apply_rules(
        self,
        heteroplasmy: HeteroplasmyLevel | None,
        gene_context: MtGeneContext,
        mitomap_entry: MitomapEntry | None,
        helix_af: float | None,
    ) -> tuple[MitoClassification, str]:
        """Apply the classification decision tree.

        Returns a (classification, reason) tuple. Rules are evaluated
        in priority order — first match wins.
        """
        is_rare = helix_af is None or helix_af < _RARE_AF_THRESHOLD
        has_high_heteroplasmy = heteroplasmy is not None and heteroplasmy.category in (
            "high",
            "homoplasmic",
        )
        has_moderate_or_high = heteroplasmy is not None and heteroplasmy.category in (
            "high",
            "homoplasmic",
            "moderate",
        )

        # Rule 1: Pathogenic
        # Confirmed in MITOMAP + high/homoplasmic heteroplasmy + rare
        if (
            mitomap_entry is not None
            and mitomap_entry.is_confirmed
            and has_high_heteroplasmy
            and is_rare
        ):
            return (
                MitoClassification.PATHOGENIC,
                f"Confirmed pathogenic in MITOMAP ({mitomap_entry.disease}), "
                f"heteroplasmy {heteroplasmy.percentage:.1f}% "  # type: ignore[union-attr]
                f"({heteroplasmy.category}), rare in population",  # type: ignore[union-attr]
            )

        # Rule 2: Benign — common haplogroup marker
        if helix_af is not None and helix_af > _BENIGN_AF_THRESHOLD:
            return (
                MitoClassification.BENIGN,
                f"Common haplogroup-defining polymorphism "
                f"(HelixMTdb AF={helix_af:.4f})",
            )

        # Rule 3: Likely Benign — moderate frequency
        if helix_af is not None and helix_af > _LIKELY_BENIGN_AF_THRESHOLD:
            return (
                MitoClassification.LIKELY_BENIGN,
                f"Moderate population frequency (HelixMTdb AF={helix_af:.4f})",
            )

        # Rule 4: Likely Pathogenic
        # Reported in MITOMAP + functional region + moderate-high heteroplasmy
        if (
            mitomap_entry is not None
            and gene_context.is_in_coding_or_trna
            and has_moderate_or_high
        ):
            return (
                MitoClassification.LIKELY_PATHOGENIC,
                f"Reported in MITOMAP ({mitomap_entry.disease}, "
                f"status={mitomap_entry.status}), "
                f"in {gene_context.gene_type} region ({gene_context.gene_name}), "
                f"heteroplasmy {heteroplasmy.percentage:.1f}%",  # type: ignore[union-attr]
            )

        # Rule 5: Pathogenic without heteroplasmy data
        # Confirmed MITOMAP + rare, but no heteroplasmy measurement
        if (
            mitomap_entry is not None
            and mitomap_entry.is_confirmed
            and is_rare
            and heteroplasmy is None
        ):
            return (
                MitoClassification.LIKELY_PATHOGENIC,
                f"Confirmed pathogenic in MITOMAP ({mitomap_entry.disease}), "
                f"rare in population, heteroplasmy unavailable",
            )

        # Default: VUS
        reason_parts = ["Novel or insufficient evidence"]
        if gene_context.gene_name is not None:
            reason_parts.append(f"in {gene_context.gene_name}")
        if helix_af is None:
            reason_parts.append("absent from HelixMTdb")
        if heteroplasmy is not None:
            reason_parts.append(f"heteroplasmy {heteroplasmy.percentage:.1f}%")

        return (MitoClassification.VUS, ", ".join(reason_parts))
