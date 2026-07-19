"""Cross-sample variant aggregation for cohort analysis.

Merges classified variants from multiple samples by genomic coordinate,
computes cohort-level recurrence frequencies, and identifies shared vs
unique variants. The merge key is (chrom, pos, ref, alt) — identical
coordinates across samples are considered the same variant.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Optional

from vartriage.models.cohort import (
    CohortConfig,
    CohortVariant,
    SampleOccurrence,
    VariantKey,
)
from vartriage.models.variant import (
    ACMGClassification,
    ClassifiedVariant,
    EvidenceTag,
    FunctionalConsequence,
)

logger = logging.getLogger(__name__)

# Severity ordering for ACMG classifications (index 0 = most severe)
_CLASSIFICATION_SEVERITY: list[ACMGClassification] = [
    ACMGClassification.PATHOGENIC,
    ACMGClassification.LIKELY_PATHOGENIC,
    ACMGClassification.VUS,
    ACMGClassification.LIKELY_BENIGN,
    ACMGClassification.BENIGN,
]

# Consequence severity (mirrors models/variant.py ordering)
_CONSEQUENCE_SEVERITY: list[FunctionalConsequence] = [
    FunctionalConsequence.FRAMESHIFT,
    FunctionalConsequence.NONSENSE,
    FunctionalConsequence.SPLICE_SITE,
    FunctionalConsequence.MISSENSE,
    FunctionalConsequence.IN_FRAME_INSERTION,
    FunctionalConsequence.IN_FRAME_DELETION,
    FunctionalConsequence.SYNONYMOUS,
    FunctionalConsequence.INTERGENIC,
]


def _most_severe_classification(
    classifications: list[ACMGClassification],
) -> ACMGClassification:
    """Return the most severe classification from a list.

    Falls back to VUS if the list is empty (should not happen in practice).
    """
    if not classifications:
        return ACMGClassification.VUS

    best_idx = len(_CLASSIFICATION_SEVERITY)
    for cls in classifications:
        try:
            idx = _CLASSIFICATION_SEVERITY.index(cls)
        except ValueError:
            continue
        if idx < best_idx:
            best_idx = idx

    if best_idx >= len(_CLASSIFICATION_SEVERITY):
        return ACMGClassification.VUS
    return _CLASSIFICATION_SEVERITY[best_idx]


def _most_severe_consequence(
    consequences: list[FunctionalConsequence],
) -> FunctionalConsequence:
    """Return the most severe consequence from a list."""
    if not consequences:
        return FunctionalConsequence.INTERGENIC

    best_idx = len(_CONSEQUENCE_SEVERITY)
    for csq in consequences:
        try:
            idx = _CONSEQUENCE_SEVERITY.index(csq)
        except ValueError:
            continue
        if idx < best_idx:
            best_idx = idx

    if best_idx >= len(_CONSEQUENCE_SEVERITY):
        return FunctionalConsequence.INTERGENIC
    return _CONSEQUENCE_SEVERITY[best_idx]


class CohortAggregator:
    """Merges per-sample classified variants into cohort-level records.

    Groups variants by genomic coordinate across all samples, then
    produces CohortVariant records with cross-sample frequency and
    merged evidence. Respects the config's min_recurrence threshold
    and max_af_threshold for filtering.

    Parameters
    ----------
    config : CohortConfig
        Cohort analysis configuration.
    """

    def __init__(self, config: CohortConfig) -> None:
        self._config = config
        self._variant_map: dict[VariantKey, list[SampleOccurrence]] = defaultdict(list)
        self._samples_added: int = 0

    @property
    def samples_added(self) -> int:
        """Number of samples that have been ingested so far."""
        return self._samples_added

    @property
    def total_distinct_variants(self) -> int:
        """Number of distinct variant coordinates seen across all samples."""
        return len(self._variant_map)

    def add_sample(
        self,
        sample_id: str,
        vcf_path: Path,
        variants: list[ClassifiedVariant],
    ) -> int:
        """Ingest classified variants from a single sample.

        Parameters
        ----------
        sample_id : str
            Human-readable sample identifier.
        vcf_path : Path
            Source VCF path for traceability.
        variants : list[ClassifiedVariant]
            All classified variants from this sample's pipeline run.

        Returns
        -------
        int
            Number of variants added from this sample (after AF filtering).
        """
        added = 0
        max_af = self._config.max_af_threshold

        for classified in variants:
            af = classified.scored.annotated.allele_frequency

            # Skip common variants that exceed cohort AF threshold
            if af is not None and af > max_af:
                continue

            key = self._variant_key(classified)
            occurrence = SampleOccurrence(
                sample_id=sample_id,
                vcf_path=vcf_path,
                classified=classified,
            )
            self._variant_map[key].append(occurrence)
            added += 1

        self._samples_added += 1
        logger.info(
            "Added sample '%s': %d variants (of %d total, %d filtered by AF)",
            sample_id,
            added,
            len(variants),
            len(variants) - added,
        )
        return added

    def aggregate(self) -> list[CohortVariant]:
        """Produce merged cohort variants from all ingested samples.

        Filtering logic:
        - Variants with sample_count >= min_recurrence always pass.
        - Singletons (sample_count == 1) pass only when include_singletons
          is True.
        - Variants with 1 < sample_count < min_recurrence are included
          (they are shared but below the recurrence highlight threshold).

        Results are sorted by sample_count descending, then by genomic
        coordinate.

        Returns
        -------
        list[CohortVariant]
            Cohort-level variant records meeting inclusion criteria.
        """
        total_samples = self._config.sample_count
        min_rec = self._config.min_recurrence
        include_singletons = self._config.include_singletons

        results: list[CohortVariant] = []

        for key, occurrences in self._variant_map.items():
            sample_count = len(occurrences)

            # Singletons are excluded unless explicitly included
            if sample_count == 1 and not include_singletons:
                continue

            # Variants below min_recurrence are excluded (unless they
            # are singletons that passed the check above)
            if sample_count < min_rec and sample_count > 1:
                continue

            cohort_variant = self._build_cohort_variant(
                key, occurrences, total_samples
            )
            results.append(cohort_variant)

        # Sort: highest recurrence first, then by coordinate for stability
        results.sort(
            key=lambda v: (-v.sample_count, v.chrom, v.pos, v.ref, v.alt)
        )

        logger.info(
            "Aggregation complete: %d variants (%d shared, %d singletons)",
            len(results),
            sum(1 for v in results if v.sample_count >= 2),
            sum(1 for v in results if v.is_singleton),
        )
        return results

    def get_recurrent_variants(self, min_count: int = 2) -> list[CohortVariant]:
        """Return only variants appearing in >= min_count samples.

        Convenience method for extracting shared variants without
        changing the config's min_recurrence permanently.

        Parameters
        ----------
        min_count : int
            Minimum sample count threshold. Default is 2.

        Returns
        -------
        list[CohortVariant]
            Filtered and sorted cohort variants.
        """
        all_variants = self.aggregate()
        return [v for v in all_variants if v.sample_count >= min_count]

    def reset(self) -> None:
        """Clear all ingested data for reuse with a fresh cohort."""
        self._variant_map.clear()
        self._samples_added = 0

    def _build_cohort_variant(
        self,
        key: VariantKey,
        occurrences: list[SampleOccurrence],
        total_samples: int,
    ) -> CohortVariant:
        """Construct a CohortVariant from merged sample occurrences."""
        chrom, pos, ref, alt = key

        # Collect classifications and consequences across samples
        classifications = [
            occ.classified.classification for occ in occurrences
        ]
        consequences = [
            occ.classified.scored.annotated.consequence for occ in occurrences
        ]

        # Union of all evidence tags
        all_tags: set[EvidenceTag] = set()
        for occ in occurrences:
            all_tags.update(occ.classified.evidence_tags)

        # Gene name: take the first non-None value
        gene_name: Optional[str] = None
        for occ in occurrences:
            gn = occ.classified.scored.annotated.gene_name
            if gn is not None:
                gene_name = gn
                break

        # Allele frequency: consensus (use first non-None)
        allele_frequency: Optional[float] = None
        for occ in occurrences:
            af = occ.classified.scored.annotated.allele_frequency
            if af is not None:
                allele_frequency = af
                break

        return CohortVariant(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            gene_name=gene_name,
            consequence=_most_severe_consequence(consequences),
            sample_count=len(occurrences),
            total_samples=total_samples,
            occurrences=tuple(occurrences),
            max_classification=_most_severe_classification(classifications),
            all_evidence_tags=frozenset(all_tags),
            allele_frequency=allele_frequency,
        )

    @staticmethod
    def _variant_key(classified: ClassifiedVariant) -> VariantKey:
        """Extract the canonical merge key from a classified variant."""
        v = classified.scored.annotated.variant
        return (v.chrom, v.pos, v.ref, v.alt)
