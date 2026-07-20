"""Cohort-level statistics computation.

Computes per-gene burden, sample-level variant counts, recurrence
distribution, and the CohortSummary from aggregated cohort variants.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from typing import Optional

from vartriage.models.cohort import (
    CohortConfig,
    CohortSummary,
    CohortVariant,
    GeneBurden,
)
from vartriage.models.variant import (
    ACMGClassification,
    CLASSIFICATION_SEVERITY_ORDER,
)

logger = logging.getLogger(__name__)

_PATHOGENIC_CLASSIFICATIONS: frozenset[ACMGClassification] = frozenset(
    {
        ACMGClassification.PATHOGENIC,
        ACMGClassification.LIKELY_PATHOGENIC,
    }
)


class CohortStatistics:
    """Compute summary statistics from aggregated cohort variants.

    Takes a list of CohortVariant records (output of CohortAggregator)
    and produces per-gene burden tables, recurrence distributions, and
    the overall CohortSummary dataclass.

    Parameters
    ----------
    config : CohortConfig
        Cohort configuration (used for cohort_name and sample metadata).
    variants : list[CohortVariant]
        Aggregated cohort variants to analyze.
    """

    def __init__(
        self, config: CohortConfig, variants: list[CohortVariant]
    ) -> None:
        self._config = config
        self._variants = variants

    @property
    def variant_count(self) -> int:
        """Total distinct variants in the cohort."""
        return len(self._variants)

    def compute_summary(
        self, samples_processed: list[str]
    ) -> CohortSummary:
        """Produce the top-level cohort summary.

        Parameters
        ----------
        samples_processed : list[str]
            Ordered list of sample identifiers that were analyzed.

        Returns
        -------
        CohortSummary
            Aggregate statistics for the cohort run.
        """
        total = len(self._variants)
        shared = sum(1 for v in self._variants if v.sample_count >= 2)
        singletons = sum(1 for v in self._variants if v.is_singleton)
        universal = sum(1 for v in self._variants if v.is_universal)

        pathogenic = sum(
            1
            for v in self._variants
            if v.max_classification == ACMGClassification.PATHOGENIC
        )
        likely_pathogenic = sum(
            1
            for v in self._variants
            if v.max_classification == ACMGClassification.LIKELY_PATHOGENIC
        )

        genes: set[str] = set()
        for v in self._variants:
            if v.gene_name is not None:
                genes.add(v.gene_name)

        top_genes = self._top_recurrent_genes(limit=10)

        return CohortSummary(
            cohort_name=self._config.cohort_name,
            total_samples=self._config.sample_count,
            total_variants=total,
            shared_variants=shared,
            singleton_variants=singletons,
            universal_variants=universal,
            pathogenic_variants=pathogenic,
            likely_pathogenic_variants=likely_pathogenic,
            genes_affected=len(genes),
            top_recurrent_genes=tuple(top_genes),
            samples_processed=tuple(samples_processed),
        )

    def compute_gene_burden(self) -> list[GeneBurden]:
        """Compute per-gene variant burden across the cohort.

        Groups variants by gene, counts pathogenic/likely_pathogenic
        hits, and tracks how many samples are affected per gene.
        Results are sorted by pathogenic_count descending, then by
        samples_affected descending.

        Returns
        -------
        list[GeneBurden]
            Per-gene burden records, sorted by severity.
        """
        # gene -> list of cohort variants
        gene_variants: dict[str, list[CohortVariant]] = defaultdict(list)

        for v in self._variants:
            if v.gene_name is None:
                continue
            gene_variants[v.gene_name].append(v)

        total_samples = self._config.sample_count
        burdens: list[GeneBurden] = []

        for gene_name, variants in gene_variants.items():
            pathogenic_count = sum(
                1
                for v in variants
                if v.max_classification in _PATHOGENIC_CLASSIFICATIONS
            )

            # Unique samples affected in this gene
            samples_in_gene: set[str] = set()
            for v in variants:
                samples_in_gene.update(v.sample_ids)

            # Most severe classification in this gene
            most_severe = self._gene_most_severe(variants)

            burdens.append(
                GeneBurden(
                    gene_name=gene_name,
                    total_variants=len(variants),
                    pathogenic_count=pathogenic_count,
                    samples_affected=len(samples_in_gene),
                    total_samples=total_samples,
                    most_severe=most_severe,
                )
            )

        burdens.sort(
            key=lambda b: (-b.pathogenic_count, -b.samples_affected, b.gene_name)
        )

        logger.info(
            "Gene burden computed: %d genes, %d with pathogenic variants",
            len(burdens),
            sum(1 for b in burdens if b.pathogenic_count > 0),
        )
        return burdens

    def recurrence_distribution(self) -> dict[int, int]:
        """Count how many variants appear in exactly N samples.

        Returns
        -------
        dict[int, int]
            Mapping of sample_count -> number of variants with that count.
            Sorted by key ascending.
        """
        counter: Counter[int] = Counter()
        for v in self._variants:
            counter[v.sample_count] += 1
        return dict(sorted(counter.items()))

    def per_sample_counts(self) -> dict[str, int]:
        """Count total variants per sample across the cohort.

        Returns
        -------
        dict[str, int]
            Mapping of sample_id -> number of cohort variants that
            include that sample. Sorted by count descending.
        """
        counts: Counter[str] = Counter()
        for v in self._variants:
            for sample_id in v.sample_ids:
                counts[sample_id] += 1
        return dict(counts.most_common())

    def classification_distribution(self) -> dict[str, int]:
        """Count variants by their most severe ACMG classification.

        Returns
        -------
        dict[str, int]
            Mapping of classification name -> count.
        """
        counter: Counter[str] = Counter()
        for v in self._variants:
            counter[v.max_classification.value] += 1
        return dict(counter.most_common())

    def consequence_distribution(self) -> dict[str, int]:
        """Count variants by functional consequence type.

        Returns
        -------
        dict[str, int]
            Mapping of consequence name -> count.
        """
        counter: Counter[str] = Counter()
        for v in self._variants:
            counter[v.consequence.value] += 1
        return dict(counter.most_common())

    def _top_recurrent_genes(self, limit: int = 10) -> list[str]:
        """Get genes with the most recurrent variants.

        Ranks genes by the sum of sample_count across all their
        variants, which captures both how many variants a gene has
        and how widely shared they are.
        """
        gene_recurrence: Counter[str] = Counter()
        for v in self._variants:
            if v.gene_name is not None:
                gene_recurrence[v.gene_name] += v.sample_count
        return [gene for gene, _ in gene_recurrence.most_common(limit)]

    @staticmethod
    def _gene_most_severe(variants: list[CohortVariant]) -> ACMGClassification:
        """Find the most severe classification across gene variants."""
        best_idx = len(CLASSIFICATION_SEVERITY_ORDER) - 1
        for v in variants:
            try:
                idx = CLASSIFICATION_SEVERITY_ORDER.index(v.max_classification)
            except ValueError:
                continue
            if idx < best_idx:
                best_idx = idx
        return CLASSIFICATION_SEVERITY_ORDER[best_idx]
