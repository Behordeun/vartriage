"""Multi-sample cohort analysis for variant prioritization.

Provides cross-sample aggregation, recurrence frequency computation,
per-gene burden analysis, and cohort-level reporting. Designed to
process multiple VCF files (one per sample) through the standard
pipeline then merge results for population-level insights.
"""

from vartriage.cohort.aggregator import CohortAggregator
from vartriage.cohort.pipeline import CohortPipeline
from vartriage.cohort.report import CohortReportGenerator
from vartriage.cohort.statistics import CohortStatistics

__all__ = [
    "CohortAggregator",
    "CohortPipeline",
    "CohortReportGenerator",
    "CohortStatistics",
]
