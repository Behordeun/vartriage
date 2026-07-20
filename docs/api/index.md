# API Reference

Auto-generated documentation from source code docstrings.

## Pipeline orchestration

- [Pipeline](pipeline.md) - Top-level orchestrator connecting all stages
- [CohortPipeline](cohort.md) - Multi-sample cohort analysis orchestrator

## Processing stages

- [VCFParser](parser.md) - VCF file streaming
- [QualityFilter](filter.md) - Quality-based variant exclusion
- [AnnotationEngine](annotation.md) - Functional consequence and population data
- [PrioritizationEngine](prioritization.md) - Frequency filtering and scoring
- [ACMGClassifier](classification.md) - Evidence tag assignment and classification
- [ReportGenerator](reporting.md) - Output serialization (JSON, CSV, PDF)

## Cohort analysis

- [CohortPipeline](cohort.md) - Multi-sample orchestrator
- [CohortAggregator](cohort.md#cohortaggregator) - Cross-sample variant merging
- [CohortStatistics](cohort.md#cohortstatistics) - Per-gene burden and recurrence stats
- [CohortReportGenerator](cohort.md#cohortreportgenerator) - Cohort report generation

## Data types

- [Data Models](models.md) - Variant, AnnotatedVariant, ScoredVariant, ClassifiedVariant, enums
- [Cohort Models](cohort.md#data-models) - CohortConfig, CohortVariant, GeneBurden, CohortSummary
- [Configuration Classes](config.md) - All config dataclasses
- [Exceptions](exceptions.md) - Error hierarchy
