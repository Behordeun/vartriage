# API Reference

Auto-generated documentation from source code docstrings.

## Pipeline orchestration

- [Pipeline](pipeline.md) - Top-level orchestrator connecting all stages

## Processing stages

- [VCFParser](parser.md) - VCF file streaming
- [QualityFilter](filter.md) - Quality-based variant exclusion
- [AnnotationEngine](annotation.md) - Functional consequence and population data
- [PrioritizationEngine](prioritization.md) - Frequency filtering and scoring
- [ACMGClassifier](classification.md) - Evidence tag assignment and classification
- [ReportGenerator](reporting.md) - Output serialization (JSON, CSV, PDF)

## Data types

- [Data Models](models.md) - Variant, AnnotatedVariant, ScoredVariant, ClassifiedVariant, enums
- [Configuration Classes](config.md) - All config dataclasses
- [Exceptions](exceptions.md) - Error hierarchy
