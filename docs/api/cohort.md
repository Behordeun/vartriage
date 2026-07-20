# Cohort Analysis

Multi-sample cohort analysis classes for cross-sample variant aggregation and reporting.

## CohortPipeline

::: vartriage.cohort.pipeline.CohortPipeline

Orchestrates per-sample Pipeline execution, aggregation, statistics computation, and report generation.

```python
from vartriage import CohortPipeline, CohortConfig

config = CohortConfig(
    sample_vcfs=[Path("a.vcf.gz"), Path("b.vcf.gz")],
    output_path=Path("output/"),
)
pipeline = CohortPipeline(cohort_config=config)
report_paths = pipeline.run()

# Access results after run()
pipeline.variants       # list[CohortVariant]
pipeline.gene_burdens   # list[GeneBurden]
pipeline.summary        # CohortSummary
```

**Parameters:**

| Name | Type | Description |
| --- | --- | --- |
| cohort_config | CohortConfig | Cohort-level settings |
| pipeline_config | PipelineConfig, optional | Base config applied to each sample |
| annotation_config | AnnotationConfig, optional | Shared annotation config (overrides pipeline_config) |
| prioritization_config | PrioritizationConfig, optional | Shared scoring config (overrides pipeline_config) |

**Methods:**

- `run() -> list[Path]` - Execute the full cohort analysis, returns paths to generated reports.

**Properties (populated after run):**

- `variants: list[CohortVariant]`
- `gene_burdens: list[GeneBurden]`
- `summary: CohortSummary | None`

---

## CohortAggregator

::: vartriage.cohort.aggregator.CohortAggregator

Merges classified variants from multiple samples by genomic coordinate.

**Methods:**

- `add_sample(sample_id, vcf_path, variants) -> int` - Ingest one sample's results. Returns count after AF filtering.
- `aggregate() -> list[CohortVariant]` - Produce merged cohort variants respecting config thresholds.
- `get_recurrent_variants(min_count=2) -> list[CohortVariant]` - Convenience filter for shared variants.
- `reset()` - Clear all data for reuse.

**Properties:**

- `samples_added: int`
- `total_distinct_variants: int`

---

## CohortStatistics

::: vartriage.cohort.statistics.CohortStatistics

Computes summary metrics from aggregated cohort variants.

**Methods:**

- `compute_summary(samples_processed) -> CohortSummary` - Top-level metrics.
- `compute_gene_burden() -> list[GeneBurden]` - Per-gene mutation burden sorted by severity.
- `recurrence_distribution() -> dict[int, int]` - Sample count histogram.
- `per_sample_counts() -> dict[str, int]` - Variants per sample.
- `classification_distribution() -> dict[str, int]` - Count by ACMG class.
- `consequence_distribution() -> dict[str, int]` - Count by consequence type.

---

## CohortReportGenerator

::: vartriage.cohort.report.CohortReportGenerator

Writes cohort analysis results to disk in JSON or CSV format.

**Methods:**

- `generate(variants, gene_burdens, summary) -> list[Path]` - Write all report files, returns paths.

---

## Data Models

### CohortConfig

Frozen dataclass. Configuration for multi-sample cohort analysis. See [Cohort Analysis Guide](../cohort-analysis.md#cohortconfig-reference) for the full field reference.

### CohortVariant

Frozen dataclass. A variant aggregated across multiple samples.

**Key fields:** `chrom`, `pos`, `ref`, `alt`, `gene_name`, `consequence`, `sample_count`, `total_samples`, `occurrences`, `max_classification`, `all_evidence_tags`, `allele_frequency`.

**Properties:** `key`, `cohort_frequency`, `is_singleton`, `is_universal`, `sample_ids`.

### SampleOccurrence

Frozen dataclass. Record of a variant's appearance in one sample.

**Fields:** `sample_id`, `vcf_path`, `classified` (ClassifiedVariant).

### GeneBurden

Frozen dataclass. Per-gene variant burden across the cohort.

**Fields:** `gene_name`, `total_variants`, `pathogenic_count`, `samples_affected`, `total_samples`, `most_severe`.

**Properties:** `penetrance` (float, fraction of cohort affected).

### CohortSummary

Frozen dataclass. Aggregate statistics for a completed cohort run.

**Fields:** `cohort_name`, `total_samples`, `total_variants`, `shared_variants`, `singleton_variants`, `universal_variants`, `pathogenic_variants`, `likely_pathogenic_variants`, `genes_affected`, `top_recurrent_genes`, `samples_processed`.
