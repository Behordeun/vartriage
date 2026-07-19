# Cohort Analysis

Analyze multiple samples together to identify shared variants, compute recurrence frequencies, and quantify per-gene mutation burden across a cohort.

## When to use cohort mode

- Identifying recurrent mutations across a patient group (e.g., 30 cardiac patients sharing the same MYBPC3 variant)
- Comparing variant landscapes between affected and unaffected individuals
- Building frequency tables for lab-internal population data
- Screening for founder mutations in geographically clustered cohorts

## CLI usage

### Direct VCF arguments

```bash
vartriage cohort \
  --vcf sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
  --output cohort_results/ \
  --cohort-name "cardiac_cohort" \
  --output-format json
```

### Manifest file

For larger cohorts, list VCF paths in a manifest file:

```bash
vartriage cohort \
  --manifest samples.tsv \
  --output cohort_results/ \
  --cohort-name "cardiac_cohort"
```

Manifest format: one VCF path per line. Optional tab-separated second column provides a human-readable label. Lines starting with `#` are comments. Blank lines are skipped.

```text
# Cardiac cohort, sequenced 2026-Q2
/data/vcfs/patient_001.vcf.gz	Patient 001
/data/vcfs/patient_002.vcf.gz	Patient 002
/data/vcfs/patient_003.vcf.gz	Patient 003
```

Relative paths in the manifest are resolved relative to the manifest file's directory.

### With annotation and filtering

```bash
vartriage cohort \
  --manifest samples.tsv \
  --output cohort_results/ \
  --cohort-name "cardiac_cohort" \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv \
  --cadd-scores cadd.tsv \
  --revel-scores revel.tsv \
  --spliceai-scores spliceai.tsv \
  --gene-list cardiac_panel.txt \
  --use-bundles
```

All reference file flags work the same as the standard single-sample pipeline. They apply uniformly across every sample in the cohort.

### Parallel processing

```bash
vartriage cohort \
  --manifest samples.tsv \
  --output cohort_results/ \
  --parallel \
  --max-workers 8
```

Parallel mode uses a thread pool. Speedup depends on I/O throughput since per-sample processing is dominated by VCF parsing and reference file reads.

### Filtering options

```bash
vartriage cohort \
  --manifest samples.tsv \
  --output cohort_results/ \
  --min-recurrence 3 \
  --max-af 0.01 \
  --no-singletons
```

| Flag | Default | Effect |
| --- | --- | --- |
| `--min-recurrence` | 2 | Variants below this sample count are still output but not flagged as recurrent |
| `--max-af` | 0.05 | Variants with gnomAD AF above this are excluded from aggregation entirely |
| `--no-singletons` | false | Drop variants appearing in only one sample from the output |
| `--output-format` | json | Choose `json` or `csv` |

## Python API

### Basic usage

```python
from pathlib import Path
from vartriage import CohortPipeline, CohortConfig

cohort_config = CohortConfig(
    sample_vcfs=[
        Path("patient_001.vcf.gz"),
        Path("patient_002.vcf.gz"),
        Path("patient_003.vcf.gz"),
    ],
    output_path=Path("cohort_results/"),
    cohort_name="cardiac_cohort",
    min_recurrence=2,
    max_af_threshold=0.05,
)

pipeline = CohortPipeline(cohort_config=cohort_config)
report_paths = pipeline.run()
```

### With shared annotation

```python
from pathlib import Path
from vartriage import (
    CohortPipeline, CohortConfig,
    AnnotationConfig, PrioritizationConfig,
)

annotation = AnnotationConfig(
    gene_annotation_path=Path("gencode.v44.gtf"),
    gnomad_path=Path("gnomad.v4.sites.tsv"),
    clinvar_path=Path("clinvar.tsv"),
)

prioritization = PrioritizationConfig(
    cadd_scores_path=Path("cadd.tsv"),
    revel_scores_path=Path("revel.tsv"),
    spliceai_scores_path=Path("spliceai.tsv"),
)

cohort_config = CohortConfig(
    sample_vcfs=[Path(f"sample_{i:03d}.vcf.gz") for i in range(1, 21)],
    output_path=Path("cohort_results/"),
    cohort_name="study_arm_a",
    min_recurrence=3,
    parallel=True,
    max_workers=8,
)

pipeline = CohortPipeline(
    cohort_config=cohort_config,
    annotation_config=annotation,
    prioritization_config=prioritization,
)
pipeline.run()
```

### Accessing results programmatically

After `run()` completes, the pipeline exposes structured results:

```python
pipeline.run()

# All aggregated variants
for v in pipeline.variants:
    print(f"{v.gene_name} {v.chrom}:{v.pos} {v.ref}>{v.alt} "
          f"in {v.sample_count}/{v.total_samples} samples "
          f"(freq={v.cohort_frequency:.2f})")

# Per-gene burden
for burden in pipeline.gene_burdens:
    print(f"{burden.gene_name}: "
          f"{burden.total_variants} variants, "
          f"{burden.pathogenic_count} pathogenic, "
          f"penetrance={burden.penetrance:.0%}")

# Summary metrics
summary = pipeline.summary
print(f"Total variants: {summary.total_variants}")
print(f"Shared (>=2 samples): {summary.shared_variants}")
print(f"Singletons: {summary.singleton_variants}")
print(f"Pathogenic: {summary.pathogenic_variants}")
print(f"Top genes: {', '.join(summary.top_recurrent_genes)}")
```

### Using the aggregator directly

For custom workflows, use `CohortAggregator` independently of `CohortPipeline`:

```python
from pathlib import Path
from vartriage import CohortAggregator, CohortConfig
from vartriage.models.variant import ClassifiedVariant

cohort_config = CohortConfig(
    sample_vcfs=[Path("a.vcf.gz"), Path("b.vcf.gz")],
    output_path=Path("out/"),
)

aggregator = CohortAggregator(cohort_config)

# Feed in pre-classified variants from your own pipeline runs
aggregator.add_sample("sample_a", Path("a.vcf.gz"), classified_a)
aggregator.add_sample("sample_b", Path("b.vcf.gz"), classified_b)

# Get merged results
cohort_variants = aggregator.aggregate()
recurrent_only = aggregator.get_recurrent_variants(min_count=2)
```

## Output files

Each cohort run produces three files in the output directory:

### Variants file (`{cohort_name}_variants.json`)

Array of objects, one per distinct variant coordinate:

```json
[
  {
    "chrom": "chr17",
    "pos": 7578406,
    "ref": "C",
    "alt": "T",
    "gene_name": "TP53",
    "consequence": "Missense",
    "sample_count": 4,
    "total_samples": 10,
    "cohort_frequency": 0.4,
    "is_singleton": false,
    "is_universal": false,
    "allele_frequency": 0.00002,
    "max_classification": "Likely_Pathogenic",
    "evidence_tags": ["PM2", "PP3", "PP5"],
    "samples": [
      {"sample_id": "patient_001", "classification": "Likely_Pathogenic", "evidence_tags": ["PM2", "PP3"]},
      {"sample_id": "patient_003", "classification": "Likely_Pathogenic", "evidence_tags": ["PM2", "PP3", "PP5"]},
      {"sample_id": "patient_007", "classification": "VUS", "evidence_tags": ["PM2"]},
      {"sample_id": "patient_009", "classification": "Likely_Pathogenic", "evidence_tags": ["PM2", "PP3"]}
    ]
  }
]
```

### Gene burden file (`{cohort_name}_gene_burden.json`)

```json
[
  {
    "gene_name": "TP53",
    "total_variants": 5,
    "pathogenic_count": 3,
    "samples_affected": 7,
    "total_samples": 10,
    "penetrance": 0.7,
    "most_severe": "Pathogenic"
  }
]
```

### Summary file (`{cohort_name}_summary.json`)

```json
{
  "cohort_name": "cardiac_cohort",
  "total_samples": 10,
  "total_variants": 1847,
  "shared_variants": 312,
  "singleton_variants": 1535,
  "universal_variants": 2,
  "pathogenic_variants": 14,
  "likely_pathogenic_variants": 38,
  "genes_affected": 423,
  "top_recurrent_genes": ["TTN", "MYBPC3", "MYH7", "KCNQ1", "SCN5A"],
  "samples_processed": ["patient_001", "patient_002", "patient_003"]
}
```

## CohortConfig reference

| Field | Type | Default | Validation |
| --- | --- | --- | --- |
| sample_vcfs | list[Path] | required | At least 2 paths |
| output_path | Path | required | Directory (created if missing) |
| cohort_name | str | "cohort" | Used in output filenames |
| min_recurrence | int | 2 | >= 1 |
| output_format | str | "json" | "json" or "csv" |
| max_af_threshold | float | 0.05 | 0.0 to 1.0 |
| include_singletons | bool | True | Include single-sample variants |
| sample_labels | dict | None | Map file stems to display names |
| parallel | bool | False | Enable thread-pool processing |
| max_workers | int | 4 | >= 1 |

## How aggregation works

1. Each sample VCF runs through the standard pipeline (parse, filter, annotate, score, classify).
2. Classified variants from all samples are grouped by coordinate key `(chrom, pos, ref, alt)`.
3. Variants with population AF above `max_af_threshold` are excluded before grouping.
4. For each coordinate group, the aggregator picks:
   - The most severe ACMG classification across samples
   - The most severe functional consequence across samples
   - The union of all evidence tags
   - The first non-null gene name and allele frequency
5. The resulting `CohortVariant` records include per-sample detail so you can trace which samples contributed what evidence.

## Performance considerations

- Memory scales with the number of variants that pass filtering across all samples. For a 20-sample cohort where each sample yields 5,000 classified variants, expect roughly 100K entries in the aggregator (many will merge to the same coordinate).
- Parallel mode helps when per-sample processing is I/O-bound (VCF parsing, reference file reads). For CPU-bound annotation with large GTFs, gains depend on GIL release patterns in pysam.
- The cohort aggregation step itself is fast (dict-based grouping). The bottleneck is per-sample pipeline execution.
- For cohorts above 50 samples, consider pre-running each sample with the standard pipeline and feeding results to `CohortAggregator` directly.
