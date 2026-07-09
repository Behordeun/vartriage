# vartriage

Variant prioritization pipeline for whole-genome sequencing data. Reads a VCF, applies quality filters, annotates functional consequence and population frequency, computes pathogenicity scores, runs ACMG/AMP evidence classification, and writes a ranked candidate list in JSON, CSV, or PDF.

Processes 4M+ variant WGS files under 2GB memory via batched iterators.

## Install

```bash
pip install vartriage
```

With faster annotation backends (polars + pyranges):

```bash
pip install vartriage[accelerated]
```

With PDF report support:

```bash
pip install vartriage[pdf]
```

All optional extras:

```bash
pip install vartriage[all]
```

## Usage

```python
from pathlib import Path
from vartriage import (
    Pipeline, PipelineConfig, AnnotationConfig,
    PrioritizationConfig, QualityFilterConfig, ReportConfig,
)

config = PipelineConfig(
    vcf_path=Path("sample.vcf.gz"),
    output_path=Path("candidates.json"),
    quality_filter=QualityFilterConfig(min_qual=30.0),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
        clinvar_path=Path("clinvar_20240101.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.01,
        cadd_scores_path=Path("cadd_scores.tsv"),
        revel_scores_path=Path("revel_scores.tsv"),
    ),
    report=ReportConfig(output_format="json"),
)

pipeline = Pipeline(config)
pipeline.run()
```

Individual stages work on their own:

```python
from vartriage import VCFParser, QualityFilter, QualityFilterConfig

with VCFParser(Path("input.vcf.gz")) as parser:
    qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
    for variant in qf.apply(iter(parser)):
        print(f"{variant.chrom}:{variant.pos} {variant.ref}>{variant.alt}")
```

## Pipeline stages

```text
VCFParser > QualityFilter > AnnotationEngine > PrioritizationEngine > ACMGClassifier > ReportGenerator
```

### Quality filtering

Drops variants where:

- `FILTER` is not `PASS` or `.`
- `QUAL` is below the threshold (default 20)
- `QUAL` field is missing (emits a warning)

Passing variants keep their original order.

### Annotation

Adds three annotations to each surviving variant:

**Functional consequence:** Looked up against gene models (GTF/GFF). Splice_Site applies within 2bp of an exon-intron boundary. When multiple transcripts disagree, the most damaging consequence wins. Severity ranking (highest first): Frameshift, Nonsense, Splice_Site, Missense, In_Frame_Insertion, In_Frame_Deletion, Synonymous, Intergenic.

**Population frequency:** Matched against gnomAD by (chrom, pos, ref, alt). Variants not found get `frequency_unknown=True` and a `MissingDataWarning`.

**ClinVar assertion:** Pathogenic, Likely_Pathogenic, VUS, Likely_Benign, or Benign when available.

### Prioritization

Two phases:

1. Frequency gate: drops variants with AF above the threshold (default 0.01). Variants marked `frequency_unknown` always pass.
2. Composite scoring: normalizes CADD Phred (divide by 99, cap at 1.0) and REVEL (already 0-1), then computes:

```text
composite = (REVEL x 0.6) + (CADD_normalized x 0.4)
```

Falls back to the single available score when only one source exists. Output sorted descending by composite rank; variants without scores go last.

### ACMG classification

Evidence tagging per ACMG/AMP 2015:

| Tag  | Strength    | Condition                                 |
| ---- | ----------- | ----------------------------------------- |
| PVS1 | Very Strong | Nonsense or Frameshift                    |
| PM2  | Moderate    | gnomAD AF < 0.0001                        |
| PP3  | Supporting  | REVEL > 0.7                               |
| PP5  | Supporting  | ClinVar Pathogenic, no conflicting Benign |

Tags combine per standard rules into: Pathogenic, Likely_Pathogenic, or VUS. If a data source is unavailable, the corresponding tag is omitted.

### Report output

Fields in all formats:

| Field                      | Description                 |
| -------------------------- | --------------------------- |
| `chromosome`             | Chromosome name             |
| `position`               | 1-based position            |
| `ref_allele`             | Reference allele            |
| `alt_allele`             | Alternate allele            |
| `functional_consequence` | Most severe consequence     |
| `allele_frequency`       | gnomAD AF (null if unknown) |
| `composite_rank`         | Pathogenicity score 0-1     |
| `clinvar_assertion`      | ClinVar significance        |
| `acmg_classification`    | Final classification        |
| `evidence_tags`          | Applied evidence codes      |

Null values: `null` in JSON, empty in CSV, `N/A` in PDF.

## Configuration

### QualityFilterConfig

| Field        | Type  | Default | Range          |
| ------------ | ----- | ------- | -------------- |
| `min_qual` | float | 20.0    | 0 to 1,000,000 |

### AnnotationConfig

| Field                    | Type | Default  | Notes                  |
| ------------------------ | ---- | -------- | ---------------------- |
| `gene_annotation_path` | Path | required | GTF/GFF                |
| `gnomad_path`          | Path | required | TSV (see format below) |
| `clinvar_path`         | Path | None     | TSV (see format below) |
| `batch_size`           | int  | 10,000   | 1,000 to 100,000       |

### PrioritizationConfig

| Field                    | Type  | Default | Range            |
| ------------------------ | ----- | ------- | ---------------- |
| `max_allele_frequency` | float | 0.01    | 0.0 to 1.0       |
| `cadd_scores_path`     | Path  | None    | CADD Phred TSV   |
| `revel_scores_path`    | Path  | None    | REVEL scores TSV |
| `batch_size`           | int   | 10,000  | 1,000 to 100,000 |

### ReportConfig

| Field             | Type | Default    | Options                          |
| ----------------- | ---- | ---------- | -------------------------------- |
| `output_format` | str  | `"json"` | `"json"`, `"csv"`, `"pdf"` |

### MissingDataConfig

| Field                 | Type | Default | Notes                         |
| --------------------- | ---- | ------- | ----------------------------- |
| `warning_threshold` | int  | 1000    | Summary warning when exceeded |

## Reference file formats

All reference files are tab-separated with a header row.

**gnomAD:**

```tsv
chrom	pos	ref	alt	af
chr1	12345	A	G	0.00032
```

**ClinVar:**

```tsv
chrom	pos	ref	alt	clinical_significance
chr1	12345	A	G	Pathogenic
```

Recognized values: `Pathogenic`, `Likely pathogenic`, `Uncertain significance`, `Likely benign`, `Benign`.

**CADD / REVEL:**

```tsv
chrom	pos	ref	alt	score
chr1	12345	A	G	28.5
```

## Missing data handling

Variants absent from gnomAD are never dropped. They get `frequency_unknown=True` and pass the frequency filter. Same for ClinVar: no match means `clinvar_unknown=True`.

A `MissingDataWarning` is emitted per lookup miss. Once the total exceeds `warning_threshold`, a summary fires with the count and contributing sources.

```python
pipeline.run()
acc = pipeline.warning_accumulator
print(f"{acc.total_count} missing data events across {acc.sources}")
```

## Dependencies

| Package   | Required | Extra             | Purpose                       |
| --------- | -------- | ----------------- | ----------------------------- |
| pysam     | yes      | n/a               | VCF streaming (htslib)        |
| numpy     | yes      | n/a               | Score normalization           |
| polars    | no       | `[accelerated]` | Batch frequency/ClinVar joins |
| pyranges  | no       | `[accelerated]` | Interval overlap queries      |
| reportlab | no       | `[pdf]`         | PDF report generation         |

Without optional extras, the library uses pure-Python fallbacks (dict-based lookups, bisect-based interval tree). Correct output either way; the accelerated path runs faster on large reference files.

## Error handling

Invalid configuration raises `ValueError` or `FileNotFoundError` at construction time, before any variants are processed.

During processing, missing reference data does not crash. The library assigns null values, sets flags, and continues. After a run, inspect `pipeline.warning_accumulator` to see how many lookup misses occurred and which sources were affected.

## Tests

```bash
pytest tests/                        # full suite
pytest tests/ -m "not slow"          # skip performance benchmarks
mypy --strict vartriage/  # type checking
```

383 tests, 0 failures. mypy strict, 0 errors.

## Project layout

```text
vartriage/
    pipeline.py           # Top-level orchestrator
    protocols.py          # Protocol interfaces for swappable backends
    io/                   # VCF parsing, exceptions
    filter/               # Quality-based exclusion
    annotation/           # Consequence, frequency, ClinVar lookups
    prioritization/       # AF gating + pathogenicity scoring
    classification/       # ACMG evidence tagging + combining
    reporting/            # JSON, CSV, PDF output
    models/               # Dataclasses, enums, configs, warnings
    _internal/            # Batch utils, interval tree, vectorized ops
```

## Requirements

- Python >= 3.10
- pysam >= 0.22.0
- numpy >= 1.24.0

## License

MIT
