# vartriage

Variant prioritization pipeline for whole-genome sequencing data. Takes a VCF, applies quality filters, annotates functional consequence and population frequency, scores pathogenicity via CADD/REVEL, runs ACMG/AMP evidence classification, and outputs a ranked candidate list.

**Benchmarks (GIAB HG002, 4,048,342 variants):**

- Peak RSS: 122 MB
- Wall time: 156 s (~26K variants/sec)
- Ti/Tv ratio: 2.10

Streaming architecture — JSON and CSV reports never buffer the full variant set in memory.

## Install

```bash
pip install vartriage
```

Optional extras:

```bash
pip install vartriage[accelerated]   # polars + pyranges backends
pip install vartriage[pdf]           # reportlab PDF reports
pip install vartriage[all]           # everything
```

## CLI

```bash
vartriage --vcf sample.vcf.gz --output candidates.json
```

Full options:

```bash
vartriage \
  --vcf sample.vcf.gz \
  --output report.json \
  --output-format json \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar_20240101.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv
```

Run `vartriage --help` for the complete list.

## Python API

Run the whole pipeline:

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

Or use stages individually:

```python
from vartriage import VCFParser, QualityFilter, QualityFilterConfig

with VCFParser(Path("input.vcf.gz")) as parser:
    qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
    for variant in qf.apply(iter(parser)):
        print(f"{variant.chrom}:{variant.pos} {variant.ref}>{variant.alt}")
```

## Pipeline stages

```
VCFParser → QualityFilter → AnnotationEngine → PrioritizationEngine → ACMGClassifier → ReportGenerator
```

**Quality filtering** — Drops variants where FILTER isn't PASS/`.`, QUAL is below threshold (default 20), or QUAL is missing entirely.

**Annotation** — Adds functional consequence (from GTF gene models), population frequency (gnomAD), and ClinVar significance. Multiple-transcript conflicts resolve to the most damaging consequence. Consequence severity: Frameshift > Nonsense > Splice_Site > Missense > In_Frame_Insertion > In_Frame_Deletion > Synonymous > Intergenic.

**Prioritization** — Two phases. First: frequency gate drops variants with AF above the threshold (default 0.01); unknown-frequency variants always pass. Second: composite scoring from normalized CADD Phred and REVEL:

```
composite = (REVEL × 0.6) + (CADD_normalized × 0.4)
```

Falls back to the single available score when only one source exists.

**ACMG classification** — Tags evidence per ACMG/AMP 2015 guidelines:

| Tag | Condition |
|------|----------------------------------------------|
| PVS1 | Nonsense or Frameshift |
| PM2 | gnomAD AF < 0.0001 |
| PP3 | REVEL > 0.7 |
| PP5 | ClinVar Pathogenic without conflicting Benign |

Tags combine into Pathogenic, Likely_Pathogenic, or VUS. Missing data sources mean the tag is simply omitted.

**Report output** — JSON and CSV stream directly from the iterator (no buffering). PDF materializes for page layout. Output fields: chromosome, position, ref/alt alleles, functional consequence, allele frequency, composite rank, ClinVar assertion, ACMG classification, evidence tags.

## Configuration

### QualityFilterConfig

| Field | Type | Default | Range |
|---|---|---|---|
| min_qual | float | 20.0 | 0–1,000,000 |

### AnnotationConfig

| Field | Type | Default | Notes |
|---|---|---|---|
| gene_annotation_path | Path | required | GTF/GFF |
| gnomad_path | Path | required | TSV |
| clinvar_path | Path | None | TSV |
| batch_size | int | 10,000 | 1,000–100,000 |

### PrioritizationConfig

| Field | Type | Default | Notes |
|---|---|---|---|
| max_allele_frequency | float | 0.01 | 0.0–1.0 |
| cadd_scores_path | Path | None | CADD Phred TSV |
| revel_scores_path | Path | None | REVEL TSV |
| batch_size | int | 10,000 | 1,000–100,000 |

### ReportConfig

| Field | Type | Default | Options |
|---|---|---|---|
| output_format | str | "json" | "json", "csv", "pdf" |

## Reference file formats

All TSV with a header row. Tab-separated.

**gnomAD** — columns: `chrom`, `pos`, `ref`, `alt`, `af`. The value `'.'` in the af column is treated as null (gnomAD compatibility).

**ClinVar** — columns: `chrom`, `pos`, `ref`, `alt`, `clinical_significance`. Values: Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign.

**CADD / REVEL** — columns: `chrom`, `pos`, `ref`, `alt`, `score`. Lines starting with `#` are skipped.

## Missing data handling

Variants absent from gnomAD are never dropped — they get `frequency_unknown=True` and pass the frequency filter. Same for ClinVar: no match means `clinvar_unknown=True`.

A `MissingDataWarning` fires per lookup miss. After a run:

```python
acc = pipeline.warning_accumulator
print(f"{acc.total_count} missing data events across {acc.sources}")
```

## Warning hierarchy

All warnings inherit from `VarTriageWarning` (a `UserWarning` subclass). Silence everything at once:

```python
import warnings
from vartriage import VarTriageWarning
warnings.filterwarnings("ignore", category=VarTriageWarning)
```

## Dependencies

| Package | Required | Extra | Purpose |
|---|---|---|---|
| pysam >=0.22,<1.0 | yes | — | VCF streaming via htslib |
| numpy >=1.24,<3.0 | yes | — | Score normalization |
| polars >=0.20,<2.0 | no | [accelerated] | Batch frequency/ClinVar joins |
| pyranges >=0.1,<1.0 | no | [accelerated] | Interval overlap queries |
| reportlab >=4.0,<5.0 | no | [pdf] | PDF report rendering |

Without optional extras, the library uses pure-Python fallbacks (dict lookups, bisect-based interval tree). Same output either way; the accelerated path is faster on large reference files.

## Type checking

The package ships a `py.typed` marker (PEP 561). All protocol return types are fully typed — no `Any` in the annotation engine interfaces.

```bash
mypy --strict vartriage/
```

## Tests

```bash
pytest tests/                     # full suite
pytest tests/ -m "not slow"       # skip benchmarks
```

## CI

GitHub Actions runs on Python 3.10, 3.11, and 3.12. PyPI publishing uses trusted publisher (no token in secrets).

## Project layout

```
vartriage/
    cli.py                # CLI entry point
    pipeline.py           # Orchestrator
    protocols.py          # Protocol interfaces (IntervalIndex, FrequencyDatabase, etc.)
    io/                   # VCF parsing
    filter/               # Quality-based exclusion
    annotation/           # Consequence, frequency, ClinVar lookups
    prioritization/       # AF gating + CADD/REVEL scoring (ScoreLoader)
    classification/       # ACMG evidence tagging
    reporting/            # JSON, CSV, PDF — streaming writers
    models/               # Dataclasses, enums, configs, warnings
    _internal/            # Batch utils, interval tree, vectorized ops
    py.typed              # PEP 561 marker
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT
