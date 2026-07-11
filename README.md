# vartriage

Variant prioritization pipeline for whole-genome sequencing data. Takes a VCF, applies quality filters, annotates functional consequence and population frequency, scores pathogenicity via CADD/REVEL/SpliceAI, runs ACMG/AMP evidence classification, and outputs a ranked candidate list.

**Benchmarks:**

| Workload | Variants | Wall time | Peak RSS | Throughput |
|---|---|---|---|---|
| GIAB HG002 (QC only, no annotation) | 4,048,342 | 156 s | 122 MB | ~26K var/sec |
| chr22 full annotation (GENCODE + 4.8M gnomAD) | 130,141 | 36.3 s | ~2 GB | ~3.6K var/sec |
| chr22 annotation (100K gnomAD subset) | 130,141 | 19.5 s | 453 MB | ~6.7K var/sec |

Streaming architecture, so JSON and CSV reports never buffer the full variant set in memory. Reference files (GTF, CADD, REVEL, SpliceAI) are cached after first parse. Subsequent runs load from cache in seconds.

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
  --revel-scores revel_scores.tsv \
  --spliceai-scores spliceai_scores.tsv
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

**Quality filtering** - Drops variants where FILTER isn't PASS/`.`, QUAL is below threshold (default 20), or QUAL is missing entirely.

**Annotation** - Adds functional consequence (from GTF gene models), population frequency (gnomAD), and ClinVar significance. Multiple-transcript conflicts resolve to the most damaging consequence. Consequence severity: Frameshift > Nonsense > Splice_Site > Missense > In_Frame_Insertion > In_Frame_Deletion > Synonymous > Intergenic.

**Prioritization** - Two phases. First: frequency gate drops variants with AF above the threshold (default 0.01); unknown-frequency variants always pass. Second: composite scoring from normalized CADD Phred, REVEL, and SpliceAI:

```
composite = (REVEL × 0.5) + (CADD_normalized × 0.3) + (SpliceAI × 0.2)
```

When only two scores are available, weights redistribute proportionally. Single available score is used directly. Falls back to the legacy two-score formula (0.6/0.4) when SpliceAI is not configured.

**ACMG classification** - Tags evidence per ACMG/AMP 2015 guidelines:

| Tag | Condition |
|------|----------------------------------------------|
| PVS1 | Nonsense, Frameshift, or Splice_Site + SpliceAI > 0.8 |
| PM2 | gnomAD AF < 0.0001 |
| PP3 | REVEL > 0.7 or SpliceAI > 0.5 on splice-adjacent |
| PP5 | ClinVar Pathogenic without conflicting Benign |

Tags combine into Pathogenic, Likely_Pathogenic, or VUS. Missing data sources mean the tag is simply omitted.

**Report output** - JSON and CSV stream directly from the iterator (no buffering). PDF materializes for page layout. Output fields: chromosome, position, ref/alt alleles, functional consequence, allele frequency, composite rank, ClinVar assertion, ACMG classification, evidence tags.

## Configuration

### QualityFilterConfig

| Field | Type | Default | Range |
|---|---|---|---|
| min_qual | float | 20.0 | 0-1,000,000 |

### AnnotationConfig

| Field | Type | Default | Notes |
|---|---|---|---|
| gene_annotation_path | Path | required | GTF/GFF |
| gnomad_path | Path | required | TSV or tabix VCF (.vcf.bgz/.vcf.gz) |
| clinvar_path | Path | None | TSV |
| batch_size | int | 10,000 | 1,000-100,000 |

### PrioritizationConfig

| Field | Type | Default | Notes |
|---|---|---|---|
| max_allele_frequency | float | 0.01 | 0.0-1.0 |
| cadd_scores_path | Path | None | CADD Phred TSV |
| revel_scores_path | Path | None | REVEL TSV |
| spliceai_scores_path | Path | None | SpliceAI TSV |
| batch_size | int | 10,000 | 1,000-100,000 |

### ReportConfig

| Field | Type | Default | Options |
|---|---|---|---|
| output_format | str | "json" | "json", "csv", "pdf" |

## Reference file formats

All TSV with a header row. Tab-separated.

**gnomAD (TSV)** - columns: `chrom`, `pos`, `ref`, `alt`, `af`. The value `'.'` in the af column is treated as null (gnomAD compatibility).

**gnomAD (tabix VCF)** - bgzipped VCF with a `.tbi` index (`.vcf.bgz` or `.vcf.gz`). When you point `gnomad_path` at a tabix-indexed file, vartriage queries it on the fly with zero memory overhead for the reference. Useful when your gnomAD file is too large to fit in RAM as a dict.

**ClinVar** - columns: `chrom`, `pos`, `ref`, `alt`, `clinical_significance`. Values: Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign.

**CADD / REVEL / SpliceAI** - columns: `chrom`, `pos`, `ref`, `alt`, `score`. Lines starting with `#` are skipped. All three use the same TSV format.

## Missing data handling

Variants absent from gnomAD are never dropped; they get `frequency_unknown=True` and pass the frequency filter. Same for ClinVar: no match means `clinvar_unknown=True`.

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
| pysam >=0.22,<1.0 | yes | - | VCF streaming via htslib |
| numpy >=1.24,<3.0 | yes | - | Score normalization |
| polars >=0.20,<2.0 | no | [accelerated] | Batch frequency/ClinVar joins |
| pyranges >=0.1,<1.0 | no | [accelerated] | Interval overlap queries |
| reportlab >=4.0,<5.0 | no | [pdf] | PDF report rendering |

Without optional extras, the library uses pure-Python fallbacks (dict lookups, bisect-based interval tree). Same output either way; the accelerated path is faster on large reference files.

## Caching

Reference files (GTF gene models, CADD scores, REVEL scores, SpliceAI scores) are parsed once and cached as pickle files adjacent to the source (with a `.vartriage.cache` suffix). On subsequent runs, the cache loads in seconds instead of re-parsing.

Cache invalidation is automatic: if the source file's mtime changes or the vartriage version changes, the cache rebuilds. Writes are atomic (temp file + rename), so a crash mid-write won't corrupt anything.

To force a fresh parse, delete the `.vartriage.cache` file next to your reference.

## Type checking

The package ships a `py.typed` marker (PEP 561). All protocol return types are fully typed, with no `Any` in the annotation engine interfaces.

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
    prioritization/       # AF gating + CADD/REVEL/SpliceAI scoring (ScoreLoader)
    classification/       # ACMG evidence tagging
    reporting/            # JSON, CSV, PDF (streaming writers)
    models/               # Dataclasses, enums, configs, warnings
    _internal/            # Batch utils, interval tree, caching, vectorized ops
    py.typed              # PEP 561 marker
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT
