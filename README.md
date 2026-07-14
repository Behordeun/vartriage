# vartriage

Clinical variant triage for gene panels. VCF in, ACMG-classified report out.

```bash
pip install vartriage[all]
vartriage --vcf patient.vcf.gz --output report.html --output-format clinical-html \
  --patient-id PAT-001 --panel-name "Cardiac Panel v3" --use-bundles
```

**What it does:** quality filtering, consequence annotation (GENCODE), population frequency lookup (gnomAD), pathogenicity scoring (CADD/REVEL/SpliceAI), ACMG/AMP classification, trio inheritance analysis, and clinical report generation with audit trail.

**Why use it:**

- Single Python package, no Java/Perl/Spark dependencies
- Streams 4M+ variant WGS files under 2 GB RAM
- Trio-aware: de novo, dominant, recessive, compound het, X-linked
- Score bundle downloader: `vartriage bundle download --bundle clinvar` fetches and prepares reference files
- Outputs: JSON, CSV, PDF, HTML clinical reports, IGV-loadable annotated VCF
- Typed API with Protocol-based backends

**Benchmarks:**

| Workload | Variants | Wall time | Peak RSS |
| --- | --- | --- | --- |
| WGS QC only | 4M | 156 s | 122 MB |
| chr22 full annotation | 130K | 36 s | ~2 GB |
| chr22 annotation (100K gnomAD) | 130K | 19.5 s | 453 MB |

Reference files are cached after first parse. Subsequent runs load from cache in seconds.

## Install

```bash
pip install vartriage
```

Optional extras:

```bash
pip install vartriage[accelerated]   # polars + pyranges backends
pip install vartriage[pdf]           # reportlab PDF reports
pip install vartriage[clinical]      # weasyprint + python-docx for clinical reports
pip install vartriage[api]           # httpx for API annotation mode
pip install vartriage[all]           # everything
```

## CLI

```bash
vartriage --vcf sample.vcf.gz --output candidates.json
```

### Score bundles (new in v0.6.0)

Download reference files automatically:

```bash
# See available bundles
vartriage bundle list

# Download ClinVar + gnomAD for chr22
vartriage bundle download --bundle clinvar
vartriage bundle download --bundle gnomad-exomes-chr22

# Run with auto-resolved reference paths
vartriage --vcf sample.vcf.gz --output results.json --use-bundles
```

### API mode (new in v0.7.0)

Annotate variants via remote APIs with zero local reference files:

```bash
# Gene panel, no downloads needed
vartriage --vcf panel.vcf --output results.json --mode api

# Hybrid: local gnomAD + API for ClinVar/CADD
vartriage --vcf panel.vcf --output results.json --mode hybrid --gnomad gnomad.tsv

# With NCBI API key for faster ClinVar queries
vartriage --vcf panel.vcf --output results.json --mode api --api-key YOUR_KEY
```

Queries Ensembl VEP, ClinVar, CADD, and SpliceAI. Responses are cached in SQLite for instant re-runs. See [docs/api-mode.md](docs/api-mode.md) for configuration and performance details.

### Full options

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
  --spliceai-scores spliceai_scores.tsv \
  --gene-list my_panel.txt \
  --regions target_regions.bed \
  --sample PROBAND_01 \
  --min-gq 20
```

Clinical report options:

```bash
vartriage \
  --vcf sample.vcf.gz \
  --output clinical_report.html \
  --output-format clinical-html \
  --patient-id PAT-2026-001 \
  --panel-name "Cardiac Panel v3" \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar_20240101.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv \
  --spliceai-scores spliceai_scores.tsv \
  --gene-list cardiac_panel.txt
```

Formats: `clinical-html` (self-contained HTML), `clinical-pdf` (requires weasyprint), `clinical-docx` (requires python-docx). Both `--patient-id` and `--panel-name` are required for clinical formats.

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

```text
VCFParser → [SampleExtractor] → [RegionFilter] → QualityFilter → AnnotationEngine → [GeneFilter] → PrioritizationEngine → ACMGClassifier → ReportGenerator
```

Stages in brackets are optional and activate based on config.

**Sample extraction** (`--sample`) - Pulls a single sample from multi-sample VCFs. Only variants where the named sample carries an alternate allele are kept. Optional `--min-gq` threshold drops low-confidence genotype calls.

**Region filtering** (`--regions`) - Restricts to variants overlapping intervals in a BED file. Useful for gene panel target regions.

**Quality filtering** - Drops variants where FILTER isn't PASS/`.`, QUAL is below threshold (default 20), or QUAL is missing entirely.

**Annotation** - Adds functional consequence (from GTF gene models), population frequency (gnomAD), and ClinVar significance. Multiple-transcript conflicts resolve to the most damaging consequence. Consequence severity: Frameshift > Nonsense > Splice_Site > Missense > In_Frame_Insertion > In_Frame_Deletion > Synonymous > Intergenic.

**Gene filtering** (`--gene-list`) - After annotation, restricts to variants in genes from a user-supplied text file. Case-insensitive matching. Logs a warning for any panel genes with zero hits (catches typos).

**Prioritization** - Two phases. First: frequency gate drops variants with AF above the threshold (default 0.01); unknown-frequency variants always pass. Second: composite scoring from normalized CADD Phred, REVEL, and SpliceAI:

```text
composite = (REVEL × 0.5) + (CADD_normalized × 0.3) + (SpliceAI × 0.2)
```

When only two scores are available, weights redistribute proportionally. Single available score is used directly. Falls back to the legacy two-score formula (0.6/0.4) when SpliceAI is not configured.

**ACMG classification** - Tags evidence per ACMG/AMP 2015 guidelines:

| Tag | Condition |
| ------ | ---------------------------------------------- |
| PVS1 | Nonsense, Frameshift, or Splice_Site + SpliceAI > 0.8 |
| PM2 | gnomAD AF < 0.0001 |
| PP3 | REVEL > 0.7 or SpliceAI > 0.5 on splice-adjacent |
| PP5 | ClinVar Pathogenic without conflicting Benign |

Tags combine into Pathogenic, Likely_Pathogenic, or VUS. Missing data sources mean the tag is simply omitted.

**Report output** - JSON and CSV stream directly from the iterator (no buffering). PDF materializes for page layout. VCF re-reads the source file, injects VARTRIAGE_* INFO fields for classified variants, and writes bgzipped output with a tabix index. Clinical formats (`clinical-html`, `clinical-pdf`, `clinical-docx`) produce structured reports with per-variant evidence narratives, an executive summary, findings table, evidence cards, limitations, methodology, and sign-off sections. A JSON audit trail sidecar (`.audit.json`) is written alongside each clinical report. Output fields: chromosome, position, ref/alt alleles, functional consequence, allele frequency, composite rank, ClinVar assertion, ACMG classification, evidence tags.

## Configuration

### QualityFilterConfig

| Field | Type | Default | Range |
| --- | --- | --- | --- |
| min_qual | float | 20.0 | 0-1,000,000 |

### AnnotationConfig

| Field | Type | Default | Notes |
| --- | --- | --- | --- |
| gene_annotation_path | Path | required | GTF/GFF |
| gnomad_path | Path | required | TSV or tabix VCF (.vcf.bgz/.vcf.gz) |
| clinvar_path | Path | None | TSV |
| batch_size | int | 10,000 | 1,000-100,000 |

### PrioritizationConfig

| Field | Type | Default | Notes |
| --- | --- | --- | --- |
| max_allele_frequency | float | 0.01 | 0.0-1.0 |
| cadd_scores_path | Path | None | CADD Phred TSV |
| revel_scores_path | Path | None | REVEL TSV |
| spliceai_scores_path | Path | None | SpliceAI TSV |
| batch_size | int | 10,000 | 1,000-100,000 |

### ReportConfig

| Field | Type | Default | Options |
| --- | --- | --- | --- |
| output_format | str | "json" | "json", "csv", "pdf", "vcf", "clinical-html", "clinical-pdf", "clinical-docx" |

### ClinicalReportConfig

| Field | Type | Default | Options |
| --- | --- | --- | --- |
| patient_id | str | required | Patient identifier |
| panel_name | str | required | Gene panel name |
| output_format | str | required | "clinical-pdf", "clinical-html", "clinical-docx" |
| report_template | str | "standard" | Template name |

Constructed automatically when `--output-format` is a `clinical-*` value. Requires `--patient-id` and `--panel-name`.

### GeneFilterConfig

| Field | Type | Default | Notes |
| --- | --- | --- | --- |
| gene_list_path | Path | required | Plain text, one gene symbol per line |

### RegionFilterConfig

| Field | Type | Default | Notes |
| --- | --- | --- | --- |
| bed_path | Path | required | BED file with target intervals |

### SampleConfig

| Field | Type | Default | Notes |
| --- | --- | --- | --- |
| sample_name | str | required | Sample name from VCF header |
| min_gq | int | None | Genotype quality threshold (0-99) |

## Reference file formats

All TSV with a header row. Tab-separated.

**Gene list** - Plain text file, one gene symbol per line. Comment lines starting with `#` are ignored. Blank lines are skipped. Symbols are matched case-insensitively.

```text
# Cardiac panel v2
BRCA1
BRCA2
TP53
MLH1
```

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
| --- | --- | --- | --- |
| pysam >=0.22,<1.0 | yes | - | VCF streaming via htslib |
| numpy >=1.24,<3.0 | yes | - | Score normalization |
| polars >=0.20,<2.0 | no | [accelerated] | Batch frequency/ClinVar joins |
| pyranges >=0.1,<1.0 | no | [accelerated] | Interval overlap queries |
| reportlab >=4.0,<5.0 | no | [pdf] | PDF report rendering |
| weasyprint >=60.0,<62.0 | no | [clinical] | Clinical PDF rendering |
| python-docx >=1.0,<2.0 | no | [clinical] | Clinical DOCX rendering |

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

```text
vartriage/
    cli.py                # CLI entry point
    pipeline.py           # Orchestrator
    protocols.py          # Protocol interfaces (IntervalIndex, FrequencyDatabase, etc.)
    io/                   # VCF parsing
    filter/               # Quality, region, sample, and gene filtering
    annotation/           # Consequence, frequency, ClinVar lookups
    prioritization/       # AF gating + CADD/REVEL/SpliceAI scoring (ScoreLoader)
    classification/       # ACMG evidence tagging
    reporting/            # JSON, CSV, PDF, VCF (streaming writers)
        clinical/         # Clinical report generation (HTML/PDF/DOCX + audit trail)
    models/               # Dataclasses, enums, configs, warnings
    _internal/            # Batch utils, interval tree, caching, vectorized ops
    py.typed              # PEP 561 marker
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT
