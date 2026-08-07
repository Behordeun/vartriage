# vartriage

[![CI](https://github.com/Behordeun/vartriage/actions/workflows/ci.yml/badge.svg)](https://github.com/Behordeun/vartriage/actions/workflows/ci.yml) [![Publish to PyPI](https://github.com/Behordeun/vartriage/actions/workflows/publish.yml/badge.svg)](https://github.com/Behordeun/vartriage/actions/workflows/publish.yml) [![CodeQL](https://github.com/Behordeun/vartriage/actions/workflows/codeql.yml/badge.svg)](https://github.com/Behordeun/vartriage/actions/workflows/codeql.yml)

Clinical variant interpretation library for gene panels and whole genomes. VCF in, ACMG-classified report out.

```bash
pip install vartriage[all]
vartriage --vcf patient.vcf.gz --output report.html --output-format clinical-html \
  --patient-id PAT-001 --panel-name "Cardiac Panel v3" --use-bundles
```

**What it does:** quality filtering, consequence annotation (GENCODE, with codon-level resolution via reference FASTA), population frequency lookup (gnomAD, population-specific via local files or API), pathogenicity scoring (CADD/REVEL/SpliceAI with ClinGen-calibrated thresholds), gene-disease linkage (OMIM/ClinGen/HPO/gnomAD constraint), phenotype-driven prioritization, ACMG/AMP classification (10 criteria: PVS1, PS1, PM2, PM5, PP3, PP5, BA1, BS1, BP4, BP7 with strength modulation), trio inheritance analysis, multi-sample cohort analysis (recurrence, gene burden), ACMG Secondary Findings screening, **structural variant triage (ClinGen 2020 framework)**, **mitochondrial variant analysis (mtDNA-specific classification with heteroplasmy, MITOMAP, and HelixMTdb)**, and clinical report generation with audit trail and computational-only disclaimer.

**Why use it:**

- Single Python package, no Java/Perl/Spark dependencies
- Streams 4M+ variant WGS files under 2 GB RAM
- Codon-level consequence calling with reference FASTA (correct missense vs synonymous)
- Benign + pathogenic ACMG criteria (10 criteria, ClinGen-calibrated): classifies variants across all 5 tiers
- Gene-disease linkage: OMIM, ClinGen validity, HPO phenotype matching, gnomAD constraint, actionability
- Phenotype-driven: `--hpo-terms` boosts variants in genes matching patient symptoms
- Trio-aware: de novo, dominant, recessive, compound het, X-linked
- ACMG Secondary Findings (SF v3.2): screens 71 medically actionable genes
- Score bundle downloader: `vartriage bundle download --bundle clinvar` fetches and prepares reference files
- API mode: annotate gene panels via Ensembl VEP + ClinVar + gnomAD API with zero local files
- Outputs: JSON, CSV, PDF, HTML clinical reports, IGV-loadable annotated VCF
- Structural variant triage: ClinGen 2020 framework for DEL/DUP/INV/INS/BND/CNV
- Mitochondrial variant analysis: automatic chrM detection, heteroplasmy quantification, MITOMAP/HelixMTdb annotation, mtDNA-specific classification
- Typed API with Protocol-based backends

**Benchmarks:**

| Workload                       | Variants | Wall time | Peak RSS |
| ------------------------------ | -------- | --------- | -------- |
| WGS QC only                    | 4M       | 156 s     | 122 MB   |
| chr22 full annotation          | 130K     | 36 s      | ~2 GB    |
| chr22 annotation (100K gnomAD) | 130K     | 19.5 s    | 453 MB   |

Reference files are cached after first parse. Subsequent runs load from cache in seconds.

## Install

```bash
pip install vartriage
```

Optional extras:

```bash
pip install vartriage[accelerated]   # polars + pyranges backends
pip install vartriage[pdf]           # reportlab PDF reports
pip install vartriage[clinical]      # weasyprint + python-docx for clinical HTML/PDF/DOCX reports
pip install vartriage[api]           # httpx for API annotation mode
pip install vartriage[all]           # everything
```

## CLI

```bash
vartriage --vcf sample.vcf.gz --output candidates.json
```

### Score bundles

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

### API mode

Annotate variants via remote APIs with zero local reference files:

```bash
# Gene panel, no downloads needed
vartriage --vcf panel.vcf --output results.json --mode api

# Hybrid: local gnomAD + API for ClinVar/CADD
vartriage --vcf panel.vcf --output results.json --mode hybrid --gnomad gnomad.tsv

# With NCBI API key for faster ClinVar queries
vartriage --vcf panel.vcf --output results.json --mode api --api-key YOUR_KEY
```

Queries Ensembl VEP, ClinVar, CADD, and SpliceAI. Responses are cached in SQLite for instant re-runs. See [API Mode Guide](https://github.com/Behordeun/vartriage/blob/main/docs/api-mode.md) for configuration and performance details.

### Cohort analysis

Analyze multiple samples together to find shared variants, compute recurrence frequencies, and generate per-gene burden reports:

```bash
# Multiple VCF files directly
vartriage cohort --vcf sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
  --output cohort_results/ --cohort-name "cardiac_cohort"

# From a manifest file (one VCF path per line, optional tab-separated labels)
vartriage cohort --manifest samples.tsv --output cohort_results/ \
  --cohort-name "cardiac_cohort" --output-format csv

# With annotation references and gene filtering
vartriage cohort --manifest samples.tsv --output cohort_results/ \
  --gene-annotation gencode.v44.gtf --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv --cadd-scores cadd.tsv --revel-scores revel.tsv \
  --gene-list cardiac_panel.txt --use-bundles

# Parallel processing with custom thresholds
vartriage cohort --manifest samples.tsv --output cohort_results/ \
  --parallel --max-workers 8 --min-recurrence 3 --max-af 0.01 --no-singletons
```

**Manifest format** - plain text, one VCF path per line. Optional tab-separated second column for sample labels:

```text
# Cardiac cohort 2026
/data/vcfs/patient_001.vcf.gz   Patient 001
/data/vcfs/patient_002.vcf.gz   Patient 002
/data/vcfs/patient_003.vcf.gz   Patient 003
```

**Output files** - three files per run:

| File                             | Contents                                                                                  |
| -------------------------------- | ----------------------------------------------------------------------------------------- |
| `{cohort_name}_variants.json`    | All cohort variants with recurrence counts, per-sample classifications, and evidence tags |
| `{cohort_name}_gene_burden.json` | Per-gene statistics: variant count, pathogenic count, penetrance, samples affected        |
| `{cohort_name}_summary.json`     | Top-level metrics: total variants, shared/singleton/universal counts, top recurrent genes |

**Key options:**

| Flag               | Default | Description                                                |
| ------------------ | ------- | ---------------------------------------------------------- |
| `--min-recurrence` | 2       | Exclude variants appearing in fewer than this many samples |
| `--max-af`         | 0.05    | Exclude variants above this population frequency           |
| `--no-singletons`  | false   | Drop variants seen in only one sample                      |
| `--parallel`       | false   | Process samples concurrently                               |
| `--max-workers`    | 4       | Thread pool size for parallel mode                         |

Parallel mode uses `ThreadPoolExecutor`. Per-sample work is I/O-bound (pysam releases the GIL during C-level VCF parsing), so threads are effective here without needing multiprocessing.

### Gene-disease linkage

Connect variants to clinical context with phenotype-driven prioritization:

```bash
# Phenotype-driven: boost epilepsy-related genes for a seizure patient
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --hpo-terms HP:0001250,HP:0001249,HP:0002197

# Filter to autosomal recessive genes
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --inheritance-mode AR

# Only actionable findings (ClinGen curations)
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --flag-actionable

# Combine all three
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --hpo-terms HP:0001250,HP:0001249 --inheritance-mode AD --flag-actionable
```

Output includes per-variant: disease associations (with MIM numbers), ClinGen validity level, gnomAD constraint metrics (pLI/LOEUF/mis_z), actionability status, and phenotype match score.

See [Gene-Disease Linkage Guide](https://github.com/Behordeun/vartriage/blob/main/docs/gene-disease-linkage.md) for data file formats, Python API usage, and validation details.

### Structural variant triage

Classify structural variants (DEL, DUP, INV, INS, BND, CNV) using the ClinGen 2020 technical standards:

```bash
# Standalone SV analysis
vartriage sv --sv-vcf sv_calls.vcf.gz --output sv_report.json \
  --gene-annotation gencode.v46.gtf \
  --dosage-sensitivity clingen_dosage.tsv \
  --gnomad-sv gnomad_sv.bed \
  --pathogenic-regions clingen_pathogenic_regions.bed

# Combined SNV + SV analysis
vartriage --vcf snv.vcf.gz --sv-vcf sv_calls.vcf.gz \
  --output report.json --gene-annotation gencode.v46.gtf
```

Pipeline: SVParser (streams from VCF) -> SVAnnotator (gene overlap, dosage sensitivity, gnomAD-SV frequency) -> SVScorer (composite pathogenicity score) -> SVClassifier (ClinGen evidence sections 1-4, 5-tier classification).

Supports Manta, DELLY, GATK-SV, GRIDSS, and LUMPY. See [Structural Variants Guide](https://github.com/Behordeun/vartriage/blob/main/docs/structural-variants.md) for the full CLI reference and Python API.

### Mitochondrial variant analysis

Automatic detection and classification of mitochondrial DNA variants using mtDNA-specific criteria:

```bash
# Automatic: chrM/MT variants in the VCF are detected and classified separately
vartriage --vcf wgs.vcf.gz --output results.json

# Custom heteroplasmy threshold (default: 1%)
vartriage --vcf wgs.vcf.gz --output results.json --mt-min-heteroplasmy 5.0

# Skip mitochondrial analysis (targeted panels without mtDNA capture)
vartriage --vcf panel.vcf.gz --output results.json --skip-mito
```

The mitochondrial pipeline uses the vertebrate mitochondrial genetic code for amino acid prediction, extracts heteroplasmy levels from AD/AF fields, queries MITOMAP for disease associations, checks HelixMTdb for population frequency, and classifies variants independently of the nuclear ACMG criteria. Results appear in a dedicated "Mitochondrial Findings" section in output reports.

See [Mitochondrial Variants Guide](https://github.com/Behordeun/vartriage/blob/main/docs/mitochondrial.md) for heteroplasmy thresholds, classification rules, and data update instructions.

### Full options

```bash
vartriage \
  --vcf sample.vcf.gz \
  --output report.json \
  --output-format json \
  --reference-fasta GRCh38.fa \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar_20240101.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv \
  --spliceai-scores spliceai_scores.tsv \
  --gene-list my_panel.txt \
  --regions target_regions.bed \
  --sample PROBAND_01 \
  --min-gq 20 \
  --secondary-findings
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

API mode (no local reference files needed):

```python
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, ReportConfig
from vartriage.api.config import APIConfig

api_config = APIConfig.load(
    mode="api",
    genome_build="grch38",
    ncbi_api_key="your-key-here",  # or set NCBI_API_KEY env var
)

config = PipelineConfig(
    vcf_path=Path("panel.vcf"),
    output_path=Path("results.json"),
    report=ReportConfig(output_format="json"),
    api=api_config,
)

pipeline = Pipeline(config)
pipeline.run()
```

Cohort analysis across multiple samples:

```python
from pathlib import Path
from vartriage import CohortPipeline, CohortConfig, PipelineConfig, AnnotationConfig, PrioritizationConfig

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
    output_format="json",
)

# Optional: shared annotation config applied to all samples
annotation = AnnotationConfig(
    gene_annotation_path=Path("gencode.v44.gtf"),
    gnomad_path=Path("gnomad.v4.sites.tsv"),
)

pipeline = CohortPipeline(
    cohort_config=cohort_config,
    annotation_config=annotation,
)
report_paths = pipeline.run()

# Access results programmatically
for variant in pipeline.variants:
    if variant.sample_count >= 2:
        print(f"{variant.gene_name}: {variant.chrom}:{variant.pos} "
              f"in {variant.sample_count}/{variant.total_samples} samples")

for burden in pipeline.gene_burdens:
    if burden.pathogenic_count > 0:
        print(f"{burden.gene_name}: {burden.pathogenic_count} pathogenic, "
              f"penetrance={burden.penetrance:.0%}")
```

Structural variant triage:

```python
from pathlib import Path
from vartriage.structural import SVTriagePipeline, SVTriageConfig

config = SVTriageConfig(
    vcf_path=Path("sv_calls.vcf.gz"),
    output_path=Path("sv_report.json"),
    gene_annotation_path=Path("gencode.v46.gtf"),
    dosage_sensitivity_path=Path("clingen_dosage.tsv"),
    gnomad_sv_path=Path("gnomad_sv.bed"),
    pathogenic_regions_path=Path("pathogenic_regions.bed"),
)

pipeline = SVTriagePipeline(config)
output = pipeline.run()
```

## Pipeline stages

```text
VCFParser → [SampleExtractor] → [RegionFilter] → QualityFilter → AnnotationEngine → [GeneFilter] → [GeneKnowledgeAnnotator] → PrioritizationEngine → [PhenotypeBoost] → ACMGClassifier → ReportGenerator
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

CADD normalization: Phred score divided by 99.0, capped at 1.0. The separate `prioritization_score` field uses Phred / 60.0 (capped at 1.0) for triage ranking. REVEL and SpliceAI are already bounded 0.0–1.0 and used directly without rescaling.

When only two scores are available, weights redistribute proportionally. Single available score is used directly. Falls back to the legacy two-score formula (0.6/0.4) when SpliceAI is not configured.

**ACMG classification** - Tags evidence per ACMG/AMP 2015 guidelines with ClinGen SVI-calibrated thresholds (Pejaver et al. 2022):

| Tag           | Strength    | Condition                                                                                                        |
| ------------- | ----------- | ---------------------------------------------------------------------------------------------------------------- |
| PVS1          | Very Strong | Nonsense, Frameshift, or Splice_Site + SpliceAI > 0.8                                                            |
| PS1           | Strong      | Same amino acid change as ClinVar Pathogenic via different nucleotide (requires protein index + reference FASTA) |
| PM2           | Moderate    | All population AFs < 0.0001 (population-aware)                                                                   |
| PM5           | Moderate    | Novel missense at amino acid position with known pathogenic missense in ClinVar (requires protein index)         |
| PP3           | Supporting  | REVEL > 0.644 or SpliceAI > 0.5 on splice-adjacent                                                               |
| PP3_MODERATE  | Moderate    | REVEL > 0.773 (ClinGen-calibrated)                                                                               |
| PP5           | Supporting  | ClinVar Pathogenic without conflicting Benign                                                                    |
| BA1           | Standalone  | Any population AF > 5% (standalone Benign)                                                                       |
| BS1           | Strong      | Any population AF > 1% (strong benign)                                                                           |
| BP4           | Supporting  | REVEL < 0.290 (missense) or CADD < 10 (non-missense)                                                             |
| BP4_MODERATE  | Moderate    | REVEL < 0.183 (ClinGen-calibrated)                                                                               |
| BP7           | Supporting  | Synonymous + SpliceAI < 0.1                                                                                      |

Tags combine into Pathogenic, Likely_Pathogenic, VUS, Likely_Benign, or Benign. Conflicting pathogenic + benign evidence yields VUS. Missing data sources mean the tag is simply omitted.

**Report output** - JSON and CSV stream directly from the iterator (no buffering). PDF materializes for page layout. VCF re-reads the source file, injects VARTRIAGE_* INFO fields for classified variants, and writes bgzipped output with a tabix index. Clinical formats (`clinical-html`, `clinical-pdf`, `clinical-docx`) produce structured reports with a computational-only disclaimer (citing ACMG/AMP 2015), per-variant evidence narratives, an executive summary, findings table, evidence cards, limitations, methodology, and sign-off sections. A JSON audit trail sidecar (`.audit.json`) is written alongside each clinical report. Output fields: chromosome, position, ref/alt alleles, gene_name, functional consequence, allele frequency, revel_score, composite rank, prioritization_score, ClinVar assertion, ACMG classification, evidence tags, disease_associations, clingen_validity, gene_constraint, is_actionable, phenotype_match_score.

## Configuration

### QualityFilterConfig

| Field    | Type  | Default | Range       |
| -------- | ----- | ------- | ----------- |
| min_qual | float | 20.0    | 0-1,000,000 |

### AnnotationConfig

| Field                | Type | Default  | Notes                                                          |
| -------------------- | ---- | -------- | -------------------------------------------------------------- |
| gene_annotation_path | Path | required | GTF/GFF                                                        |
| gnomad_path          | Path | required | TSV or tabix VCF (.vcf.bgz/.vcf.gz)                            |
| clinvar_path         | Path | None     | TSV                                                            |
| reference_fasta_path | Path | None     | Indexed FASTA (.fa + .fai) for codon-level consequence calling |
| batch_size           | int  | 10,000   | 1,000-100,000                                                  |

### PrioritizationConfig

| Field                | Type  | Default | Notes          |
| -------------------- | ----- | ------- | -------------- |
| max_allele_frequency | float | 0.01    | 0.0-1.0        |
| cadd_scores_path     | Path  | None    | CADD Phred TSV |
| revel_scores_path    | Path  | None    | REVEL TSV      |
| spliceai_scores_path | Path  | None    | SpliceAI TSV   |
| batch_size           | int   | 10,000  | 1,000-100,000  |

### ReportConfig

| Field         | Type | Default | Options                                                                       |
| ------------- | ---- | ------- | ----------------------------------------------------------------------------- |
| output_format | str  | "json"  | "json", "csv", "pdf", "vcf", "clinical-html", "clinical-pdf", "clinical-docx" |

### ClinicalReportConfig

| Field           | Type | Default    | Options                                          |
| --------------- | ---- | ---------- | ------------------------------------------------ |
| patient_id      | str  | required   | Patient identifier                               |
| panel_name      | str  | required   | Gene panel name                                  |
| output_format   | str  | required   | "clinical-pdf", "clinical-html", "clinical-docx" |
| report_template | str  | "standard" | Template name                                    |

Constructed automatically when `--output-format` is a `clinical-*` value. Requires `--patient-id` and `--panel-name`.

### GeneFilterConfig

| Field          | Type | Default  | Notes                                |
| -------------- | ---- | -------- | ------------------------------------ |
| gene_list_path | Path | required | Plain text, one gene symbol per line |

### RegionFilterConfig

| Field    | Type | Default  | Notes                          |
| -------- | ---- | -------- | ------------------------------ |
| bed_path | Path | required | BED file with target intervals |

### SampleConfig

| Field       | Type | Default  | Notes                             |
| ----------- | ---- | -------- | --------------------------------- |
| sample_name | str  | required | Sample name from VCF header       |
| min_gq      | int  | None     | Genotype quality threshold (0-99) |

### CohortConfig

| Field              | Type       | Default  | Notes                                     |
| ------------------ | ---------- | -------- | ----------------------------------------- |
| sample_vcfs        | list[Path] | required | At least 2 VCF file paths                 |
| output_path        | Path       | required | Output directory for reports              |
| cohort_name        | str        | "cohort" | Identifier for output filenames           |
| min_recurrence     | int        | 2        | Minimum samples for recurrence            |
| output_format      | str        | "json"   | "json" or "csv"                           |
| max_af_threshold   | float      | 0.05     | Max population AF for inclusion (0.0-1.0) |
| include_singletons | bool       | True     | Include variants in only 1 sample         |
| sample_labels      | dict       | None     | Map file stems to display labels          |
| parallel           | bool       | False    | Process samples concurrently              |
| max_workers        | int        | 4        | Thread pool size (>= 1)                   |

### KnowledgeBaseConfig

| Field            | Type            | Default  | Notes                                                   |
| ---------------- | --------------- | -------- | ------------------------------------------------------- |
| data_dir         | Path \| None    | None     | Custom TSV directory (defaults to bundled package data) |
| hpo_terms        | frozenset[str]  | empty    | Patient HPO terms (HP:NNNNNNN format)                   |
| inheritance_mode | str \| None     | None     | Filter: AD, AR, XL, XLD, XLR, MT                        |
| flag_actionable  | bool            | False    | Filter to ClinGen actionable genes only                 |

Any non-default field activates the gene-disease linkage pipeline stage.

### MitoConfig

| Field            | Type         | Default | Notes                                                    |
| ---------------- | ------------ | ------- | -------------------------------------------------------- |
| enabled          | bool         | True    | Enable/disable mitochondrial analysis                    |
| min_heteroplasmy | float        | 1.0     | Minimum heteroplasmy % for reporting (0.0-100.0)         |
| gene_map_path    | Path \| None | None    | Custom mt_gene_map.tsv (defaults to bundled)             |
| mitomap_path     | Path \| None | None    | Custom mitomap_pathogenic.tsv (defaults to bundled)      |
| helixmtdb_path   | Path \| None | None    | Custom helixmtdb_frequency.tsv (defaults to bundled)     |

Mitochondrial analysis is auto-enabled when chrM/MT variants are present in the VCF. Use `--skip-mito` (CLI) or `MitoConfig(enabled=False)` to disable.

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

| Package                 | Required | Extra         | Purpose                            |
| ----------------------- | -------- | ------------- | ---------------------------------- |
| pysam >=0.22,<1.0       | yes      | -             | VCF streaming via htslib           |
| numpy >=1.24,<3.0       | yes      | -             | Score normalization                |
| polars >=0.20,<2.0      | no       | [accelerated] | Batch frequency/ClinVar joins      |
| pyranges >=0.1,<1.0     | no       | [accelerated] | Interval overlap queries           |
| reportlab >=4.0,<5.0    | no       | [pdf]         | PDF report rendering               |
| weasyprint >=60.0,<63.0 | no       | [clinical]    | Clinical PDF rendering             |
| python-docx >=1.0,<2.0  | no       | [clinical]    | Clinical DOCX rendering            |
| httpx >=0.27,<1.0       | no       | [api]         | Remote API annotation              |

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
    knowledge/            # Gene-disease linkage (OMIM, ClinGen, HPO, constraint)
    prioritization/       # AF gating + CADD/REVEL/SpliceAI scoring (ScoreLoader)
    classification/       # ACMG evidence tagging
    structural/           # SV triage: parser, annotator, scorer, classifier, report
    mito/                 # Mitochondrial: genetic code, heteroplasmy, MITOMAP, classifier
    reporting/            # JSON, CSV, PDF, VCF (streaming writers)
        clinical/         # Clinical report generation (HTML/PDF/DOCX + audit trail)
    cohort/               # Multi-sample cohort analysis
    models/               # Dataclasses, enums, configs, warnings
    bundle/               # Score bundle downloader and transformer
    api/                  # Remote API annotation (VEP, ClinVar, CADD, SpliceAI)
    data/knowledge/       # Bundled gene-disease TSV files
    data/mito/            # Bundled mtDNA reference data (gene map, MITOMAP, HelixMTdb)
    _internal/            # Batch utils, interval tree, caching, vectorized ops
    py.typed              # PEP 561 marker
```

## Contributing

See [CONTRIBUTING.md](https://github.com/Behordeun/vartriage/blob/main/CONTRIBUTING.md).

## License

MIT
