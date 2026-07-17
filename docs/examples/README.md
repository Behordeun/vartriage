# Sample Output Files

This directory contains representative output files from vartriage runs.
Use these to understand the format and structure of each output type
without needing to run the full pipeline.

## Files

| File                                     | Description                                                           |
| ---------------------------------------- | --------------------------------------------------------------------- |
| `sample_clinical_report.html`            | Self-contained clinical HTML report for a hereditary cancer panel     |
| `sample_clinical_report.html.audit.json` | Audit trail sidecar for the clinical report                           |
| `sample_pipeline_output.json`            | Standard JSON output with gene_name, revel_score, and all ACMG fields |
| `sample_pipeline_output.csv`             | CSV equivalent (12 columns: chromosome through evidence_tags)         |

## How these were generated

All samples use synthetic variant data modeled on real-world exome
sequencing patterns. Patient identifiers are fictional. The pipeline
configuration and reference file checksums are representative of a
typical clinical-grade setup.

### Clinical report

Using installed bundles (recommended):

```bash
vartriage bundle download --bundle clinvar
vartriage bundle download --bundle gnomad-exomes-chr22
vartriage bundle download --bundle gencode

vartriage \
  --vcf synthetic_exome.vcf.gz \
  --output sample_clinical_report.html \
  --output-format clinical-html \
  --patient-id "DEMO-2026-001" \
  --panel-name "Hereditary Cancer Panel v2" \
  --use-bundles \
  --gene-list refs/hereditary_cancer_panel.txt
```

Or with explicit paths:

```bash
vartriage \
  --vcf synthetic_exome.vcf.gz \
  --output sample_clinical_report.html \
  --output-format clinical-html \
  --patient-id "DEMO-2026-001" \
  --panel-name "Hereditary Cancer Panel v2" \
  --reference-fasta refs/GRCh38.fa \
  --gene-annotation refs/gencode.v44.gtf \
  --gnomad refs/gnomad.v4.exomes.tsv \
  --clinvar refs/clinvar.tsv \
  --cadd-scores refs/cadd_v1.7.tsv \
  --revel-scores refs/revel_v1.3.tsv \
  --spliceai-scores refs/spliceai_scores.tsv \
  --gene-list refs/hereditary_cancer_panel.txt
```

### Standard JSON/CSV

```bash
vartriage \
  --vcf synthetic_exome.vcf.gz \
  --output sample_pipeline_output.json \
  --output-format json \
  --reference-fasta refs/GRCh38.fa \
  --gene-annotation refs/gencode.v44.gtf \
  --gnomad refs/gnomad.v4.exomes.tsv \
  --clinvar refs/clinvar.tsv \
  --cadd-scores refs/cadd_v1.7.tsv \
  --revel-scores refs/revel_v1.3.tsv
```

### API mode (no local reference files)

```bash
pip install vartriage[api]

vartriage \
  --vcf synthetic_exome.vcf.gz \
  --output sample_api_output.json \
  --output-format json \
  --mode api
```

API mode queries Ensembl VEP, ClinVar, and CADD via HTTP. Responses
are cached locally in SQLite (`~/.vartriage/api_cache.db`) so repeated
runs complete instantly. See [docs/api-mode.md](../api-mode.md) for
full configuration and performance details.

### Hybrid mode (local gnomAD + API for the rest)

```bash
vartriage \
  --vcf synthetic_exome.vcf.gz \
  --output sample_hybrid_output.json \
  --output-format json \
  --mode hybrid \
  --gnomad refs/gnomad.v4.exomes.tsv \
  --revel-scores refs/revel_v1.3.tsv
```

Local files take priority where provided. The pipeline queries APIs
only for sources without a local path.
