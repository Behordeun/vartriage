# Sample Output Files

This directory contains representative output files from vartriage runs.
Use these to understand the format and structure of each output type
without needing to run the full pipeline.

## Files

| File                                     | Description                                                       |
| ---------------------------------------- | ----------------------------------------------------------------- |
| `sample_clinical_report.html`            | Self-contained clinical HTML report for a hereditary cancer panel |
| `sample_clinical_report.html.audit.json` | Audit trail sidecar for the clinical report                       |
| `sample_pipeline_output.json`            | Standard JSON output from a prioritized exome run                 |
| `sample_pipeline_output.csv`             | CSV equivalent of the JSON output                                 |

## How these were generated

All samples use synthetic variant data modeled on real-world exome
sequencing patterns. Patient identifiers are fictional. The pipeline
configuration and reference file checksums are representative of a
typical clinical-grade setup.

### Clinical report

```bash
vartriage \
  --vcf synthetic_exome.vcf.gz \
  --output sample_clinical_report.html \
  --output-format clinical-html \
  --patient-id "DEMO-2026-001" \
  --panel-name "Hereditary Cancer Panel v2" \
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
  --gene-annotation refs/gencode.v44.gtf \
  --gnomad refs/gnomad.v4.exomes.tsv \
  --clinvar refs/clinvar.tsv \
  --cadd-scores refs/cadd_v1.7.tsv \
  --revel-scores refs/revel_v1.3.tsv
```
