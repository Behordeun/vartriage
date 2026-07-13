# GIAB Validation

This document describes how to validate vartriage against Genome in a Bottle (GIAB) benchmark data, the gold standard for germline variant calling evaluation.

## Overview

The validation uses GIAB sample HG002 (Ashkenazi Jewish trio son) with the v4.2.1 benchmark set on GRCh38. We run vartriage on the benchmark VCF and compare classification results against ClinVar assertions to measure concordance.

This is not a variant *calling* benchmark (vartriage does not call variants). It validates the annotation, prioritization, and ACMG classification stages against known-truth clinical significance.

## What we measure

| Metric                      | Definition                                                                                                                  |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Sensitivity vs ClinVar      | Fraction of ClinVar Pathogenic/Likely Pathogenic variants that vartriage also classifies as Pathogenic or Likely Pathogenic |
| Classification distribution | Breakdown of P / LP / VUS across the variant set                                                                            |
| Evidence tag coverage       | How often each ACMG criterion (PVS1, PM2, PP3, PP5) fires                                                                   |
| Missing data rate           | Fraction of variants lacking gnomAD, ClinVar, or predictor scores                                                           |

## Quick start

```bash
# Install dependencies
pip install vartriage[all]
brew install bcftools htslib wget  # macOS
# or: apt-get install bcftools tabix wget  # Ubuntu/Debian

# Run the validation script
./scripts/validate_giab.sh --output-dir validation_results
```

The script handles downloading, extraction, running the pipeline, and computing metrics. Expect ~20 GB of downloads on first run (gnomAD chr22 is the largest file). Subsequent runs reuse cached downloads.

## Data sources

| Resource                          | Size    | Source                                                                                                                          |
| --------------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------- |
| GIAB HG002 benchmark VCF (GRCh38) | ~250 MB | [NIST FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/) |
| GIAB high-confidence regions BED  | ~2 MB   | Same FTP                                                                                                                        |
| ClinVar VCF (GRCh38)              | ~80 MB  | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)                                                                 |
| gnomAD v4.1.1 exomes chr22        | ~4.7 GB | [gnomAD Downloads](https://gnomad.broadinstitute.org/downloads)                                                                  |

The script downloads chr22 only for a faster validation cycle. For full-genome validation, modify the script to download the complete gnomAD exomes VCF and remove the `-r chr22` filter.

## Expected results (chr22)

Based on GIAB HG002 chr22 with ClinVar and gnomAD annotation:

- ~85,000 variants pass quality filters in the high-confidence regions
- 5-15 variants have ClinVar Pathogenic or Likely Pathogenic assertions
- Pipeline sensitivity vs ClinVar should be >90% for P/LP concordance
- Most variants classify as VUS (expected for a healthy individual's germline)

## Interpreting validation output

The script produces four output files:

### validation_metrics.json

```json
{
  "validation_summary": {
    "total_variants_processed": 84231,
    "classification_distribution": {
      "VUS": 84210,
      "Likely_Pathogenic": 14,
      "Pathogenic": 7
    }
  },
  "concordance": {
    "clinvar_actionable_count": 12,
    "pipeline_actionable_count": 21,
    "concordant_count": 11,
    "sensitivity_vs_clinvar": 0.917
  }
}
```

### Key considerations

**Why sensitivity might be below 100%:**

- ClinVar assertions based on functional studies that the pipeline cannot replicate computationally
- Variants where CADD/REVEL scores are below threshold despite clinical evidence
- Variants in regions not covered by the predictor score files

**Why pipeline might flag variants ClinVar does not:**

- PM2 fires for rare variants absent from gnomAD, even without ClinVar annotation
- PP3 fires on high predictor scores regardless of ClinVar status
- These are *candidates* for clinical review, not false positives per se

## Full-genome validation

For comprehensive validation across all chromosomes:

```bash
# Modify the script to use full gnomAD and CADD/REVEL scores
# Warning: requires ~500 GB disk space and several hours
export FULL_GENOME=true
./scripts/validate_giab.sh --output-dir full_validation
```

You will also need:

- CADD whole-genome pre-scored file (~350 GB uncompressed)
- REVEL scores file (~2 GB)
- SpliceAI pre-computed scores (~30 GB)
- GENCODE v44 GTF (~1.5 GB)

## Trio validation

To validate inheritance pattern classification, use the full Ashkenazi trio:

```bash
vartriage \
  --vcf giab_trio_merged.vcf.gz \
  --proband HG002 \
  --mother HG004 \
  --father HG003 \
  --output trio_validation.json \
  --output-format json \
  --gnomad refs/gnomad.v4.exomes.vcf.bgz \
  --clinvar refs/clinvar.tsv \
  --cadd-scores refs/cadd_v1.7.tsv \
  --revel-scores refs/revel_v1.3.tsv
```

Expected: de novo variants in HG002 should be correctly identified when comparing against parental genotypes.

## CI integration

The validation script can run in CI with a cached reference data volume:

```yaml
# .github/workflows/validation.yml (manual trigger)
name: GIAB Validation
on: workflow_dispatch
jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - run: pip install -e .[all]
      - run: sudo apt-get install -y bcftools tabix
      - run: ./scripts/validate_giab.sh
      - uses: actions/upload-artifact@v4
        with:
          name: validation-results
          path: validation_results/
```

This runs on-demand rather than on every PR (too slow and too much data).
