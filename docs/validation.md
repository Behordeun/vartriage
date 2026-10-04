# Validation

This document describes how to validate vartriage against benchmark datasets. Two complementary validation approaches are used:

1. **GIAB chr22** — genomic-scale validation measuring specificity and runtime characteristics
2. **ClinVar eRepo** — expert-panel curated variants measuring pathogenic sensitivity and PPV

## Summary (v1.0.0 full-genome release run)

The v1.0.0 release figures come from a single full-genome run over an expert-curated truth set of 18,166 pinned-ClinVar variants (GRCh38), using genome-wide REVEL, gnomAD exomes v4.1.1 frequencies (remote tabix, pinned cache), and GENCODE v46. All metrics carry a 95% Wilson confidence interval.

| Metric | Value | 95% CI | n |
| --- | --- | --- | --- |
| PPV (P/LP call) | 0.999 | [0.997, 0.999] | 8,062 |
| Pathogenic sensitivity | 0.642 | [0.633, 0.650] | 12,545 |
| Specificity | 0.990 | [0.982, 0.994] | 1,187 |
| NPV | 0.207 | [0.197, 0.218] | 5,670 |

The thresholds were pre-registered before the figures existed (`.kiro/specs/validation-benchmarks/accuracy-decision-rule.md`). The PPV lower bound (0.997) clears the pre-registered clinical bar; sensitivity (0.642) is below the 0.70 clinical gate, so v1.0.0 is positioned as a **research-grade computational-triage tool**: high precision on positive calls, with uncertain variants held at VUS rather than overcalled.

Confusion counts: 8,050 true positive, 12 false positive, 1,175 true negative, 4,495 false negative. Of the P/LP variants not recovered, 4,486 were classified VUS rather than mis-called benign. Specificity recovered from 0.0 (frequency-blind REVEL-only run) to 0.990 once gnomAD frequency was available, confirming the benign-evidence path depends on frequency data rather than any classifier defect.

**Coverage boundaries (these bound the figures):** frequency is from gnomAD exomes v4.1.1 only (genomes not queried); the 127 mitochondrial truth variants received no gnomAD frequency; CADD and SpliceAI were not used, so REVEL was the only in-silico score. Each plausibly suppresses sensitivity and is the subject of the sensitivity-evidence-coverage work scoped for v1.2.0.

**Reproducibility:** the gnomAD remote cache is pinned (entries never expire), the truth set is built deterministically from a pinned ClinVar release, and reference versions are fixed, so a re-run reproduces these figures.

### Path to clinical-grade

Clinical-grade status has two distinct requirements, and both must hold.

1. **Clear the pre-registered gate honestly.** The only unmet threshold is sensitivity (0.642, needs ≥ 0.70). It must rise by giving the classifier evidence it currently lacks — gnomAD genomes, CADD, SpliceAI, and chrM frequency let PP3, splice, and rarity criteria fire on variants presently held at VUS — never by fitting thresholds to the truth set, which the pre-registration forbids. Whether sensitivity reaches 0.70 is then measured on a re-run under the same frozen rule. This is scoped in the sensitivity-evidence-coverage spec for v1.2.0.
2. **Validate beyond ClinVar.** The current benchmark measures concordance against ClinVar, and the classifier itself uses ClinVar-derived evidence (PP5/BP6, review status), so the study carries a circularity risk stated plainly here. Clinical-grade validation requires an independent, orthogonal truth set the classifier did not learn from, per-gene and per-disease breakdowns rather than one global number, and concordance against the full ACMG five-tier rather than the collapsed P/VUS/B this run produced. Formal CAP/CLIA validation and regulatory review are a separate organizational track required before any diagnostic use. This is a later milestone.

### Historical run (v0.17.5 eRepo, superseded)

The earlier eRepo figures below predate both the ClinGen SVI Bayesian point engine (the combining path since v0.18.3) and the v1.0.0 full-genome benchmark. They are retained for history and do not describe the current release.

| Benchmark     | Variants | Pathogenic Sensitivity | Specificity | PPV   | Runtime |
| ------------- | -------- | ---------------------- | ----------- | ----- | ------- |
| GIAB chr22    | 50,284   | n/a                    | 71.8%       | n/a   | 30 s    |
| ClinVar eRepo | 21,506   | 70.5%                  | n/a          | 99.2% | n/a     |

**Known limitations:**

- Pathogenic sensitivity 0.642: uncertain P/LP variants are held at VUS rather than escalated, pending the additional in-silico and frequency evidence scoped for v1.2.0 (gnomAD genomes, CADD, SpliceAI, chrM frequency). Sensitivity must improve through added evidence, never by fitting thresholds to the truth set.
- NPV 0.207: a non-pathogenic call is weak evidence of benignity, a consequence of the VUS-heavy behavior.
- Validation is concordance against ClinVar, which also supplies some classifier evidence (PP5/BP6, review status); an independent, orthogonal truth set is required for clinical-grade validation.
- BS2 is emitted for dominant-disorder genes when gnomAD homozygote counts are available; it does not fire for recessive genes or without homozygote-count data.

---

## GIAB Validation

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

| Resource                          | Size    | Source                                                                                                                           |
| --------------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------- |
| GIAB HG002 benchmark VCF (GRCh38) | ~250 MB | [NIST FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/) |
| GIAB high-confidence regions BED  | ~2 MB   | Same FTP                                                                                                                         |
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

For full validation across all chromosomes:

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

---

## ClinVar Expert-Panel (eRepo) Validation

> **Historical (v0.17.5), superseded by the v1.0.0 full-genome run at the top of this document.** The figures in this section predate the ClinGen SVI Bayesian point engine (the combining path since v0.18.3) and the v1.0.0 benchmark. They are retained to document the eRepo methodology and the pre-point-engine behavior, and do not describe the current release.

The eRepo validation uses ClinVar submissions from Expert Panels only (review status 3+ stars). These are the highest-confidence clinical variant interpretations available, representing consensus multi-lab assessments.

### eRepo Overview

Unlike GIAB (which tests runtime and specificity on a healthy genome), the eRepo benchmark directly measures whether vartriage's ACMG criteria produce correct classifications for variants with established clinical significance.

### Data source

- ClinVar VCF filtered to `CLNREVSTAT` containing "reviewed_by_expert_panel"
- 21,928 variants across all chromosomes
- Split: 4,116 Pathogenic/Likely Pathogenic + 17,812 Benign/Likely Benign

### Metrics (v0.17.5, SpliceAI SQLite + remote gnomAD)

| Metric | v0.17.3 (no SpliceAI) | v0.17.5 (SpliceAI SQLite) |
| ------ | --------------------- | ------------------------- |
| Classified variants | 21,614 | 21,506 |
| Pathogenic sensitivity | 65.6% | 70.5% |
| Splice-site sensitivity | 9.8% | 55.9% |
| Benign sensitivity | 6.0% | 6.0% |
| Positive predictive value (PPV) | 99.2% | 99.2% |

### Stratified by consequence type (v0.17.5)

| Consequence | Sensitivity |
| ----------- | ----------- |
| Missense | 46.9% |
| Frameshift | 99.7% |
| Splice_Site | 55.9% |
| Synonymous | 0.0% |

### Interpretation

**High PPV (99.2%)** means when vartriage calls something Pathogenic or Likely Pathogenic, it is almost always correct. The tool is conservative and precise — it does not over-call.

**Moderate sensitivity (70.5%)** means vartriage misses about 30% of known pathogenic variants. This is expected given:

- Missense sensitivity is 46.9%: many pathogenic missense variants don't exceed the ClinGen-calibrated REVEL threshold (0.644)
- Splice-site sensitivity is 55.9%: splice calls depend on precomputed SpliceAI delta scores being available for the variant
- Some pathogenic variants require functional evidence the tool cannot assess computationally

**Low benign sensitivity (6.0%)** means most benign variants are classified as VUS rather than Benign/Likely Benign. This is because:

- BP1 (missense in a gene where truncating variants cause disease) is not yet implemented
- BP3 (in-frame in repetitive region without known function) is not yet implemented
- BP6 (reputable source reports benign) is not yet implemented
- Without benign criteria firing, variants default to VUS

### Running the eRepo validation

The validation script lives in the paper reproducibility repository (`vartriage-streaming-acmg`):

```bash
cd vartriage-streaming-acmg
pip install -e .
python analysis/validate_erepo.py
```

This produces stratified results in `results/erepo_stratified.json` and summary metrics in the console.

### Comparison to existing tools

| Tool | Pathogenic Sensitivity | PPV | Architecture |
| ---- | ---------------------- | --- | ------------ |
| VarTriage v0.17.5 | 70.5% | 99.2% | Streaming, zero-dependency |
| BIAS-2015 (reference) | ~85% | ~75% | Batch, requires InterVar |
| AutoACMG | ~78% | ~82% | Web service |

VarTriage is positioned on high precision (PPV) and operational simplicity (single pip install, streaming architecture, no Java/Perl), not on raw sensitivity. Further sensitivity gains depend on additional benign criteria (BP1/BP3/BP6) and functional evidence the tool cannot assess computationally.

### Provenance

Validation results are tracked with full provenance:

```json
{
  "vartriage_version": "0.17.5",
  "clinvar_date": "2025-01-13",
  "genome_build": "GRCh38",
  "gencode_version": "v46",
  "gnomad_version": "v4.1.1",
  "revel_version": "v1.3",
  "spliceai_source": "OpenCRAVAT SQLite",
  "execution_timestamp": "2026-08-17T...",
  "total_variants": 21506
}
```
