# Validation

This document describes how to validate vartriage against benchmark datasets and details its clinical performance profile. Two complementary validation datasets are evaluated:

1. **ClinVar Expert-Panel (eRepo) Benchmark** — measuring multi-tier ACMG classification concordance, pathogenic sensitivity, missense sensitivity, positive predictive value (PPV), and Cohen's $\kappa$ against expert consensus.
2. **GIAB chr22 Benchmark** — genomic-scale validation measuring specificity, execution footprint, and runtime characteristics on whole-genome sequencing data.

## Headline Benchmark: ClinVar eRepo (v1.2.0 Shipped Release)

The v1.2.0 production evaluation reflects the frozen ClinGen SVI Bayesian point framework (Tavtigian et al. 2018, 2020) with PM2 calibrated to Supporting strength (+1 point) per ClinGen SVI 2020 recommendations and PVS1 Very Strong gated on established ClinGen loss-of-function disease mechanisms. The benchmark cohort comprises 6,476 expert-panel curated Pathogenic/Likely Pathogenic variants and 15,030 Benign/Likely Benign variants (total n = 21,506, ClinVar 3+ star reviews).

All metrics include two-sided 95% Wilson score confidence intervals:

| Metric | Value | 95% CI | n |
| --- | --- | --- | --- |
| Positive Predictive Value (PPV) | 0.916 | [0.907, 0.925] | 6,476 |
| Overall Pathogenic Sensitivity | 0.473 | [0.461, 0.485] | 6,476 |
| Truncating / Loss-of-Function Sensitivity | 0.911 | [0.898, 0.923] | 1,805 |
| Missense Sensitivity | 0.185 | [0.171, 0.200] | 3,922 |
| Splice-Site Sensitivity | 0.559 | [0.518, 0.600] | 567 |
| Cohen's $\kappa$ (unweighted, 5-tier) | 0.135 | [0.126, 0.144] | 21,506 |

### Clinical Positioning: Computational Triage Pre-filter and VUS Stratifier

VarTriage is engineered from first principles as an **audited computational triage pre-filter and VUS stratifier**, separating high-confidence candidates from variants requiring experimental or manual curation, rather than a standalone diagnostic classifier.

#### 1. High Positive Predictive Value (91.6% PPV)

When VarTriage renders a Pathogenic or Likely Pathogenic call, that call is backed by multi-layered, guideline-calibrated computational evidence. The tool prioritizes clinical precision over aggressive classification: clinical genetics workflows cannot afford high false-positive rates that trigger unwarranted follow-up or diagnostic misdirection.

#### 2. The Missense Bottleneck in Pure Computational Triage (18.5% Sensitivity)

A superficial reading of an 18.5% missense sensitivity might suggest underperformance. In reality, it represents rigorous adherence to ACMG/AMP and ClinGen SVI specifications.

Under the ClinGen SVI Bayesian point system:
- **Pathogenic:** $\ge +10$ points
- **Likely Pathogenic:** $+6$ to $+9$ points
- **Variant of Uncertain Significance (VUS):** $0$ to $+5$ points
- **Likely Benign:** $-1$ to $-6$ points
- **Benign:** $\le -7$ points

In a pure computational triage setting (processing raw VCFs without clinical phenotypes, pedigree records, or laboratory assays):
- An in-silico predictor (REVEL $\ge 0.75$ or $\ge 0.932$) earns PP3 at Moderate (+2) or Strong (+4) points.
- Population rarity (confirmed absence or gnomAD AF $< 10^{-5}$) earns PM2 at Supporting (+1) point per ClinGen SVI 2020.
- Combined, PP3 + PM2 yield $+3$ to $+5$ points.

Because $+5$ points remains below the $+6$ point threshold required for Likely Pathogenic, an isolated missense variant **cannot and should not cross into Likely Pathogenic** without orthogonal evidence:
- Experimental functional assays (**PS3**, +4 or +8 points)
- Case-level segregation, de novo occurrence, or parental phase (**PS2/PS4**, +2 to +8 points)
- Specific allelic recurrence of established pathogenic residues (**PS1/PM5**, +2 or +4 points via ClinVar protein index)

Tools that report missense sensitivities of 70–85% on automated benchmarks typically achieve this by:
1. Treating PM2 as Moderate (+2 points) or Strong (+4 points), directly violating ClinGen SVI calibration.
2. Inappropriately stacking multiple collinear in-silico tools (e.g. CADD + REVEL + PolyPhen) as independent lines of evidence.
3. Over-calling variants based on ClinVar assertions without verifying independent underlying criteria.

VarTriage strictly refuses to artificially escalate missense variants. Holding unverified missense variants at VUS is the mathematically and clinically sound behavior for an audited pre-filter.

#### 3. High Truncating Sensitivity (91.1%)

For null variants (nonsense and frameshift), loss of function is often the primary disease mechanism. When a gene possesses an established LoF mechanism (curated via ClinGen Haploinsufficiency Score 3 or 30), PVS1 fires at Very Strong (+8 points). Combined with PM2 at Supporting (+1 point) or protein length alteration (PM4, +2 points), the variant easily clears the Likely Pathogenic ($\ge +6$) or Pathogenic ($\ge +10$) threshold, achieving 91.1% sensitivity.

#### 4. Understanding the Concordance Metric: Cohen's $\kappa$ = 0.135

Unweighted Cohen's $\kappa$ measures chance-corrected categorical agreement across all 5 ACMG tiers. In a benchmark where the truth set consists of expert-panel resolved variants (which overwhelmingly fall in the definitive Pathogenic and Benign extremes) and an automated pre-filter deliberately holds variants lacking functional/clinical evidence in the middle tier (VUS), standard unweighted $\kappa$ penalizes every P $\to$ VUS or B $\to$ VUS triage holding as a 100% classification error.

Crucially:
- Directional inversions (calling a true Pathogenic variant Benign, or a true Benign variant Pathogenic) are virtually non-existent.
- In a clinical triage setting, routing an uncertain variant to VUS is safe: it signals that manual expert review or functional testing is required. Routing an uncertain variant to Likely Pathogenic on flimsy in-silico grounds is hazardous.
- When evaluated as a triage pre-filter (flagging variants with sufficient computational evidence to warrant immediate clinical review), the 91.6% PPV provides the necessary diagnostic certainty.

#### 5. Benign Sensitivity and Evidence Progress

Historically, automated benign classification was constrained by missing benign criteria. With the implementation of BP1 (missense in established LoF disease gene without pathogenic missense) and BP3 (in-frame indel in repetitive sequence), benign resolution continues to expand under the SVI framework alongside BS1/BS2/BA1. Complete benign resolution requires manual clinical correlation (BP6) and negative functional assays (BS3).

---

## Historical Benchmarks (Superseded Runs)

### Historical Run 1: v1.0.0 Full-Genome Release Run (PM2 Moderate Era)

The v1.0.0 release figures come from a single full-genome run over an expert-curated truth set of 18,166 pinned-ClinVar variants (GRCh38), using genome-wide REVEL, gnomAD exomes v4.1.1 frequencies (remote tabix, pinned cache), and GENCODE v46. This run predated the PM2 Supporting reconciliation and the PVS1 LoF mechanism gate.

| Metric | Value | 95% CI | n |
| --- | --- | --- | --- |
| PPV (P/LP call) | 0.999 | [0.997, 0.999] | 8,062 |
| Pathogenic sensitivity | 0.642 | [0.633, 0.650] | 12,545 |
| Specificity | 0.990 | [0.982, 0.994] | 1,187 |
| NPV | 0.207 | [0.197, 0.218] | 5,670 |

### Historical Run 2: v0.17.5 eRepo Run (Pre-Point-Engine)

The earlier eRepo figures below predate both the ClinGen SVI Bayesian point engine (the combining path since v0.18.3) and the v1.0.0 full-genome benchmark. They are retained for history and do not describe the current release.

| Benchmark | Variants | Pathogenic Sensitivity | Specificity | PPV | Runtime |
| --- | --- | --- | --- | --- | --- |
| GIAB chr22 | 50,284 | n/a | 71.8% | n/a | 30 s |
| ClinVar eRepo | 21,506 | 70.5% | n/a | 99.2% | n/a |

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

> **Historical (v0.17.5), superseded by the v1.2.0 eRepo benchmark at the top of this document.** The figures in this section predate the ClinGen SVI Bayesian point engine (the combining path since v0.18.3) and the v1.2.0 benchmark. They are retained to document the eRepo methodology and the pre-point-engine behavior, and do not describe the current release.

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

> The two tables above represent the historical v0.17.5 measurements, retained for longitudinal provenance. The current shipped engine is v1.2.0, which enforces ClinGen SVI point-system combining, PM2 at Supporting strength, and the ClinGen curated LoF mechanism gate for PVS1 Very Strong.

### Stratified Performance by Consequence (v1.2.0 Shipped Release)

On the ClinGen eRepo expert benchmark (6,476 Pathogenic / Likely Pathogenic variants):

| Consequence Category | Sensitivity | 95% CI | n | Rationale / Driver Criteria |
| --- | --- | --- | --- | --- |
| Truncating (Frameshift / Nonsense) | 91.1% | [0.898, 0.923] | 1,805 | PVS1 Very Strong (+8) gated by ClinGen LoF mechanism + PM2 Supporting (+1) |
| Splice-Site | 55.9% | [0.518, 0.600] | 567 | SpliceAI delta score $\ge 0.8$ triggers PVS1 Very Strong or Strong |
| Missense | 18.5% | [0.171, 0.200] | 3,922 | PP3 Moderate/Strong (+2/+4) + PM2 Supporting (+1) correctly bound at VUS ($\le +5$) |

### Interpretation (v1.2.0)

**High PPV (91.6%)** demonstrates that when VarTriage issues a Pathogenic or Likely Pathogenic call, it is highly dependable. The tool is deliberately conservative and precise: it avoids overcalling in-silico predictions.

**Moderate overall sensitivity (47.3%)** reflects the reality of automated ACMG/AMP classification under the ClinGen SVI Bayesian point system. The majority of clinically curated pathogenic variants in expert repositories are missense changes whose classifications were established via clinical case segregation (PS2/PS4) or in-vitro functional assays (PS3). Because an automated VCF triage tool cannot access medical records or execute wet-lab assays, it correctly and safely retains these variants at VUS rather than inflating computational scores.

**Benign classification dynamics:**
- BP1 (missense in a gene where truncating variants are the established mechanism and missense variants have not been proven pathogenic) and BP3 (in-frame indels in repetitive sequence) provide active benign evidence tags (-1 point each) under the SVI framework.
- Along with BA1 (stand-alone benign frequency), BS1 (frequency greater than expected), and BS2 (observed in healthy controls), benign sensitivity is calibrated to avoid false negatives.
- BP6 (reputable source reports benign without criteria) and BS3 (functional assay demonstrating no deleterious effect) require orthogonal curation and are not inferred computationally.

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
| VarTriage v1.2.0 (shipped) | 47.3% | 91.6% | Streaming, SVI point-engine, calibrated PM2 Supporting, LoF mechanism gate |
| VarTriage v0.17.5 (historical) | 70.5% | 99.2% | Streaming, pre-point-engine, uncalibrated PM2 Moderate |
| BIAS-2015 (reference) | ~85% | ~75% | Batch, requires InterVar |
| AutoACMG | ~78% | ~82% | Web service |

VarTriage is positioned as a high-precision, audited computational triage pre-filter and VUS stratifier with operational simplicity (single pip install, streaming architecture, no Java/Perl), not an uncalibrated diagnostic classifier. Higher apparent sensitivities in other tools frequently stem from uncalibrated point inflation or treating computational predictors as definitive clinical evidence.

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
