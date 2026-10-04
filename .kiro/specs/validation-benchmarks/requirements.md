# Validation and Benchmarks (v1.0.0) - Requirements

## Problem Statement

A variant classification tool used in clinical settings must demonstrate quantified accuracy against expert-curated truth sets before any laboratory can adopt it for patient reporting. Without validated sensitivity and specificity figures, clinical geneticists cannot assess whether the tool meets their minimum acceptable error rates, and laboratory accreditation bodies (CAP, CLIA) require documented analytical validation for any software used in the diagnostic workflow.

VarTriage currently passes 2,700+ internal tests that verify code consistency, but these do not measure whether the tool produces classifications that match expert human curators. The gap between "code works as designed" and "code produces clinically correct results" is what this spec closes.

## User Stories

### US-1: Clinical lab adoption decision

As a clinical lab director evaluating variant classification tools for CLIA-regulated reporting, I need published sensitivity, specificity, and PPV against expert-curated variants (ClinGen Evidence Repository), so I can determine whether vartriage meets our lab's minimum 95% PPV threshold before putting it into the diagnostic workflow.

### US-2: Per-criterion accuracy for troubleshooting

As a clinical bioinformatician investigating a discordant classification, I need to know which ACMG criteria have validated accuracy and which are known weak points, so I can assess whether the tool's output for a specific variant category is trustworthy or requires manual review.

### US-3: Consequence annotation accuracy

As a laboratory developing a validation plan for vartriage, I need the consequence caller (frameshift, missense, splice-site) verified against a truth set, so I can document that the upstream annotation feeding ACMG criteria is reliable.

### US-4: Performance relative to alternatives

As a genetic epidemiologist selecting tools for a 10,000-patient cohort, I need vartriage benchmarked against InterVar and BIAS-2015 on the same dataset, so I can make an evidence-based decision about which tool best fits my accuracy, throughput, and infrastructure requirements.

### US-5: Reproducible validation for regulatory compliance

As a quality officer maintaining ISO 15189 accreditation, I need validation scripts committed to the repository that can be re-executed on demand, so our annual revalidation audit can confirm that tool accuracy has not regressed since initial adoption.

---

## Functional Requirements

### FR-1: Truth set construction script

Script that:
1. Downloads ClinVar VCF (genome-wide, latest release).
2. Filters to review_status in ("reviewed by expert panel", "practice guideline").
3. Excludes conflicting interpretations.
4. Splits into Pathogenic/Likely_Pathogenic (positive) and Benign/Likely_Benign (negative).
5. Normalizes variants using the same normalizer as the pipeline.
6. Outputs: truth_set.vcf + truth_labels.tsv (chrom, pos, ref, alt, gene, clinvar_class).

### FR-2: BenchmarkRunner class

Automated benchmark execution:
1. Run vartriage on the truth set VCF with standard configuration.
2. Match pipeline output classifications to truth labels by coordinate.
3. Compute: sensitivity, specificity, PPV, NPV (with Wilson confidence intervals).
4. Compute Cohen's kappa for multi-class agreement (P, LP, VUS, LB, B).
5. Generate a 5x5 confusion matrix (predicted vs truth).

### FR-3: Per-criterion concordance

For each evidence tag:
1. Count how often it fires on truth-positive (P/LP) variants.
2. Count how often it fires on truth-negative (B/LB) variants.
3. Report: precision (tag fires AND variant is truly pathogenic / all times tag fires).
4. Identify criteria that fire falsely on benign variants.

### FR-4: Consequence accuracy comparison

1. Run the same variant set through local mode (v0.8.0 codon resolver) and API mode (VEP).
2. Tabulate concordance by consequence type.
3. Identify systematic misclassification patterns.
4. Report: overall concordance %, per-type concordance %.

### FR-5: Tool comparison (InterVar)

1. Install InterVar locally (or use its web output for a subset).
2. Run the same truth set through InterVar.
3. Compare classification output: vartriage vs InterVar vs ClinVar expert.
4. Report: classification agreement between tools, unique correct/incorrect calls per tool.

### FR-6: Validation report document

Fold the release validation report into `docs/validation.md` (single source of truth; a standalone `docs/validation_report.md` was intentionally not created to avoid figure drift):
1. Truth set description (source, size, composition, date).
2. Pipeline configuration used for benchmarking.
3. All metrics with 95% confidence intervals.
4. Per-criterion breakdown table.
5. Known limitations and their measured impact.
6. Comparison with other tools (if FR-5 completed).

### FR-7: REVEL coverage documentation

Quantify for the truth set:
1. What fraction of missense variants have REVEL scores.
2. PP3 fire rate with vs without REVEL.
3. Impact on sensitivity when REVEL is unavailable.
4. Recommendation for users without REVEL data.

---

## Non-Functional Requirements

### NFR-1: Reproducibility

The entire validation pipeline (download, filter, run, compute metrics) is scripted in `scripts/` and runnable with a single command. Pin the ClinVar release date for reproducibility.

### NFR-2: CI integration

Add a `make validate` target (or equivalent) that runs the benchmark on a small subset (~1000 variants) as a regression check. Full validation (~50K variants) runs on-demand or nightly.

### NFR-3: Publication quality

Metrics must be computed with proper statistical methods (Wilson CI, not normal approximation). Figures are generated as publication-ready plots (matplotlib or equivalent).
