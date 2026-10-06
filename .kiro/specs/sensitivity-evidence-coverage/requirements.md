# Sensitivity Evidence Coverage (v1.2.0) - Requirements

## Problem Statement

The v1.0.0 full-genome validation established a strong PPV (0.999, 95% CI [0.997, 0.999]) but a pathogenic sensitivity of 0.642 (95% CI [0.633, 0.650]), below the pre-registered clinical-grade gate of 0.70. The validation report shows the shortfall is overwhelmingly P/LP variants held at VUS for lack of evidence, not variants driven to the wrong side: of the P/LP variants not recovered, 4,486 were VUS and only a handful were called benign.

The v1.0.0 run used a deliberately narrow evidence set, and each omission plausibly suppresses sensitivity:

- Allele frequency from gnomAD exomes v4.1.1 only; genomes were not queried, so deep-intronic and regulatory variants received no frequency evidence.
- No CADD and no SpliceAI; REVEL was the only in-silico score, so splice and non-missense in-silico evidence did not contribute to PP3.
- The 127 mitochondrial truth variants received no gnomAD frequency.

This spec adds that missing evidence so the frequency-driven and in-silico criteria can fire on variants presently held at VUS. The goal is a higher, honestly-earned sensitivity measured on a re-run of the frozen validation-benchmarks truth set.

## Hard constraint (carried from the validation pre-registration)

Sensitivity must improve only by supplying the classifier with evidence it currently lacks. Combining-path thresholds must never be fitted to the truth set. The acceptance measurement re-runs the existing benchmark under the pre-registered `accuracy-decision-rule.md` unchanged; this spec does not edit that rule. Every reported metric carries a confidence interval, and any evidence source that remains unavailable is recorded as "not performed," never estimated.

## User Stories

### US-1: Splice and non-missense in-silico evidence

As a clinical scientist, I need splice-region and non-missense variants to receive in-silico pathogenicity evidence (SpliceAI, CADD) in addition to REVEL, so that truncating-adjacent and splice variants currently held at VUS can escalate when the evidence supports it.

### US-2: Whole-genome frequency coverage

As an analyst running whole-genome input, I need allele frequencies for deep-intronic and regulatory variants that gnomAD exomes does not cover, by adding gnomAD genomes, so the rarity (PM2) and benign-frequency (BA1/BS1/BS2) criteria apply genome-wide rather than exome-only.

### US-3: Mitochondrial frequency

As a mitochondrial-disease specialist, I need chrM variants to receive population frequency (from a mitochondrial frequency source, since gnomAD exomes publishes no chrM file), so the 127 chrM truth variants are no longer frequency-blind.

### US-4: Re-measured sensitivity under the frozen rule

As the release owner, I need the full-genome benchmark re-run with the expanded evidence and the result read against the unchanged pre-registered rule, so I can see whether sensitivity legitimately reaches the 0.70 clinical gate without any threshold change.

## Acceptance

- gnomAD genomes v4.1.1 queryable alongside exomes through the existing remote-tabix backend and population cache, with the exome/genome source combination recorded in the run provenance.
- CADD and SpliceAI contribute to PP3/BP4 and splice evidence through the existing score backends.
- A mitochondrial frequency source supplies chrM frequencies.
- The frozen validation-benchmarks truth set is re-run with the expanded evidence; the report states the new sensitivity with its CI and the exome-vs-genome, CADD, SpliceAI, and chrM coverage deltas against the v1.0.0 run.
- No change to `accuracy-decision-rule.md` and no threshold fitting against the truth set.

## Out of scope

Independent (non-ClinVar) truth-set validation and per-gene/per-disease breakdowns are a separate, later milestone, as is any formal CAP/CLIA or regulatory work. This spec is evidence coverage and the honest re-measurement of sensitivity only.
