# vartriage

A streaming computational triage pre-filter and VUS stratification pipeline for identifying, prioritizing, and classifying genetic variants from VCF data. Processes whole-genome scale files (4M+ variants) under 2GB memory via batched iterators.

VarTriage operates as an audited computational pre-filter for clinical and research genomics pipelines. By implementing the ClinGen SVI Bayesian point system, it maintains high positive predictive value (91.6% PPV) on positive calls while conservatively holding variants lacking experimental functional assays or patient-level segregation data at VUS rather than risking clinical overcalls.

Reads a VCF, applies quality filters, annotates functional consequence and population frequency, computes pathogenicity scores, runs ACMG/AMP evidence classification, and writes a ranked candidate list in JSON, CSV, or PDF.

## Install

```bash
pip install vartriage
```

## Quick start

```python
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, AnnotationConfig

config = PipelineConfig(
    vcf_path=Path("sample.vcf.gz"),
    output_path=Path("candidates.json"),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
    ),
)

pipeline = Pipeline(config)
pipeline.run()
```

See [Getting Started](getting-started.md) for installation options and a full walkthrough.

For zero-config annotation without local reference files, see [API Mode](api-mode.md).

For multi-sample cohort analysis (shared variants, gene burden), see [Cohort Analysis](cohort-analysis.md).

For mitochondrial DNA variant analysis (heteroplasmy, MITOMAP, mtDNA-specific classification), see [Mitochondrial Variants](mitochondrial.md).

For pre-flight sample quality control (Ti/Tv, het/hom, variant count sanity checks with a strict gate), see [Quality Control](quality-control.md).

For the full ACMG/AMP classification reference (criteria, thresholds, PS1/PM5 setup, combining rules), see [ACMG Criteria](acmg-criteria.md).
