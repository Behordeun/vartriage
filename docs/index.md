# vartriage

A streaming pipeline for identifying and classifying pathogenic genetic variants from VCF data. Processes whole-genome scale files (4M+ variants) under 2GB memory via batched iterators.

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
