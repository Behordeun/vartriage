# Getting Started

## Requirements

- Python 3.10+
- A VCF file (`.vcf` or `.vcf.gz` with `.tbi` index)
- Reference files for annotation (gnomAD, gene models, optionally ClinVar)

## Installation

Base install (pure-Python backends, no optional dependencies):

```bash
pip install vartriage
```

With faster annotation backends (polars for batch joins, pyranges for interval queries):

```bash
pip install vartriage[accelerated]
```

With PDF report support:

```bash
pip install vartriage[pdf]
```

Everything:

```bash
pip install vartriage[all]
```

## Quick example

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

The output is a JSON file at `candidates.json` with variants ranked by composite pathogenicity score, classified per ACMG/AMP 2015 guidelines.

## What you need

| Item | Purpose | Where to get it |
|------|---------|-----------------|
| VCF file | Input variants | Your sequencing pipeline (GATK, DeepVariant, etc.) |
| Gene annotation (GTF/GFF) | Functional consequence assignment | [GENCODE](https://www.gencodegenes.org/human/) |
| gnomAD frequency file | Population frequency filtering | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) |
| ClinVar file (optional) | Clinical significance lookup | [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/) |
| CADD scores (optional) | Pathogenicity scoring | [CADD download](https://cadd.gs.washington.edu/download) |
| REVEL scores (optional) | Pathogenicity scoring | [REVEL download](https://sites.google.com/site/revelgenomics/) |

See [Reference Files](reference-files.md) for column format requirements and preparation steps.
