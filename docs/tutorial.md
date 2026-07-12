# Tutorial

A complete analysis walkthrough. Exome sequencing on a patient sample, identifying pathogenic variants that may explain a rare genetic disorder.

## Scenario

You have a proband with suspected hereditary cardiomyopathy. The clinical lab ran exome sequencing and delivered a VCF file. You want to identify candidate pathogenic variants, ranked by severity.

## 1. Load the VCF

```python
from pathlib import Path
from vartriage import VCFParser

vcf_path = Path("proband_exome.vcf.gz")

with VCFParser(vcf_path) as parser:
    count = sum(1 for _ in parser)
    print(f"Total variants: {count}")
```

The parser streams records one at a time. Memory stays flat regardless of file size.

## 2. Configure the pipeline

```python
from vartriage import (
    Pipeline, PipelineConfig, AnnotationConfig,
    PrioritizationConfig, QualityFilterConfig, ReportConfig,
)

config = PipelineConfig(
    vcf_path=Path("proband_exome.vcf.gz"),
    output_path=Path("results/candidates.json"),
    quality_filter=QualityFilterConfig(min_qual=30.0),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("references/gencode.v44.gtf"),
        gnomad_path=Path("references/gnomad.v4.exomes.tsv"),
        clinvar_path=Path("references/clinvar_20240101.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.001,  # stringent for rare disease
        cadd_scores_path=Path("references/cadd_v1.7.tsv"),
        revel_scores_path=Path("references/revel_v1.3.tsv"),
    ),
    report=ReportConfig(output_format="json"),
)
```

Configuration validates at construction. If a reference file path does not exist, you get a `FileNotFoundError` immediately.

## 3. Run the pipeline

```python
pipeline = Pipeline(config)
output_path = pipeline.run()
print(f"Report written to: {output_path}")
```

The pipeline wires stages sequentially:

```text
VCFParser → QualityFilter → AnnotationEngine → PrioritizationEngine → ACMGClassifier → ReportGenerator
```

## 4. Interpret the output

The JSON output is an array of classified variants, sorted by composite pathogenicity rank (highest first):

```json
[
  {
    "chromosome": "chr11",
    "position": 47332400,
    "ref_allele": "G",
    "alt_allele": "A",
    "functional_consequence": "Nonsense",
    "allele_frequency": 0.000012,
    "composite_rank": 0.92,
    "clinvar_assertion": "Pathogenic",
    "acmg_classification": "Pathogenic",
    "evidence_tags": ["PVS1", "PM2", "PP5"]
  },
  {
    "chromosome": "chr1",
    "position": 237778000,
    "ref_allele": "C",
    "alt_allele": "T",
    "functional_consequence": "Missense",
    "allele_frequency": null,
    "composite_rank": 0.78,
    "clinvar_assertion": null,
    "acmg_classification": "VUS",
    "evidence_tags": ["PP3"]
  }
]
```

Key fields:

- `composite_rank`: 0 to 1, higher means more likely pathogenic
- `evidence_tags`: which ACMG criteria were satisfied
- `acmg_classification`: final call (Pathogenic, Likely_Pathogenic, or VUS)
- `allele_frequency`: null means the variant wasn't found in gnomAD

## 5. Check the warning accumulator

After a run, inspect what data was missing:

```python
acc = pipeline.warning_accumulator
print(f"Total missing data events: {acc.total_count}")
print(f"Sources with missing data: {acc.sources}")
```

The counts show how many variants lacked gnomAD frequency or ClinVar assertions. A high gnomAD miss count may point to an incomplete reference file for the target regions.

## 6. CLI alternative

The same analysis from the command line:

```bash
vartriage \
  --vcf proband_exome.vcf.gz \
  --output results/candidates.json \
  --output-format json \
  --gene-annotation references/gencode.v44.gtf \
  --gnomad references/gnomad.v4.exomes.tsv \
  --clinvar references/clinvar_20240101.tsv \
  --cadd-scores references/cadd_v1.7.tsv \
  --revel-scores references/revel_v1.3.tsv \
  --spliceai-scores references/spliceai_scores.tsv
  --gene-list references/cardiac_panel.txt
  --spliceai-scores references/spliceai_scores.tsv
```

Output is identical to the Python API: same JSON structure, same ranking.

## 7. Working with the output

### Load and inspect results

```python
import json
from pathlib import Path

results = json.loads(Path("results/candidates.json").read_text())
print(f"{len(results)} classified variants")
```

### Filter to actionable variants

```python
actionable = [
    v for v in results
    if v["acmg_classification"] in ("Pathogenic", "Likely_Pathogenic")
]
print(f"{len(actionable)} pathogenic/likely pathogenic variants")
```

### Export top candidates to CSV

```python
import csv

fields = [
    "chromosome", "position", "ref_allele", "alt_allele",
    "functional_consequence", "allele_frequency",
    "composite_rank", "acmg_classification",
]

with open("results/top_candidates.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(actionable)
```

## 8. Without annotation (basic QC mode)

You don't always need the full reference stack. If you just want to parse a VCF and run quality filtering (maybe to validate a file before a bigger run), skip the annotation flags entirely:

```bash
vartriage --vcf sample.vcf.gz --output qc_output.json
```

Or in Python:

```python
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, ReportConfig

config = PipelineConfig(
    vcf_path=Path("sample.vcf.gz"),
    output_path=Path("qc_output.json"),
    report=ReportConfig(output_format="json"),
)

pipeline = Pipeline(config)
pipeline.run()
```

Without annotation references, variants pass through with `Intergenic` consequence and null frequencies. Downstream stages still run. You get ACMG classifications based on available evidence (which will be minimal). Useful for verifying your VCF parses cleanly and checking variant counts before investing in the full annotation step.

## Using individual stages

Each stage works independently:

```python
from vartriage import VCFParser, QualityFilter, QualityFilterConfig

with VCFParser(Path("input.vcf.gz")) as parser:
    qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
    passing = list(qf.apply(iter(parser)))
    print(f"{len(passing)} variants passed quality filter")
```

## 9. Generate a clinical report

When you need a structured, sign-off-ready document for clinical review, use the clinical report output format. This produces a report with per-variant evidence narratives, ACMG criteria explanations, and a JSON audit trail.

### CLI

```bash
vartriage \
  --vcf proband_exome.vcf.gz \
  --output results/clinical_report.html \
  --output-format clinical-html \
  --patient-id PAT-2025-001 \
  --panel-name "Cardiac Panel v3" \
  --gene-annotation references/gencode.v44.gtf \
  --gnomad references/gnomad.v4.exomes.tsv \
  --clinvar references/clinvar_20240101.tsv \
  --cadd-scores references/cadd_v1.7.tsv \
  --revel-scores references/revel_v1.3.tsv \
  --spliceai-scores references/spliceai_scores.tsv \
  --gene-list references/cardiac_panel.txt
```

This produces two files:

- `results/clinical_report.html`: self-contained HTML report viewable in any browser without network access
- `results/clinical_report.html.audit.json`: machine-parseable audit trail with run manifest and decision log

### Python API

```python
from pathlib import Path
from vartriage import (
    Pipeline, PipelineConfig, AnnotationConfig,
    PrioritizationConfig, QualityFilterConfig, ReportConfig,
)
from vartriage.models.config import ClinicalReportConfig

config = PipelineConfig(
    vcf_path=Path("proband_exome.vcf.gz"),
    output_path=Path("results/clinical_report.html"),
    quality_filter=QualityFilterConfig(min_qual=30.0),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("references/gencode.v44.gtf"),
        gnomad_path=Path("references/gnomad.v4.exomes.tsv"),
        clinvar_path=Path("references/clinvar_20240101.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.001,
        cadd_scores_path=Path("references/cadd_v1.7.tsv"),
        revel_scores_path=Path("references/revel_v1.3.tsv"),
        spliceai_scores_path=Path("references/spliceai_scores.tsv"),
    ),
    report=ReportConfig(output_format="clinical-html"),
    clinical_report=ClinicalReportConfig(
        patient_id="PAT-2025-001",
        panel_name="Cardiac Panel v3",
        output_format="clinical-html",
    ),
)

pipeline = Pipeline(config)
pipeline.run()
```

### PDF output

For a printable PDF, change the format and install weasyprint:

```bash
pip install weasyprint

vartriage \
  --vcf proband_exome.vcf.gz \
  --output results/clinical_report.pdf \
  --output-format clinical-pdf \
  --patient-id PAT-2025-001 \
  --panel-name "Cardiac Panel v3" \
  --gene-annotation references/gencode.v44.gtf \
  --gnomad references/gnomad.v4.exomes.tsv \
  --clinvar references/clinvar_20240101.tsv \
  --cadd-scores references/cadd_v1.7.tsv \
  --revel-scores references/revel_v1.3.tsv \
  --spliceai-scores references/spliceai_scores.tsv \
  --gene-list references/cardiac_panel.txt
```

### DOCX output

For a Word document compatible with institutional templates:

```bash
pip install python-docx

vartriage \
  --vcf proband_exome.vcf.gz \
  --output results/clinical_report.docx \
  --output-format clinical-docx \
  --patient-id PAT-2025-001 \
  --panel-name "Cardiac Panel v3" \
  --gene-annotation references/gencode.v44.gtf \
  --gnomad references/gnomad.v4.exomes.tsv \
  --clinvar references/clinvar_20240101.tsv \
  --cadd-scores references/cadd_v1.7.tsv \
  --revel-scores references/revel_v1.3.tsv \
  --spliceai-scores references/spliceai_scores.tsv \
  --gene-list references/cardiac_panel.txt
```

### Inspecting the audit trail

The `.audit.json` sidecar captures everything needed to reproduce the analysis:

```python
import json
from pathlib import Path

audit = json.loads(
    Path("results/clinical_report.html.audit.json").read_text()
)

# Run manifest: config, references, timestamps
manifest = audit["run_manifest"]
print(f"Patient: {manifest['patient_id']}")
print(f"Pipeline version: {manifest['pipeline_version']}")
print(f"Execution: {manifest['execution_timestamp']}")

# Decision log: one entry per classified variant
for entry in audit["decision_log"]:
    print(
        f"  {entry['gene_name']} {entry['chromosome']}:{entry['position']} "
        f"-> {entry['classification']} "
        f"(tags: {entry['evidence_tags_assigned']})"
    )
```
