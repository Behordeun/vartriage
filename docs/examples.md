# Examples

Real-world usage patterns beyond the basic tutorial.

## Whole-genome QC run (no annotation)

Parse and export without any reference files. Good for validating a VCF before committing to a full annotated run.

```python
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, ReportConfig

config = PipelineConfig(
    vcf_path=Path("wgs_sample.vcf.gz"),
    output_path=Path("wgs_qc.json"),
    report=ReportConfig(output_format="json"),
)

pipeline = Pipeline(config)
pipeline.run()

# Check how many variants made it through quality filtering
import json
results = json.loads(Path("wgs_qc.json").read_text())
print(f"{len(results)} variants passed QC")
```

No annotation means no frequency lookups, no ClinVar hits — just parsing and quality filtering. Gives you variant counts and confirms the file is well-formed.

## Exome panel with strict AF filtering

For rare disease analysis, tighten the allele frequency cutoff to exclude anything seen in more than 1 in 1000 people:

```python
from pathlib import Path
from vartriage import (
    Pipeline, PipelineConfig, AnnotationConfig,
    PrioritizationConfig, QualityFilterConfig, ReportConfig,
)

config = PipelineConfig(
    vcf_path=Path("patient_exome.vcf.gz"),
    output_path=Path("rare_candidates.json"),
    quality_filter=QualityFilterConfig(min_qual=30.0),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("refs/gencode.v44.gtf"),
        gnomad_path=Path("refs/gnomad.v4.exomes.tsv"),
        clinvar_path=Path("refs/clinvar.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.001,  # 0.1% — strict rare disease threshold
        cadd_scores_path=Path("refs/cadd_v1.7.tsv"),
        revel_scores_path=Path("refs/revel_v1.3.tsv"),
    ),
    report=ReportConfig(output_format="json"),
)

pipeline = Pipeline(config)
pipeline.run()
```

Anything with gnomAD AF > 0.001 gets deprioritized. Combined with CADD and REVEL scores, this surfaces genuinely rare, computationally predicted-damaging variants.

## CLI batch processing

Loop over a directory of VCFs — one output per sample:

```bash
for vcf in samples/*.vcf.gz; do
  sample=$(basename "$vcf" .vcf.gz)
  vartriage \
    --vcf "$vcf" \
    --output "results/${sample}.json" \
    --gene-annotation refs/gencode.v44.gtf \
    --gnomad refs/gnomad.v4.exomes.tsv \
    --clinvar refs/clinvar.tsv \
    --cadd-scores refs/cadd_v1.7.tsv \
    --revel-scores refs/revel_v1.3.tsv
  echo "Done: $sample"
done
```

Each run is independent — safe to parallelize with `xargs` or GNU `parallel` if your machine has the RAM for it.

## Accessing individual stages

You don't have to run the full pipeline. Each stage is a standalone component.

### VCFParser + QualityFilter only

```python
from pathlib import Path
from vartriage import VCFParser, QualityFilter, QualityFilterConfig

qf = QualityFilter(QualityFilterConfig(min_qual=30.0))

with VCFParser(Path("input.vcf.gz")) as parser:
    passing = list(qf.apply(iter(parser)))

print(f"{len(passing)} variants passed QUAL >= 30")
```

### AnnotationEngine standalone

```python
from pathlib import Path
from vartriage import (
    VCFParser, QualityFilter, QualityFilterConfig,
    AnnotationEngine, AnnotationConfig,
)

ann = AnnotationEngine(AnnotationConfig(
    gene_annotation_path=Path("refs/gencode.v44.gtf"),
    gnomad_path=Path("refs/gnomad.v4.exomes.tsv"),
    clinvar_path=Path("refs/clinvar.tsv"),
))

with VCFParser(Path("input.vcf.gz")) as parser:
    qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
    annotated = list(ann.annotate(qf.apply(iter(parser))))

# Now you have AnnotatedVariant objects with consequence, AF, ClinVar
for v in annotated[:5]:
    print(f"{v.variant.chromosome}:{v.variant.position} "
          f"{v.consequence.value} AF={v.allele_frequency}")
```

## Custom filtering on output

Load the JSON output from a previous run and filter down to actionable variants:

```python
import csv
import json
from pathlib import Path

results = json.loads(Path("results/candidates.json").read_text())

# Keep only Pathogenic and Likely Pathogenic
pathogenic = [
    v for v in results
    if v["acmg_classification"] in ("Pathogenic", "Likely_Pathogenic")
]

print(f"{len(pathogenic)} actionable variants out of {len(results)} total")

# Export to CSV for sharing with clinicians
fields = [
    "chromosome", "position", "ref_allele", "alt_allele",
    "functional_consequence", "allele_frequency",
    "composite_rank", "acmg_classification", "clinvar_assertion",
]

with open("results/pathogenic_only.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(pathogenic)
```

## Memory monitoring

Track peak RSS during a large run. Useful for verifying the streaming pipeline stays within bounds on your hardware:

```python
import os
import resource
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, ReportConfig

config = PipelineConfig(
    vcf_path=Path("large_wgs.vcf.gz"),
    output_path=Path("large_output.json"),
    report=ReportConfig(output_format="json"),
)

pipeline = Pipeline(config)
pipeline.run()

# Peak RSS in MB (Linux/macOS)
peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
# macOS reports bytes, Linux reports KB
if os.uname().sysname == "Darwin":
    peak_mb = peak_kb / (1024 * 1024)
else:
    peak_mb = peak_kb / 1024

print(f"Peak RSS: {peak_mb:.1f} MB")
```

The pipeline streams variants in batches, so RSS should stay well under 2GB even for 4M+ variant WGS files.

## Suppressing warnings

vartriage emits warnings when variants lack gnomAD frequency data or ClinVar assertions. For large runs where you expect many missing annotations, suppress them:

```python
import warnings
from vartriage import VarTriageWarning

# Silence all vartriage warnings
warnings.filterwarnings("ignore", category=VarTriageWarning)
```

Or target a specific warning subclass:

```python
import warnings
from vartriage import MissingDataWarning

# Only suppress missing-data warnings, keep everything else
warnings.filterwarnings("ignore", category=MissingDataWarning)
```

You can also check what was suppressed after the run via the warning accumulator:

```python
pipeline.run()
acc = pipeline.warning_accumulator
print(f"Suppressed {acc.total_count} missing-data events")
print(f"Sources: {acc.sources}")
```
