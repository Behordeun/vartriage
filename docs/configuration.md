# Configuration

All configuration classes are frozen dataclasses. They validate parameters at construction and raise `ValueError` for out-of-range values.

## QualityFilterConfig

Controls which variants pass the quality gate.

| Field | Type | Default | Valid range |
|-------|------|---------|-------------|
| `min_qual` | `float` | `20.0` | 0 to 1,000,000 |

```python
from vartriage import QualityFilterConfig

# Default: QUAL >= 20
config = QualityFilterConfig()

# Stringent clinical threshold
config = QualityFilterConfig(min_qual=50.0)
```

## AnnotationConfig

Paths to reference data and batch size for the annotation engine.

| Field | Type | Default | Valid range / notes |
|-------|------|---------|-------------------|
| `gene_annotation_path` | `Path` | required | GTF or GFF file |
| `gnomad_path` | `Path` | required | gnomAD TSV |
| `clinvar_path` | `Optional[Path]` | `None` | ClinVar TSV, or None to skip |
| `batch_size` | `int` | `10_000` | 1,000 to 100,000 |

```python
from pathlib import Path
from vartriage import AnnotationConfig

config = AnnotationConfig(
    gene_annotation_path=Path("gencode.v44.gtf"),
    gnomad_path=Path("gnomad.v4.sites.tsv"),
    clinvar_path=Path("clinvar_20240101.tsv"),
    batch_size=50_000,  # larger batches for faster processing
)
```

## PrioritizationConfig

Controls frequency filtering and pathogenicity scoring.

| Field | Type | Default | Valid range |
|-------|------|---------|-------------|
| `max_allele_frequency` | `float` | `0.01` | 0.0 to 1.0 |
| `cadd_scores_path` | `Optional[Path]` | `None` | CADD Phred TSV |
| `revel_scores_path` | `Optional[Path]` | `None` | REVEL scores TSV |
| `batch_size` | `int` | `10_000` | 1,000 to 100,000 |

```python
from vartriage import PrioritizationConfig

# Rare disease: very stringent frequency cutoff
config = PrioritizationConfig(
    max_allele_frequency=0.0001,
    cadd_scores_path=Path("cadd_v1.7.tsv"),
    revel_scores_path=Path("revel_v1.3.tsv"),
)

# Research: relaxed frequency, no score files
config = PrioritizationConfig(
    max_allele_frequency=0.05,
)
```

## ReportConfig

Output format selection.

| Field | Type | Default | Options |
|-------|------|---------|---------|
| `output_format` | `Literal["json", "csv", "pdf"]` | `"json"` | `"json"`, `"csv"`, `"pdf"` |

```python
from vartriage import ReportConfig

config = ReportConfig(output_format="csv")
```

## InheritanceConfig

Trio-based inheritance pattern classification settings.

| Field | Type | Default | Notes |
|-------|------|---------|-------|
| `proband` | `str` | required | Proband sample name |
| `mother` | `str` | required | Mother sample name |
| `father` | `str` | required | Father sample name |
| `patterns` | `list[str]` | all five | Patterns to evaluate |

Supported patterns: `de_novo`, `dominant`, `recessive`, `compound_het`, `x_linked`.

```python
from vartriage import InheritanceConfig

# All patterns (default)
config = InheritanceConfig(proband="CHILD", mother="MOM", father="DAD")

# Only de novo and recessive
config = InheritanceConfig(
    proband="CHILD", mother="MOM", father="DAD",
    patterns=["de_novo", "recessive"],
)
```

Raises `ValueError` if sample names are empty, patterns list is empty, or any pattern is not in the supported set.

## MissingDataConfig

Controls the missing data warning threshold.

| Field | Type | Default | Notes |
|-------|------|---------|-------|
| `warning_threshold` | `int` | `1000` | Summary warning fires when exceeded |

```python
from vartriage import MissingDataConfig

config = MissingDataConfig(warning_threshold=500)
```

## PipelineConfig

Top-level configuration aggregating all sub-configs.

| Field | Type | Default | Notes |
|-------|------|---------|-------|
| `vcf_path` | `Path` | required | `.vcf` or `.vcf.gz` |
| `output_path` | `Path` | required | Output report path |
| `quality_filter` | `QualityFilterConfig` | default instance | |
| `annotation` | `Optional[AnnotationConfig]` | `None` | None skips annotation |
| `prioritization` | `PrioritizationConfig` | default instance | |
| `report` | `ReportConfig` | default instance | |
| `missing_data` | `MissingDataConfig` | default instance | |
| `inheritance` | `Optional[InheritanceConfig]` | `None` | None skips trio analysis |

## Example configurations

### Stringent clinical filtering

For rare Mendelian disease panels:

```python
config = PipelineConfig(
    vcf_path=Path("patient.vcf.gz"),
    output_path=Path("clinical_report.json"),
    quality_filter=QualityFilterConfig(min_qual=50.0),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
        clinvar_path=Path("clinvar.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.0001,
        cadd_scores_path=Path("cadd.tsv"),
        revel_scores_path=Path("revel.tsv"),
    ),
    report=ReportConfig(output_format="pdf"),
)
```

### Relaxed research filtering

For exploratory variant discovery:

```python
config = PipelineConfig(
    vcf_path=Path("cohort_merged.vcf.gz"),
    output_path=Path("research_candidates.csv"),
    quality_filter=QualityFilterConfig(min_qual=10.0),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.05,
    ),
    report=ReportConfig(output_format="csv"),
)
```

### Annotation-only mode

Run annotation without full pipeline scoring. Use individual stages:

```python
from vartriage import VCFParser, QualityFilter, AnnotationEngine

with VCFParser(Path("input.vcf.gz")) as parser:
    qf = QualityFilter(QualityFilterConfig(min_qual=20.0))
    engine = AnnotationEngine(AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
    ))
    for annotated in engine.annotate(qf.apply(iter(parser))):
        print(annotated.consequence, annotated.allele_frequency)
```
