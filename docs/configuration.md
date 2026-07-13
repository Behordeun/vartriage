# Configuration

All configuration classes are frozen dataclasses. They validate parameters at construction and raise `ValueError` for out-of-range values.

## QualityFilterConfig

Controls which variants pass the quality gate.

| Field | Type | Default | Valid range |
| ------- | ------ | --------- | ------------- |
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
| ------- | ------ | --------- | --------- | ------------------- |
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
| ------- | ------ | --------- | ------------- |
| `max_allele_frequency` | `float` | `0.01` | 0.0 to 1.0 |
| `cadd_scores_path` | `Optional[Path]` | `None` | CADD Phred TSV |
| `revel_scores_path` | `Optional[Path]` | `None` | REVEL scores TSV |
| `spliceai_scores_path` | `Optional[Path]` | `None` | SpliceAI scores TSV |
| `batch_size` | `int` | `10_000` | 1,000 to 100,000 |

```python
from vartriage import PrioritizationConfig

# Rare disease: very stringent frequency cutoff
config = PrioritizationConfig(
    max_allele_frequency=0.0001,
    cadd_scores_path=Path("cadd_v1.7.tsv"),
    revel_scores_path=Path("revel_v1.3.tsv"),
    spliceai_scores_path=Path("spliceai_scores.tsv"),
)

# Research: relaxed frequency, no score files
config = PrioritizationConfig(
    max_allele_frequency=0.05,
)
```

## ReportConfig

Output format selection.

| Field | Type | Default | Options |
| ------- | ------ | --------- | --------- |
| `output_format` | `Literal[...]` | `"json"` | `"json"`, `"csv"`, `"pdf"`, `"vcf"`, `"clinical-html"`, `"clinical-pdf"`, `"clinical-docx"` |

```python
from vartriage import ReportConfig

config = ReportConfig(output_format="csv")
```

When the format is `clinical-html`, `clinical-pdf`, or `clinical-docx`, the pipeline constructs a `ClinicalReportConfig` instead and delegates to the clinical report generator. See the ClinicalReportConfig section below.

## InheritanceConfig

Trio-based inheritance pattern classification settings.

| Field | Type | Default | Notes |
| ------- | ------ | --------- | ---- | ------- |
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

## ClinicalReportConfig

Configuration for structured clinical report generation. Required when `--output-format` is `clinical-html`, `clinical-pdf`, or `clinical-docx`.

| Field | Type | Default | Notes |
| ------- | ------ | --------- | ------- |
| `patient_id` | `str` | required | Patient identifier (non-empty) |
| `panel_name` | `str` | required | Gene panel name (non-empty) |
| `output_format` | `Literal` | required | `"clinical-pdf"`, `"clinical-html"`, or `"clinical-docx"` |
| `report_template` | `str` | `"standard"` | Report template name |

```python
from vartriage.models.config import ClinicalReportConfig

config = ClinicalReportConfig(
    patient_id="PAT-2025-001",
    panel_name="Cardiac Panel v3",
    output_format="clinical-html",
)

# With custom template
config = ClinicalReportConfig(
    patient_id="PAT-2025-002",
    panel_name="Hereditary Cancer Panel",
    output_format="clinical-pdf",
    report_template="standard",
)
```

Raises `ValueError` at construction if `patient_id` or `panel_name` is empty or whitespace-only.

The clinical report produces:

- Self-contained HTML (no external resources, no JavaScript)
- PDF via WeasyPrint (install with `pip install weasyprint`)
- DOCX via python-docx (install with `pip install python-docx`)

A JSON audit trail sidecar (`.audit.json`) is written alongside every clinical report. It contains the run manifest (config, reference checksums, timestamps) and a per-variant decision log.

## GeneFilterConfig

Restricts analysis to variants in a user-supplied gene list.

| Field | Type | Default | Notes |
| ------- | ------ | --------- | ------- |
| `gene_list_path` | `Path` | required | Plain text file, one gene symbol per line |

```python
from pathlib import Path
from vartriage import GeneFilterConfig

config = GeneFilterConfig(gene_list_path=Path("cardiac_panel.txt"))
```

The gene list file format: one symbol per line, blank lines and lines starting with `#` are skipped, matching is case-insensitive.

## MissingDataConfig

Controls the missing data warning threshold.

| Field | Type | Default | Notes |
| ------- | ------ | --------- | ------- |
| `warning_threshold` | `int` | `1000` | Summary warning fires when exceeded |

```python
from vartriage import MissingDataConfig

config = MissingDataConfig(warning_threshold=500)
```

## PipelineConfig

Top-level configuration aggregating all sub-configs.

| Field | Type | Default | Notes |
| ------- | ------ | --------- | ------- |
| `vcf_path` | `Path` | required | `.vcf` or `.vcf.gz` |
| `output_path` | `Path` | required | Output report path |
| `quality_filter` | `QualityFilterConfig` | default instance | |
| `annotation` | `Optional[AnnotationConfig]` | `None` | None skips annotation |
| `prioritization` | `PrioritizationConfig` | default instance | |
| `report` | `ReportConfig` | default instance | |
| `missing_data` | `MissingDataConfig` | default instance | |
| `inheritance` | `Optional[InheritanceConfig]` | `None` | None skips trio analysis |
| `gene_filter` | `Optional[GeneFilterConfig]` | `None` | None skips gene filtering |
| `region_filter` | `Optional[RegionFilterConfig]` | `None` | None skips region filtering |
| `sample` | `Optional[SampleConfig]` | `None` | None skips sample extraction |
| `clinical_report` | `Optional[ClinicalReportConfig]` | `None` | Required for clinical formats |

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

### Clinical report output

For sign-off-ready clinical reporting with audit trail:

```python
from vartriage.models.config import ClinicalReportConfig

config = PipelineConfig(
    vcf_path=Path("patient.vcf.gz"),
    output_path=Path("clinical_report.html"),
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
        spliceai_scores_path=Path("spliceai.tsv"),
    ),
    report=ReportConfig(output_format="clinical-html"),
    clinical_report=ClinicalReportConfig(
        patient_id="PAT-2025-001",
        panel_name="Cardiac Panel v3",
        output_format="clinical-html",
    ),
)
```

This writes `clinical_report.html` (self-contained, viewable offline) and `clinical_report.html.audit.json` (machine-parseable decision log).
