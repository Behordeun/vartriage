# Examples

Usage examples covering every vartriage capability, organized by version.

## Sample output files

| File | Description |
| ------ | ------------- |
| `sample_pipeline_output.json` | Standard JSON output with ACMG classification, evidence tags, and `prioritization_score` |
| `sample_pipeline_output.csv` | CSV equivalent (15 columns including `prioritization_score`) |
| `sample_clinical_report.html` | Self-contained clinical HTML report for a hereditary cancer panel |
| `sample_clinical_report.html.audit.json` | Audit trail sidecar for the clinical report |

All samples use synthetic variant data. Patient identifiers are fictional.

**Output fields:**

- `prioritization_score`: literature-validated ranking (REVEL for missense, SpliceAI for splice-adjacent, CADD Phred/99 capped at 1.0 for others). Recommended for triage.
- `composite_rank`: legacy weighted average. Deprecated, will be removed in v1.0.0.
- `evidence_tags`: ACMG criteria the classifier emits (PVS1, PS1, PM1, PM2, PM4, PM5, PP3, PP5, BA1, BS1, BP4, BP7, plus the strength-modulated PVS1_Strong, PP3_Moderate, BP4_Moderate). BS2 is defined but not emitted.
- `acmg_classification`: final 5-tier call (Pathogenic, Likely_Pathogenic, VUS, Likely_Benign, Benign)

---

## 1. Basic pipeline (v0.1.0)

```bash
pip install vartriage

vartriage --vcf sample.vcf.gz --output candidates.json
```

With annotation references for consequence and frequency:

```bash
vartriage \
  --vcf sample.vcf.gz \
  --output candidates.json \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv
```

Python API equivalent:

```python
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, AnnotationConfig, PrioritizationConfig

config = PipelineConfig(
    vcf_path=Path("sample.vcf.gz"),
    output_path=Path("candidates.json"),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
    ),
    prioritization=PrioritizationConfig(
        cadd_scores_path=Path("cadd_scores.tsv"),
        revel_scores_path=Path("revel_scores.tsv"),
    ),
)

pipeline = Pipeline(config)
pipeline.run()
```

---

## 2. Multi-sample VCF and region filtering (v0.3.0)

```bash
vartriage \
  --vcf multi_sample.vcf.gz \
  --sample PROBAND_01 \
  --min-gq 20 \
  --regions target_regions.bed \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --output panel_results.json
```

---

## 3. Trio inheritance analysis (v0.4.0)

```bash
vartriage \
  --vcf trio.vcf.gz \
  --proband CHILD --mother MOM --father DAD \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --spliceai-scores spliceai_scores.tsv \
  --output trio_results.json
```

---

## 4. Gene panel filtering and VCF output (v0.4.0)

```bash
vartriage \
  --vcf exome.vcf.gz \
  --gene-list cardiac_panel.txt \
  --output-format vcf \
  --output annotated.vcf.gz \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv \
  --spliceai-scores spliceai_scores.tsv
```

The output VCF carries `VARTRIAGE_CONSEQUENCE`, `VARTRIAGE_AF`,
`VARTRIAGE_RANK`, `VARTRIAGE_ACMG`, and `VARTRIAGE_TAGS` INFO fields.
Includes a tabix index for fast region queries.

---

## 5. Clinical report generation (v0.5.0)

```bash
vartriage \
  --vcf exome.vcf.gz \
  --output clinical_report.html \
  --output-format clinical-html \
  --patient-id "PAT-2026-001" \
  --panel-name "Hereditary Cancer Panel v3" \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv \
  --spliceai-scores spliceai_scores.tsv \
  --gene-list hereditary_cancer_panel.txt
```

Available formats: `clinical-html`, `clinical-pdf` (requires weasyprint),
`clinical-docx` (requires python-docx).

---

## 6. Score bundle management (v0.6.0)

```bash
# List available bundles
vartriage bundle list

# Download reference data
vartriage bundle download --bundle clinvar
vartriage bundle download --bundle gnomad-exomes-chr22
vartriage bundle download --bundle gencode
vartriage bundle download --bundle revel

# Run with auto-resolved paths
vartriage --vcf sample.vcf.gz --output results.json --use-bundles

# Check bundle status
vartriage bundle status
```

---

## 7. API mode (v0.7.0)

```bash
pip install vartriage[api]

# Pure API mode
vartriage --vcf panel.vcf --output results.json --mode api

# Hybrid: local gnomAD + remote ClinVar/CADD
vartriage --vcf panel.vcf --output results.json \
  --mode hybrid --gnomad gnomad.tsv

# With NCBI API key for faster ClinVar queries
vartriage --vcf panel.vcf --output results.json \
  --mode api --api-key YOUR_NCBI_KEY
```

Queries Ensembl VEP, ClinVar, CADD, and SpliceAI. Responses cached in
SQLite (`~/.vartriage/api_cache.db`).

---

## 8. Codon-level consequence calling (v0.8.0)

```bash
vartriage \
  --vcf exome.vcf.gz \
  --output results.json \
  --reference-fasta GRCh38.fa \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --use-bundles
```

Without `--reference-fasta`, the pipeline uses a positional heuristic
(any CDS SNV = Missense). With FASTA, it performs actual codon
translation and correctly calls synonymous, missense, and nonsense.

---

## 9. ACMG Secondary Findings screening (v0.10.0)

Screen against the 71-gene ACMG SF v3.2 list:

```bash
vartriage \
  --vcf wgs.vcf.gz \
  --output results.json \
  --secondary-findings \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --use-bundles
```

---

## 10. Multi-sample cohort analysis (v0.11.0)

```bash
# From a manifest file
vartriage cohort \
  --manifest samples.tsv \
  --output cohort_results/ \
  --cohort-name "cardiac_cohort" \
  --min-recurrence 2 \
  --max-af 0.01 \
  --parallel --max-workers 8 \
  --use-bundles

# Or pass VCFs directly
vartriage cohort \
  --vcf patient_001.vcf.gz patient_002.vcf.gz patient_003.vcf.gz \
  --output cohort_results/ \
  --cohort-name "trio_cohort" \
  --output-format csv
```

Produces three output files: variants (with recurrence counts),
gene burden table, and cohort summary.

Python API:

```python
from pathlib import Path
from vartriage import CohortPipeline, CohortConfig, AnnotationConfig

cohort_config = CohortConfig(
    sample_vcfs=[Path("p1.vcf.gz"), Path("p2.vcf.gz"), Path("p3.vcf.gz")],
    output_path=Path("cohort_results/"),
    cohort_name="cardiac_cohort",
    min_recurrence=2,
    parallel=True,
    max_workers=4,
)

pipeline = CohortPipeline(
    cohort_config=cohort_config,
    annotation_config=AnnotationConfig(
        gene_annotation_path=Path("gencode.v44.gtf"),
        gnomad_path=Path("gnomad.v4.sites.tsv"),
    ),
)
pipeline.run()
```

---

## 11. Gene-disease linkage and phenotype-driven prioritization (v0.12.0)

```bash
# Phenotype-driven: boost genes matching patient symptoms
vartriage \
  --vcf patient.vcf.gz \
  --output results.json \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --hpo-terms HP:0001250,HP:0001249,HP:0002197 \
  --use-bundles

# Filter to autosomal recessive genes only
vartriage \
  --vcf patient.vcf.gz \
  --output results.json \
  --inheritance-mode AR \
  --use-bundles

# Only ClinGen actionable findings
vartriage \
  --vcf patient.vcf.gz \
  --output results.json \
  --flag-actionable \
  --use-bundles

# Combined: phenotype + inheritance + actionability
vartriage \
  --vcf patient.vcf.gz \
  --output results.json \
  --hpo-terms HP:0001250,HP:0001249 \
  --inheritance-mode AD \
  --flag-actionable \
  --use-bundles
```

Output includes per-variant: disease associations (MIM numbers),
ClinGen validity, gnomAD constraint (pLI/LOEUF/mis_z), actionability,
and phenotype match score.

---

## 12. Structural variant triage (v0.13.0)

```bash
# Standalone SV analysis
vartriage sv \
  --sv-vcf sv_calls.vcf.gz \
  --output sv_report.json \
  --gene-annotation gencode.v46.gtf \
  --dosage-sensitivity clingen_dosage.tsv \
  --gnomad-sv gnomad_sv.bed \
  --pathogenic-regions clingen_pathogenic_regions.bed \
  --benign-regions clingen_benign_regions.bed

# Combined SNV + SV in one command
vartriage \
  --vcf snv.vcf.gz \
  --sv-vcf sv_calls.vcf.gz \
  --output report.json \
  --gene-annotation gencode.v46.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --use-bundles

# CSV output with benign SVs included
vartriage sv \
  --sv-vcf sv_calls.vcf.gz \
  --output sv_report.csv \
  --output-format csv \
  --include-benign \
  --gene-annotation gencode.v46.gtf \
  --dosage-sensitivity clingen_dosage.tsv
```

Python API:

```python
from pathlib import Path
from vartriage.structural import SVTriagePipeline, SVTriageConfig

config = SVTriageConfig(
    vcf_path=Path("sv_calls.vcf.gz"),
    output_path=Path("sv_report.json"),
    gene_annotation_path=Path("gencode.v46.gtf"),
    dosage_sensitivity_path=Path("clingen_dosage.tsv"),
    gnomad_sv_path=Path("gnomad_sv.bed"),
    pathogenic_regions_path=Path("pathogenic_regions.bed"),
    benign_regions_path=Path("benign_regions.bed"),
)

pipeline = SVTriagePipeline(config)
output = pipeline.run()
```

Supports Manta, DELLY, GATK-SV, GRIDSS, and LUMPY. See
[Structural Variants Guide](../structural-variants.md) for full
CLI reference and output format details.

---

## 13. Mitochondrial variant analysis (v0.15.0)

```bash
# Automatic: chrM/MT variants detected and classified separately
vartriage --vcf wgs.vcf.gz --output results.json --use-bundles

# Custom heteroplasmy threshold
vartriage --vcf wgs.vcf.gz --output results.json --mt-min-heteroplasmy 5.0

# Skip mitochondrial analysis (targeted panels without mtDNA capture)
vartriage --vcf panel.vcf.gz --output results.json --skip-mito
```

The mitochondrial pipeline uses the vertebrate mitochondrial genetic code,
extracts heteroplasmy from AD/AF fields, queries MITOMAP for disease
associations, checks HelixMTdb for population frequency, and classifies
independently of nuclear ACMG criteria. Results appear in a dedicated
"Mitochondrial Findings" section.

---

## 14. Remote tabix scoring (v0.16.0)

```bash
# CADD scores via HTTP byte-range (no 80 GB download)
vartriage --vcf panel.vcf --output results.json \
  --gene-annotation gencode.gtf --cadd-remote cadd-v1.7-grch38

# gnomAD frequencies via remote tabix
vartriage --vcf panel.vcf --output results.json \
  --gene-annotation gencode.gtf --gnomad-remote gnomad-exomes-v4-grch38

# Both remote, with pinned cache for clinical reproducibility
vartriage --vcf panel.vcf --output results.json \
  --gene-annotation gencode.gtf \
  --cadd-remote cadd-v1.7-grch38 \
  --gnomad-remote gnomad-exomes-v4-grch38 \
  --remote-cache-ttl -1

# List available presets
vartriage remote list-presets
```

Scores are cached in SQLite (`~/.vartriage/remote_cache.db`, 30-day TTL).
Local files always take priority. A circuit breaker prevents stalling on
network failures.

---

## 15. SpliceAI SQLite backend (v0.17.5)

Query precomputed SpliceAI delta scores directly from the OpenCRAVAT SQLite
database instead of pre-filtering scores into a per-analysis TSV:

```bash
vartriage \
  --vcf exome.vcf.gz \
  --output results.json \
  --gene-annotation gencode.v46.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --spliceai-db spliceai_scores.db
```

`--spliceai-db` (SQLite) and `--spliceai-scores` (TSV) are mutually
exclusive. The database backend returns the max of the four delta scores
(ds_ag, ds_al, ds_dg, ds_dl) per variant and covers all chromosomes without
a pre-filtering step.

---

## 16. VCF quality control (v0.18.0)

```bash
# Standalone QC check, no annotation
vartriage qc --vcf sample.vcf.gz --sample SAMPLE1 --assay-type wes

# Write a machine-readable QC report
vartriage qc --vcf sample.vcf.gz --sample SAMPLE1 --assay-type wgs \
  --output-json qc_report.json

# Gate the full pipeline: halt before annotation on any FAIL
vartriage --vcf sample.vcf.gz --output results.json \
  --gene-annotation gencode.v46.gtf --gnomad gnomad.v4.sites.tsv \
  --assay-type wgs --strict-qc

# Skip QC for a pre-validated file
vartriage --vcf sample.vcf.gz --output results.json \
  --gene-annotation gencode.v46.gtf --gnomad gnomad.v4.sites.tsv --skip-qc
```

QC computes Ti/Tv, het/hom, variant count, and ins/del ratios in a single
streaming pass, validates them against the assay-specific ranges (`wgs`,
`wes`, `panel`), and prints a PASS/WARN/FAIL table to stderr. With
`--strict-qc`, a FAIL halts the pipeline before annotation and exits with
code 3. See [Quality Control Guide](../quality-control.md) for metric
definitions and thresholds.

Python API:

```python
from pathlib import Path
from vartriage.qc.config import QCConfig
from vartriage.qc.metrics import compute_qc_metrics
from vartriage.qc.validator import QCValidator

metrics = compute_qc_metrics(Path("sample.vcf.gz"), sample_id="SAMPLE1")
report = QCValidator(QCConfig(assay_type="wes")).validate(metrics)
print(report.overall_status)
```

---

## Complete pipeline example

```bash
vartriage \
  --vcf patient_exome.vcf.gz \
  --sv-vcf patient_sv.vcf.gz \
  --output clinical_report.html \
  --output-format clinical-html \
  --patient-id "PAT-2026-042" \
  --panel-name "Comprehensive Genomics Panel" \
  --reference-fasta GRCh38.fa \
  --gene-annotation gencode.v46.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv \
  --cadd-scores cadd_scores.tsv \
  --revel-scores revel_scores.tsv \
  --spliceai-scores spliceai_scores.tsv \
  --gene-list panel_genes.txt \
  --hpo-terms HP:0001250,HP:0001249 \
  --secondary-findings \
  --assay-type wes \
  --strict-qc \
  --use-bundles
```
