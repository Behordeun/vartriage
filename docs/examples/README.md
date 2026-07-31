# Examples

Usage examples covering every vartriage capability, organized by version.

## Sample output files

| File | Description |
|------|-------------|
| `sample_pipeline_output.json` | Standard JSON output with ACMG classification and evidence tags |
| `sample_pipeline_output.csv` | CSV equivalent (19 columns) |
| `sample_clinical_report.html` | Self-contained clinical HTML report for a hereditary cancer panel |
| `sample_clinical_report.html.audit.json` | Audit trail sidecar for the clinical report |

All samples use synthetic variant data. Patient identifiers are fictional.

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
  --use-bundles
```
