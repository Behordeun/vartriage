# Gene-Disease Linkage Guide

Added in v0.12.0. Connects variants to their clinical context via gene-level annotations from OMIM, ClinGen, HPO, and gnomAD constraint data.

## Overview

When gene-disease linkage is active, each variant receives a `gene_context` containing:

- **Disease associations** from OMIM (disease name, MIM number, inheritance mode)
- **ClinGen validity** level (Definitive, Strong, Moderate, Limited, etc.)
- **gnomAD constraint** metrics (pLI, LOEUF, mis_z)
- **Actionability** status from ClinGen curations
- **Phenotype match score** (0.0-1.0) when patient HPO terms are provided

## Quick start

```bash
# Phenotype-driven prioritization for an epilepsy patient
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --hpo-terms HP:0001250,HP:0001249,HP:0002197

# Filter to autosomal recessive genes only
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --inheritance-mode AR

# Show only medically actionable findings
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --flag-actionable

# Custom knowledge data directory
vartriage --vcf patient.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --hpo-terms HP:0001250 --knowledge-dir /path/to/custom/data/
```

## CLI flags

| Flag                | Description                                                                 |
| ------------------- | --------------------------------------------------------------------------- |
| `--hpo-terms`       | Comma-separated HPO term IDs (HP:NNNNNNN format). Activates phenotype boost |
| `--inheritance-mode`| Filter to genes matching this mode: AD, AR, XL, XLD, XLR, MT                |
| `--flag-actionable` | Only pass variants in ClinGen actionable genes                              |
| `--knowledge-dir`   | Path to custom TSV data directory (defaults to bundled package data)        |

Any of these flags activates the gene-disease linkage pipeline stage.

## Phenotype-driven prioritization

When `--hpo-terms` is provided, variants in genes whose HPO annotations overlap with the patient's phenotype receive a ranking boost:

```text
boosted_score = prioritization_score * (1 + overlap)
```

Where `overlap` = (number of patient HPO terms matching the gene) / (total patient HPO terms).

The boost factor ranges from 1.0 (no match) to 2.0 (perfect match). This happens after prioritization scoring but before ACMG classification, so classification tiers are never violated by phenotype matching.

### Example

Patient presenting with seizures (HP:0001250) and intellectual disability (HP:0001249):

- **SCN1A** (Dravet syndrome): gene has both terms in HPO annotations. Overlap = 2/2 = 1.0. Score boosted by 2x.
- **BRCA1** (breast cancer): gene HPO terms don't overlap with epilepsy. Overlap = 0/2 = 0.0. No boost.

Result: SCN1A variants surface above BRCA1 in the prioritized output within the same classification tier.

## Inheritance mode filtering

When `--inheritance-mode AR` is set:

- Variants in genes with **only AD** inheritance in OMIM are dropped
- Variants in genes with **AR** (among potentially multiple modes) pass through
- Intergenic variants (no gene) pass through
- Genes **absent from OMIM** pass through (benefit of the doubt)

This is useful for consanguineous families where you expect homozygous AR variants to be causative.

## Actionability filtering

When `--flag-actionable` is set:

- Only variants in genes with ClinGen actionability curations pass through
- Intergenic variants pass through
- Genes without actionability data are filtered out

This surfaces findings where established medical interventions exist (surveillance, prophylactic surgery, therapeutic options).

## Output fields

### JSON

```json
{
  "chromosome": "chr2",
  "position": 166848884,
  "gene_name": "SCN1A",
  "acmg_classification": "Pathogenic",
  "disease_associations": [
    {
      "disease_name": "Dravet syndrome",
      "mim_number": "607208",
      "inheritance_mode": "AD"
    }
  ],
  "clingen_validity": "Definitive",
  "gene_constraint": {
    "pli": 1.0,
    "loeuf": 0.07,
    "mis_z": 5.44
  },
  "is_actionable": false,
  "phenotype_match_score": 1.0
}
```

### CSV

New columns appended after `evidence_tags`:

| Column                  | Description                                               |
| ----------------------- | --------------------------------------------------------- |
| `disease_associations`  | Semicolon-separated: `"Dravet syndrome [MIM:607208] (AD)"`|
| `clingen_validity`      | Definitive, Strong, Moderate, Limited, or empty           |
| `gene_constraint_pli`   | pLI score (0.0-1.0)                                       |
| `gene_constraint_loeuf` | LOEUF score                                               |
| `gene_constraint_mis_z` | Missense Z-score                                          |
| `is_actionable`         | True/False                                                |
| `phenotype_match_score` | Overlap fraction (0.0-1.0)                                |

## Python API

```python
from pathlib import Path
from vartriage import Pipeline, PipelineConfig, AnnotationConfig, ReportConfig
from vartriage.knowledge.config import KnowledgeBaseConfig

knowledge = KnowledgeBaseConfig(
    hpo_terms=frozenset({"HP:0001250", "HP:0001249"}),
    inheritance_mode="AD",
)

config = PipelineConfig(
    vcf_path=Path("patient.vcf.gz"),
    output_path=Path("results.json"),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("gencode.gtf"),
        gnomad_path=Path("gnomad.tsv"),
    ),
    report=ReportConfig(output_format="json"),
    knowledge=knowledge,
)

pipeline = Pipeline(config)
pipeline.run()
```

Access the registry directly for gene lookups:

```python
from vartriage.knowledge import GeneKnowledgeRegistry, KnowledgeBaseConfig

config = KnowledgeBaseConfig(
    hpo_terms=frozenset({"HP:0001250", "HP:0001249"}),
)
registry = GeneKnowledgeRegistry(config)

# Gene-level annotations
ann = registry.annotate_gene("SCN1A")
print(ann.disease_associations)   # tuple of DiseaseAssociation
print(ann.constraint.pli)         # 1.0
print(ann.clingen_validity)       # "Definitive"

# Phenotype overlap
overlap = registry.phenotype_overlap("SCN1A")  # 1.0
```

## Data files

Bundled TSV files live at `vartriage/data/knowledge/`:

| File                          | Columns                                                  | Content                        |
| ----------------------------- | -------------------------------------------------------- | ------------------------------ |
| `omim_gene_disease.tsv`       | gene_symbol, disease_name, mim_number, inheritance_mode  | OMIM gene-disease map          |
| `hpo_gene_annotations.tsv`    | gene_symbol, hpo_terms (semicolon-separated)             | HPO phenotype terms            |
| `clingen_validity.tsv`        | gene_symbol, validity_level                              | ClinGen validity               |
| `gnomad_constraint.tsv`       | gene_symbol, pli, loeuf, mis_z                           | gnomAD constraint              |
| `clingen_actionability.tsv`   | gene_symbol, intervention_type                           | ClinGen actionability          |

### Custom data

To use your own gene-disease data, create a directory with the same TSV format and pass it via `--knowledge-dir`. Missing files are handled gracefully (logged warning, empty lookups).

### Updating bundled data

The bundled dataset covers 22 clinically relevant genes and is current as of the release date. For broader coverage, prepare your own TSVs from:

- [ClinGen Gene-Disease Validity](https://search.clinicalgenome.org/kb/gene-validity)
- [HPO Gene Annotations](https://hpo.jax.org/data/annotations)
- [gnomAD Gene Constraint](https://gnomad.broadinstitute.org/downloads)
- [ClinGen Actionability](https://actionability.clinicalgenome.org/)

## Pipeline integration

```text
AnnotationEngine (consequence + frequency + ClinVar)
  → [GeneFilter]
  → GeneKnowledgeAnnotator (gene_context + inheritance filter + actionability filter)
  → PrioritizationEngine (frequency gate + scoring)
  → Phenotype Boost (score * (1 + overlap))
  → ACMGClassifier
  → ReportGenerator
```

The `GeneKnowledgeAnnotator` runs after the `AnnotationEngine` because it needs the resolved `gene_name` from consequence calling. It's constructed once at `Pipeline.__init__` and reused across runs to avoid re-loading TSV files.

## Validation

The `scripts/validate_affected_patient.py` script demonstrates a full clinical validation by:

1. Taking the healthy GIAB HG002 chr22 VCF (50,284 variants)
2. Spiking in 5 known pathogenic NF2 variants from ClinVar
3. Running the pipeline with NF2-matching HPO terms
4. Verifying that NF2 variants are correctly annotated with disease associations, constraint metrics, actionability, and phenotype match score

Run it:

```bash
python scripts/validate_affected_patient.py
```

## Limitations

- Phenotype matching uses exact HPO term overlap. Semantic similarity (HPO ontology traversal) is not yet implemented.
- Bundled data covers 22 genes. Broader coverage requires custom TSV files.
- Gene-disease data is static per release. Use `--knowledge-dir` for more current data.
- Inheritance mode filtering uses OMIM annotations only (not patient genotype).
