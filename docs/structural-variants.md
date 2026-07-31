# Structural Variant Triage

Parse, annotate, score, and classify structural variants (DEL, DUP, INV, INS, BND, CNV) using the ClinGen 2020 technical standards for CNV interpretation.

## Quick start

```bash
vartriage sv --sv-vcf sv_calls.vcf.gz --output sv_report.json \
  --gene-annotation gencode.v46.gtf \
  --dosage-sensitivity clingen_dosage.tsv \
  --gnomad-sv gnomad_sv.bed \
  --pathogenic-regions clingen_pathogenic_regions.bed \
  --benign-regions clingen_benign_regions.bed
```

## Pipeline stages

```text
SVParser → SVAnnotator → SVScorer → SVClassifier → Report
```

1. **SVParser** - Streams SV records from VCF, extracting SVTYPE, END/SVLEN, confidence intervals, BND notation. Supports Manta, DELLY, GATK-SV, GRIDSS, and LUMPY INFO field conventions.

2. **SVAnnotator** - Determines gene overlap (whole-gene vs partial), attaches ClinGen dosage sensitivity (HI/TS scores), and matches against gnomAD-SV population frequencies via reciprocal overlap.

3. **SVScorer** - Computes a composite pathogenicity score (0-1) from four weighted components:
   - Gene impact severity (35%)
   - Dosage sensitivity (30%)
   - Population frequency rarity (20%)
   - SV size (15%)

4. **SVClassifier** - Applies ClinGen 2020 evidence framework (Sections 1-4) to produce a 5-tier classification: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign.

## CLI reference

### Standalone SV triage

```bash
vartriage sv [OPTIONS]
```

| Flag | Default | Description |
| ------ | --------- | ------------- |
| `--sv-vcf` | required | Input VCF with SV calls |
| `--output` | required | Output report path |
| `--output-format` | json | `json` or `csv` |
| `--gene-annotation` | None | GTF/GFF for gene overlap |
| `--dosage-sensitivity` | None | ClinGen dosage TSV |
| `--gnomad-sv` | None | gnomAD-SV BED (chrom, start, end, type, af) |
| `--pathogenic-regions` | None | Known pathogenic CNV regions BED |
| `--benign-regions` | None | Known benign CNV regions BED |
| `--min-sv-size` | 50 | Minimum SV size (bp) |
| `--max-sv-size` | 0 | Maximum SV size (0 = no limit) |
| `--sv-types` | all | Comma-separated: DEL,DUP,INV,INS,BND,CNV |
| `--max-af` | 0.01 | Maximum population frequency |
| `--min-quality` | 20.0 | Minimum QUAL threshold |
| `--reciprocal-overlap` | 0.5 | Minimum reciprocal overlap for frequency match |
| `--whole-gene-threshold` | 0.8 | Fraction to classify as whole-gene event |
| `--include-benign` | false | Include Benign/Likely_Benign in output |

### Combined SNV + SV analysis

Add `--sv-vcf` to the main vartriage command to run both pipelines:

```bash
vartriage --vcf snv.vcf.gz --sv-vcf sv_calls.vcf.gz \
  --output report.json --gene-annotation gencode.v46.gtf \
  --gnomad gnomad.tsv
```

This runs the standard SNV pipeline and the SV triage pipeline, producing two output files: `report.json` (SNVs) and `report_sv.json` (SVs).

## Reference data

### ClinGen dosage sensitivity

Tab-separated file with columns: `gene_symbol`, `hi_score`, `ts_score`

Score scale (ClinGen curation levels):

- 0 = No evidence for dosage sensitivity
- 1 = Little evidence
- 2 = Emerging evidence
- 3 = Sufficient evidence (gene is established HI or TS)

A bundled subset ships with the package at `vartriage/data/clingen_dosage.tsv`. To download the full dataset:

```bash
python scripts/download_clingen_dosage.py --output vartriage/data/clingen_dosage.tsv
```

### gnomAD-SV frequency database

Tab-separated BED format: `chrom`, `start`, `end`, `sv_type`, `allele_frequency`

Available as a bundle:

```bash
vartriage bundle download --bundle gnomad-sv
```

### Pathogenic/benign region BED files

Four-column BED: `chrom`, `start`, `end`, `syndrome_name`

Bundled files include curated ClinGen regions for common microdeletion/microduplication syndromes (22q11.2, Williams-Beuren, Prader-Willi/Angelman, etc.).

## Classification framework

The classifier implements the ACMG/ClinGen Technical Standards (Riggs et al. 2020):

**Section 1** - Genomic content assessment

- 1A: Contains protein-coding genes
- 1B: Contains established HI/TS gene

**Section 2** - Overlap with established regions

- 2A: Complete overlap with pathogenic region
- 2B-2C: Partial overlap
- 2D-2F: Benign region overlap

**Section 3** - Gene-level evaluation (losses)

- 3A: HI gene fully contained
- 3B: HI gene partially deleted
- 3C: Breakpoint within gene

**Section 4** - Duplication-specific

- 4F: TS gene fully contained
- 4G: Gene disrupted by dup breakpoint
- 4H: Intragenic dup without disruption

Evidence points accumulate and map to classification:

- >= 0.99: Pathogenic
- 0.90 - 0.98: Likely Pathogenic
- -0.89 to 0.89: VUS
- -0.98 to -0.90: Likely Benign
- <= -0.99: Benign

## Python API

```python
from pathlib import Path
from vartriage.structural import SVTriagePipeline, SVTriageConfig

config = SVTriageConfig(
    vcf_path=Path("sv_calls.vcf.gz"),
    output_path=Path("sv_report.json"),
    gene_annotation_path=Path("gencode.v46.gtf"),
    dosage_sensitivity_path=Path("clingen_dosage.tsv"),
    gnomad_sv_path=Path("gnomad_sv.bed"),
    pathogenic_regions_path=Path("clingen_pathogenic_regions.bed"),
    benign_regions_path=Path("clingen_benign_regions.bed"),
)

pipeline = SVTriagePipeline(config)
output = pipeline.run()
```

### Using individual components

```python
from vartriage.structural import SVParser, SVAnnotator, SVScorer, SVClassifier

# Parse
with SVParser(Path("sv_calls.vcf.gz"), min_size=50) as parser:
    for sv in parser:
        print(sv.chrom, sv.start, sv.end, sv.sv_type.value)

# Annotate
annotator = SVAnnotator(
    gene_annotation_path=Path("gencode.v46.gtf"),
    dosage_sensitivity_path=Path("clingen_dosage.tsv"),
)
annotated_stream = annotator.annotate(iter(parser))

# Score
scorer = SVScorer(max_allele_frequency=0.01)
scored_stream = scorer.score(annotated_stream)

# Classify
classifier = SVClassifier(
    pathogenic_regions=[("chr22", 18916842, 21465659)],
)
for classified in classifier.classify(scored_stream):
    print(classified.classification.value, classified.syndrome_name)
```

## Output format

### JSON

```json
{
  "pipeline": "structural_variant_triage",
  "version": "0.13.0",
  "total_variants": 3,
  "variants": [
    {
      "chrom": "chr22",
      "start": 18916842,
      "end": 21465659,
      "sv_type": "DEL",
      "length": 2548818,
      "consequence": "Whole_Gene_Deletion",
      "classification": "Pathogenic",
      "pathogenicity_score": 0.92,
      "evidence_score": 1.35,
      "evidence_categories": ["1A", "1B", "2A", "3A"],
      "genes_affected": 5,
      "hi_genes_affected": 1,
      "syndrome_name": "22q11.2 deletion syndrome (DiGeorge)",
      "gene_details": [...],
      "scoring": {"dosage_score": 1.0, "size_score": 0.9, "frequency_score": 1.0}
    }
  ]
}
```

### CSV

Flat format with columns: chrom, start, end, sv_type, length, consequence, classification, pathogenicity_score, evidence_score, genes_affected, hi_genes_affected, population_frequency, evidence_categories, gene_symbols.

## Supported SV callers

The parser handles caller-specific INFO field conventions:

| Caller | SVTYPE | END field | Length field | Copy number | Mate ID |
| -------- | -------- | ----------- | ------------- | ------------- | --------- |
| Manta | SVTYPE | END | SVLEN, INSLEN | - | MATEID |
| DELLY | SVTYPE | END | SVLEN | - | MATEID |
| GATK-SV | SVTYPE | END | SVLEN | CN, CNVAL | - |
| GRIDSS | SVTYPE | END, END2 | SVLEN | - | PARID |
| LUMPY | SVTYPE | END | SVLEN | - | MATEID |

All standard VCF 4.3 SV fields are supported. The parser auto-detects the available fields and resolves coordinates accordingly.
