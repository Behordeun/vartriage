# Reference Files

All reference files are tab-separated with a header row. Column names must match exactly.

## gnomAD

Population allele frequency data. Used by the annotation engine for frequency lookup and by the prioritization engine for the frequency gate.

**Required columns:**

| Column | Type | Description |
| -------- | ------ | ------------- |
| `chrom` | str | Chromosome name (e.g., "chr1" or "1") |
| `pos` | int | 1-based genomic position |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele |
| `af` | float | Global allele frequency (0.0 to 1.0) |

**Example:**

```tsv
chrom    pos      ref    alt    af
chr1     12345    A      G      0.00032
chr1     54321    CT     C      0.15
chr2     99999    G      T      0.0
```

**Where to download:**

[gnomAD Downloads](https://gnomad.broadinstitute.org/downloads). Use the sites VCF, then extract the required columns with bcftools or a script.

**Preparation from gnomAD VCF:**

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' \
    gnomad.genomes.v4.sites.vcf.bgz \
    | awk 'BEGIN{print "chrom\tpos\tref\talt\taf"}{print}' \
    > gnomad.v4.sites.tsv
```

## ClinVar

Clinical significance assertions. Used during annotation and ACMG classification (PP5 criterion).

**Required columns:**

| Column | Type | Description |
| -------- | ------ | ------------- |
| `chrom` | str | Chromosome name |
| `pos` | int | 1-based genomic position |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele |
| `clinical_significance` | str | One of the accepted values below |

**Accepted values for `clinical_significance`:**

- `Pathogenic`
- `Likely pathogenic`
- `Uncertain significance`
- `Likely benign`
- `Benign`

**Example:**

```tsv
chrom    pos          ref     alt    clinical_significance
chr1     12345        A       G      Pathogenic
chr7     117559590    ATCT    A      Likely pathogenic
chr17    7674220      G       A      Uncertain significance
```

**Where to download:**

[ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/). Use `variant_summary.txt.gz` and extract the relevant columns.

**Preparation from ClinVar variant_summary:**

```bash
zcat variant_summary.txt.gz \
    | awk -F'\t' 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt","clinical_significance"} \
    NR>1 && $17=="GRCh38" {print "chr"$19, $20, $22, $23, $7}' \
    > clinvar_20240101.tsv
```

## CADD

Combined Annotation Dependent Depletion scores. Phred-scaled scores where higher means more likely deleterious.

**Required columns:**

| Column | Type | Description |
| -------- | ------ | ------------- |
| `chrom` | str | Chromosome name |
| `pos` | int | 1-based genomic position |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele |
| `score` | float | CADD Phred score (typically 0 to 99) |

**Example:**

```tsv
chrom    pos      ref    alt    score
chr1     12345    A      G      28.5
chr1     54321    C      T      15.2
chr2     99999    G      A      35.0
```

**Where to download:**

[CADD Downloads](https://cadd.gs.washington.edu/download). Get the whole-genome or exome pre-scored file, then extract columns.

**Interpretation:**

- Score >= 20: top 1% most deleterious substitutions
- Score >= 30: top 0.1%
- The pipeline normalizes by dividing by 99 and capping at 1.0

## REVEL

Rare Exome Variant Ensemble Learner. An ensemble score for missense variants.

**Required columns:**

| Column | Type | Description |
| -------- | ------ | ------------- |
| `chrom` | str | Chromosome name |
| `pos` | int | 1-based genomic position |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele |
| `score` | float | REVEL score (0.0 to 1.0) |

**Example:**

```tsv
chrom    pos          ref    alt    score
chr1     12345        A      G      0.85
chr1     54321        C      T      0.12
chr7     117559590    A      G      0.95
```

**Where to download:**

[REVEL Downloads](https://sites.google.com/site/revelgenomics/downloads). The distribution file uses `chr`, `hg19_pos`, `grch38_pos` columns; extract and reformat.

**Interpretation:**

- Score > 0.7: likely pathogenic (PP3 threshold used by this library)
- Score < 0.15: likely benign
- REVEL is already normalized to 0-1, no additional scaling applied

## SpliceAI

Deep-learning splicing impact predictions. Scores represent the maximum delta score across all splice-altering events (acceptor gain/loss, donor gain/loss).

**Required columns:**

| Column | Type | Description |
| -------- | ------ | ------------- |
| `chrom` | str | Chromosome name |
| `pos` | int | 1-based genomic position |
| `ref` | str | Reference allele |
| `alt` | str | Alternate allele |
| `score` | float | SpliceAI max delta score (0.0 to 1.0) |

**Example:**

```tsv
chrom    pos          ref    alt    score
chr1     12345        A      G      0.85
chr1     54321        C      T      0.02
chr7     117559590    A      G      0.15
```

**Where to download:**

[SpliceAI pre-computed scores](https://basespace.illumina.com/s/otSPW8hnhaZR) from Illumina. Extract the max delta score from the VCF INFO field and reformat to the TSV layout above.

**Interpretation:**

- Score > 0.8: high splicing impact (triggers PVS1 for splice-site variants)
- Score > 0.5: moderate splicing impact (triggers PP3 for splice-adjacent variants)
- Score > 0.2: suggests possible splicing effect
- Scores outside [0.0, 1.0] are rejected with a warning and treated as unavailable

## Chromosome naming

The library performs exact string matching on chromosome names. Your VCF and reference files must use the same naming convention.

Common conventions:

- UCSC style: `chr1`, `chr2`, ..., `chrX`, `chrY`
- Ensembl style: `1`, `2`, ..., `X`, `Y`

If your VCF uses `chr1` but your gnomAD file uses `1`, lookups will return no matches. Normalize before running the pipeline.
