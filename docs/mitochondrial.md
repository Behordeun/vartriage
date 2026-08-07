# Mitochondrial Variant Analysis

Classify mitochondrial DNA (chrM/MT) variants using mtDNA-specific criteria distinct from the nuclear ACMG/AMP 2015 framework. Automatic detection, no special flags required.

## Quick start

```bash
# mtDNA variants are automatically detected from chrM/MT in the VCF
vartriage --vcf wgs.vcf.gz --output results.json

# Override the minimum heteroplasmy threshold (default: 1%)
vartriage --vcf wgs.vcf.gz --output results.json --mt-min-heteroplasmy 5.0

# Skip mitochondrial analysis (targeted panels without mtDNA capture)
vartriage --vcf panel.vcf.gz --output results.json --skip-mito
```

## How it works

When the pipeline encounters chrM/MT variants in the input VCF, it routes them to a separate mitochondrial sub-pipeline while nuclear variants proceed through the standard ACMG flow. Results from both pipelines merge at the report stage.

```text
VCFParser
  ├── chrM/MT variants → MitochondrialPipeline
  │     HeteroplasmyExtractor → GeneMapAnnotator → MITOMAP/HelixMTdb → Classifier
  │
  └── Nuclear variants → Standard ACMG Pipeline
        QualityFilter → AnnotationEngine → Prioritization → ACMGClassifier

Both streams → ReportGenerator (merged output)
```

## Mitochondrial genetic code

The vertebrate mitochondrial genome uses a different translation table at 4 codons:

| Codon | Standard Code | Mitochondrial Code |
|-------|--------------|-------------------|
| TGA   | Stop (*)     | Trp (W)           |
| ATA   | Ile (I)      | Met (M)           |
| AGA   | Arg (R)      | Stop (*)          |
| AGG   | Arg (R)      | Stop (*)          |

VarTriage automatically selects the mitochondrial translation table when resolving amino acid changes for chrM/MT variants. This prevents false stop-gain calls (e.g., a TGA codon in MT-ATP6 is Trp, not a premature stop).

## Heteroplasmy

Mitochondrial DNA exists in hundreds to thousands of copies per cell. A variant can be present in a fraction of copies rather than the binary het/hom states of nuclear diploid variants.

### Extraction

Heteroplasmy is extracted from VCF FORMAT fields:

1. **AD field** (primary): `ALT_depth / (REF_depth + ALT_depth)`
2. **AF field** (fallback): Used directly when AD is unavailable (e.g., Mutect2 mitochondrial mode)

### Level classification

| Category       | Range      | Clinical relevance                          |
|---------------|------------|---------------------------------------------|
| Homoplasmic   | >= 95%     | Functionally equivalent to fixed mutation   |
| High          | 60-95%     | Strong phenotypic effect expected           |
| Moderate      | 20-60%     | Variable expressivity, tissue-dependent     |
| Low           | 1-20%      | Usually subclinical, monitor over time      |
| Sub-threshold | < 1%       | Likely sequencing noise, filtered by default|

### Filtering

Variants below the `--mt-min-heteroplasmy` threshold (default 1.0%) are excluded from output. This filters sequencing artifacts and ultra-low-level noise.

## Classification criteria

mtDNA variants are classified using a rule-based system independent of ACMG/AMP 2015:

### Pathogenic

All three conditions must be met:

- Confirmed pathogenic in MITOMAP (status = "Cfrm")
- High or homoplasmic heteroplasmy (>= 60%)
- Rare in HelixMTdb (AF < 0.01%)

### Likely Pathogenic

- Reported in MITOMAP (any status) AND
- Located in a protein-coding gene or tRNA AND
- Moderate or higher heteroplasmy (>= 20%)

Also assigned when a variant is confirmed in MITOMAP and rare, but heteroplasmy data is unavailable.

### Benign

- Common haplogroup-defining polymorphism: AF > 5% in HelixMTdb

### Likely Benign

- Moderate population frequency: AF > 0.1% in HelixMTdb

### VUS (Variant of Uncertain Significance)

- Default for variants without strong evidence in either direction
- Includes novel variants absent from both MITOMAP and HelixMTdb

## Data sources

### MITOMAP

[MITOMAP](https://www.mitomap.org/) is the primary mtDNA disease database. The bundled dataset includes ~60 confirmed and reported pathogenic mutations covering:

- MELAS (m.3243A>G, m.3271T>C)
- MERRF (m.8344A>G, m.8356T>C)
- LHON (m.3460G>A, m.11778G>A, m.14484T>C)
- NARP/Leigh syndrome (m.8993T>G/C)
- Aminoglycoside-induced deafness (m.1555A>G)

Update the bundled data:

```bash
python scripts/download_mitomap.py --output vartriage/data/mito/mitomap_pathogenic.tsv
```

### HelixMTdb

Population allele frequencies from [HelixMTdb](https://www.helix.com/pages/mitochondrial-variant-database), a large-scale mitochondrial database (~200k samples). Used to distinguish common haplogroup markers from rare potentially pathogenic variants.

Update the bundled data:

```bash
python scripts/download_helixmtdb.py --output vartriage/data/mito/helixmtdb_frequency.tsv
```

### MT gene map

The bundled gene map covers all 37 mitochondrial genes:

- 13 protein-coding genes (MT-ND1 through MT-ND6, MT-CO1-3, MT-ATP6, MT-ATP8, MT-CYB)
- 22 tRNA genes (MT-TF, MT-TL1, MT-TL2, MT-TK, etc.)
- 2 rRNA genes (MT-RNR1, MT-RNR2)
- D-loop/control region (positions 1-576 and 16024-16569)

Coordinates follow the revised Cambridge Reference Sequence (rCRS, NC_012920).

## Maternal inheritance check

When trio data is available (proband + mother + father), the pipeline verifies maternal inheritance:

- **Maternal**: Present in mother, absent in father (expected for mtDNA)
- **De novo**: Absent in both parents (rare but clinically significant)
- **Paternal unexpected**: Detected in father (indicates NUMTs, contamination, or data error)

This check is automatic when inheritance analysis is configured with `--proband`, `--mother`, and `--father` flags.

## Output format

### JSON output

When mitochondrial variants are present, JSON output includes a `mitochondrial_findings` array:

```json
{
  "variants": [...],
  "mitochondrial_findings": [
    {
      "chromosome": "chrM",
      "position": 3243,
      "ref_allele": "A",
      "alt_allele": "G",
      "mt_classification": "Likely_Pathogenic",
      "heteroplasmy_level": 85.2,
      "heteroplasmy_category": "high",
      "heteroplasmy_depth": 1200,
      "mitomap_disease": "MELAS / Leigh syndrome / Diabetes-Deafness",
      "mitomap_status": "Cfrm",
      "mt_gene": "MT-TL1",
      "mt_gene_type": "tRNA",
      "helix_af": 0.0002,
      "classification_reason": "Reported in MITOMAP ..."
    }
  ],
  "metadata": {
    "mitochondrial_note": "Mitochondrial variants classified using mtDNA-specific criteria"
  }
}
```

### CSV output

CSV output appends a mitochondrial findings section after the nuclear variants with dedicated columns:

- `mt_classification`, `heteroplasmy_level`, `heteroplasmy_category`
- `mitomap_disease`, `mitomap_status`
- `mt_gene`, `mt_gene_type`, `helix_af`

### Clinical reports

Clinical reports (HTML/PDF/DOCX) include a "Mitochondrial Findings" section with a table showing position, nucleotide change, gene, classification, heteroplasmy level, and MITOMAP disease association.

## CLI reference

| Flag | Default | Description |
|------|---------|-------------|
| `--skip-mito` | disabled | Skip mitochondrial analysis entirely |
| `--mt-min-heteroplasmy` | 1.0 | Minimum heteroplasmy % for reporting (0.0-100.0) |

## Module structure

```
vartriage/mito/
├── __init__.py          # Package marker
├── genetic_code.py      # MT codon table + is_mitochondrial() helper
├── gene_map.py          # 37-gene interval lookup (binary search)
├── heteroplasmy.py      # AD/AF extraction + level classification
├── mitomap.py           # MITOMAP disease database (O(1) lookup)
├── frequency.py         # HelixMTdb population frequency (O(1) lookup)
├── classifier.py        # 5-rule classification decision tree
├── maternal.py          # Trio maternal inheritance verification
├── pipeline.py          # Sub-pipeline orchestrator
├── config.py            # MitoConfig dataclass
└── report.py            # Output serialization for JSON/CSV

vartriage/data/mito/
├── mt_gene_map.tsv          # Gene intervals (rCRS coordinates)
├── mitomap_pathogenic.tsv   # Curated pathogenic mutations
└── helixmtdb_frequency.tsv  # Population allele frequencies
```
