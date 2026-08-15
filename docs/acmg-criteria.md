# ACMG/AMP Variant Classification

How vartriage evaluates ACMG/AMP 2015 evidence criteria and combines them into a five-tier classification. This document covers every criterion the classifier implements, the data sources each one requires, ClinGen-calibrated thresholds, combining rules, and how to enable PS1/PM5 via the ClinVar protein index.

## Overview

The `ACMGClassifier` evaluates each scored variant against 10 evidence criteria (6 pathogenic, 4 benign) and applies ACMG/AMP 2015 combining rules to produce a final classification: Pathogenic, Likely Pathogenic, VUS, Likely Benign, or Benign. PP3 and BP4 each fire at two strength levels (supporting or moderate) based on ClinGen-calibrated thresholds, but they are single criteria with strength modulation, not separate criteria.

When a required data source is unavailable for a criterion, that criterion is skipped and the source is recorded in `missing_data_sources` on the output. This makes the system additive: providing more reference data enables more criteria without changing behavior for the criteria you already have.

## Pathogenic criteria

### PVS1 (Very Strong / Strong)

Null variant in a gene where loss-of-function is a known mechanism of disease.

| Condition | Fires when | Strength |
| --------- | ---------- | -------- |
| Nonsense or Frameshift in LoF-intolerant gene (pLI > 0.9) | Always | Very Strong |
| Nonsense or Frameshift in LoF-tolerant gene (pLI < 0.9) | Always | Strong (PVS1_Strong) |
| Nonsense or Frameshift without constraint data | Always | Very Strong (benefit of doubt) |
| Splice site | SpliceAI max delta > 0.8 | Very Strong |

When an explicit `lof_gene_list` is provided to the classifier, genes on the list always receive Very Strong PVS1. Genes not on the list receive Strong (downgraded), regardless of pLI. This allows labs to curate a trusted list of LoF-mechanism genes.

If a splice-site variant lacks SpliceAI data, PVS1 is not assigned and "SpliceAI" is recorded as a missing source.

**Required data:** Functional consequence (from GTF annotation). SpliceAI scores for splice-site variants. gnomAD constraint data (pLI) for strength determination.

**Limitation:** PVS1 does not check the following per the 2018 ClinGen PVS1 specification (Tayoun et al. 2018): variants near the 3' end of the transcript, alternatively spliced exons, or single-exon genes. These refinements are planned for a future release.

### PS1 (Strong)

Same amino acid change as an established pathogenic variant, arising from a different nucleotide change.

This criterion applies only to missense variants. It fires when:

1. The variant produces a specific amino acid substitution (e.g., M1775R in BRCA1)
2. ClinVar has a Pathogenic variant at that exact position with that exact substitution
3. The ClinVar variant uses a different codon (different nucleotide change)

The logic: if a different DNA change producing the same protein-level effect is already known to be pathogenic, the query variant is expected to have the same functional impact.

**Required data:**

- Reference FASTA (for codon resolution to determine amino acid change)
- ClinVar protein index TSV (pre-processed list of pathogenic missense variants with protein annotations)

**Does not fire when:** The query variant is the exact same nucleotide change as the ClinVar entry (that would be PP5, not PS1).

### PM2 (Moderate)

Absent from controls (or at extremely low frequency if recessive).

Fires when all population-specific allele frequencies are below 0.0001. Uses gnomAD per-population data when available (AFR, AMR, ASJ, EAS, FIN, NFE, SAS). If any single population exceeds the threshold, PM2 does not fire.

Falls back to global allele frequency when per-population data is absent. When all frequency fields are None, "gnomAD" is recorded as missing.

**Required data:** gnomAD allele frequencies (local TSV/tabix or API).

### PM5 (Moderate)

Novel missense change at an amino acid residue where a different pathogenic missense change has been observed.

This criterion applies only to missense variants. It fires when:

1. The variant introduces a missense change at position X
2. ClinVar has a Pathogenic missense variant at position X
3. The ClinVar variant produces a *different* amino acid substitution than the query

The logic: the position is sensitive to missense variation, so novel missense changes there carry elevated risk.

**Required data:** Same as PS1 (reference FASTA + ClinVar protein index).

**Does not fire when:** PS1 already fired for this variant. PS1 is stronger evidence (same change) and subsumes PM5 (different change at same position).

### PP3 (Supporting or Moderate)

Computational evidence supports a deleterious effect.

Thresholds are calibrated per ClinGen SVI recommendations (Pejaver et al., Am J Hum Genet, 2022). The calibration study validated REVEL thresholds against ClinVar truth sets and established that different score ranges correspond to different evidence strengths.

| Strength | Threshold | Rationale |
| -------- | --------- | --------- |
| Moderate | REVEL > 0.773 | Odds of pathogenicity > 4.33 (moderate evidence) |
| Supporting | REVEL > 0.644 | Odds of pathogenicity > 2.08 (supporting evidence) |
| Supporting | SpliceAI > 0.5 on splice-adjacent variant | Moderate splice disruption predicted |

Only the highest applicable strength fires. When REVEL > 0.773, PP3 fires at moderate strength (not both moderate and supporting).

SpliceAI-based PP3 applies only to variants with SPLICE_SITE or MISSENSE consequence near splice junctions.

**Required data:** REVEL scores and/or SpliceAI scores. When neither is available, both are recorded as missing.

**Previous behavior (v0.13.0):** PP3 fired at supporting level only, with REVEL > 0.7. The ClinGen-calibrated thresholds in v0.14.0 lower the supporting threshold to 0.644 and add a moderate tier at 0.773.

### PP5 (Supporting)

Reputable source (ClinVar) reports the variant as Pathogenic.

Fires when the variant has a ClinVar Pathogenic assertion and no conflicting Benign or Likely Benign assertion exists in ClinVar for the same variant.

**Required data:** ClinVar annotations.

## Benign criteria

### BA1 (Standalone)

Allele frequency is above 5% in any gnomAD population.

This is standalone evidence: BA1 alone classifies a variant as Benign without needing any other criteria. It represents a hard frequency ceiling above which pathogenicity is implausible for Mendelian disease.

Uses population-specific frequencies. If any single population (AFR, AMR, ASJ, EAS, FIN, NFE, SAS) or global AF exceeds 5%, BA1 fires.

**Required data:** gnomAD allele frequencies.

### BS1 (Strong)

Allele frequency is above 1% in any population.

Only evaluated when BA1 did not fire (BA1 already covers frequencies above 5%). BS1 catches variants that are too common for rare disease but below the standalone threshold.

**Required data:** gnomAD allele frequencies.

### BP4 (Supporting or Moderate)

Computational evidence supports no impact on gene or gene product.

ClinGen-calibrated thresholds (Pejaver et al., 2022):

| Strength | Condition | Rationale |
| -------- | --------- | --------- |
| Moderate | Missense + REVEL < 0.183 | Odds of pathogenicity < 0.23 (moderate benign) |
| Supporting | Missense + REVEL < 0.290 | Odds of pathogenicity < 0.48 (supporting benign) |
| Supporting | Non-missense + CADD Phred < 10 | Low overall deleteriousness |

Only the strongest applicable level fires. Does not apply to nonsense or frameshift variants (low computational scores on a truncating variant are not meaningful benign evidence).

**Required data:** REVEL scores (for missense) or CADD scores (for non-missense).

**Previous behavior (v0.13.0):** BP4 fired at supporting level only, with REVEL < 0.15. The ClinGen-calibrated thresholds raise the supporting boundary to 0.290 and add a moderate tier at 0.183.

### BP7 (Supporting)

Synonymous variant with no predicted splice impact.

Fires when:

1. Consequence is SYNONYMOUS
2. SpliceAI max delta < 0.1

The logic: a synonymous change that doesn't affect splicing has no plausible mechanism of pathogenicity.

**Required data:** SpliceAI scores. When unavailable, BP7 is not assigned and "SpliceAI" is recorded as missing.

## Combining rules

Evidence tags combine into a final classification following ACMG/AMP 2015 Table 5.

### Pathogenic

| Rule | Example |
| ---- | ------- |
| 1 Very Strong + 1 Strong | PVS1 + PS1 |
| 2 Strong + 1 Supporting | PS1 + (hypothetical second Strong) + PP3 |
| 1 Very Strong + 2 Supporting | PVS1 + PP3 + PP5 |

### Likely Pathogenic

| Rule | Example |
| ---- | ------- |
| 1 Very Strong + 1 Moderate | PVS1 + PM2 |
| 1 Strong + 1 or more Moderate | PS1 + PM2, or PS1 + PM2 + PM5 |
| 1 Strong + 2 Supporting | PS1 + PP3 + PP5 |

Note: having more moderate evidence than the minimum (e.g., 1S + 3M) still qualifies. The threshold is a floor, not a bounded range.

### Benign

| Rule | Example |
| ---- | ------- |
| 1 BA (Standalone) | BA1 alone |
| 2 Strong benign | BS1 + BS2 |

Note: BS2 (observed in a healthy adult with full penetrance for a recessive condition) exists as a tag in the combining rules but has no evaluator implemented yet. It requires gnomAD homozygote count data that is not currently parsed. The 2-Strong-benign rule activates only if BS2 is manually assigned or added in a future release.

### Likely Benign

| Rule | Example |
| ---- | ------- |
| 1 Strong benign + 1 Supporting benign | BS1 + BP4 |
| 1 Strong benign + 1 Moderate benign | BS1 + BP4_Moderate |
| 2 Moderate benign | BP4_Moderate + BP4_Moderate |
| 1 Moderate benign + 2 Supporting benign | BP4_Moderate + BP4 + BP7 |

### Conflicting evidence

When both pathogenic and benign tags are present on the same variant, the result is VUS. The system does not attempt to weigh conflicting evidence against each other.

### VUS (default)

Any tag combination that doesn't meet the thresholds above results in VUS. This includes:

- Single moderate pathogenic evidence (PM2 alone)
- Single supporting benign evidence (BP7 alone)
- No evidence at all (when all data sources are missing)

## Enabling PS1/PM5

PS1 and PM5 require two things: a reference FASTA (for codon resolution) and a ClinVar protein index file.

### Step 1: Reference FASTA

Provide an indexed reference genome FASTA. This enables codon-level consequence calling, which determines the exact amino acid change for each variant.

```bash
# Download GRCh38 reference (if you don't have one)
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Step 2: ClinVar protein index

The protein index is a TSV of ClinVar pathogenic missense variants annotated with their amino acid changes. It's keyed by (gene, amino acid position) for fast lookup.

**Generate the index:**

```bash
python scripts/prepare_clinvar_protein_index.py \
  --clinvar-vcf clinvar_20240101.vcf.gz \
  --reference-fasta GRCh38.fa \
  --gene-annotation gencode.v44.gtf \
  --output clinvar_protein_index.tsv
```

The script filters ClinVar for Pathogenic/Likely Pathogenic missense variants, resolves their amino acid changes using the reference FASTA, and writes the index.

**TSV format:**

<!-- markdownlint-disable MD010 -->
```text
gene	position	ref_aa	alt_aa	chrom	pos	ref	alt	significance
BRCA1	1775	M	R	chr17	43091429	T	G	Pathogenic
TP53	175	R	H	chr17	7675088	G	A	Pathogenic
```
<!-- markdownlint-enable MD010 -->

### Step 3: Run with both

```bash
vartriage \
  --vcf patient.vcf.gz \
  --output report.json \
  --reference-fasta GRCh38.fa \
  --gene-annotation gencode.v44.gtf \
  --gnomad gnomad.v4.sites.tsv \
  --clinvar clinvar.tsv \
  --revel-scores revel.tsv \
  --clinvar-protein-index clinvar_protein_index.tsv
```

### Python API

```python
from pathlib import Path
from vartriage import ACMGClassifier
from vartriage.annotation.clinvar_protein_index import ClinVarProteinIndex

# Load the protein index
protein_index = ClinVarProteinIndex()
protein_index.load(Path("clinvar_protein_index.tsv"))

# Create classifier with PS1/PM5 enabled
classifier = ACMGClassifier(protein_index=protein_index)

# Classify variants (from your pipeline's scored variant stream)
for classified in classifier.classify(scored_variants):
    print(f"{classified.scored.annotated.variant}: "
          f"{classified.classification.value} "
          f"[{', '.join(t.value for t in classified.evidence_tags)}]")
```

### Without the protein index

When no protein index is provided, PS1 and PM5 are simply omitted. The classifier records "ClinVar_protein_index" in `missing_data_sources` for any missense variant that could have been evaluated. All other criteria function identically to v0.13.0.

## gnomAD API for population frequencies

In API mode or when local gnomAD files are unavailable, the gnomAD GraphQL client provides per-population allele frequencies directly from gnomAD v4.

```bash
# API mode: gnomAD frequencies via GraphQL
vartriage --vcf panel.vcf --output results.json --mode api
```

The client queries `gnomad.broadinstitute.org/api` with a GraphQL variant lookup. It returns frequencies for all 7 ancestry groups (AFR, AMR, ASJ, EAS, FIN, NFE, SAS) plus global AF. Responses are cached in the local SQLite database at `~/.vartriage/api_cache.db`.

See [API Mode Guide](api-mode.md) for full configuration.

## Migration from v0.13.0

### Threshold changes

| Criterion | v0.13.0 | v0.14.0 |
| --------- | ------- | ------- |
| PP3 (Supporting) | REVEL > 0.7 | REVEL > 0.644 |
| PP3 (Moderate) | not available | REVEL > 0.773 |
| BP4 (Supporting) | REVEL < 0.15 | REVEL < 0.290 |
| BP4 (Moderate) | not available | REVEL < 0.183 |

### Impact on existing analyses

Variants with REVEL between 0.644 and 0.7 now receive PP3 (previously they did not). This means some VUS variants may shift to Likely Pathogenic if they already had other moderate or strong pathogenic evidence.

Variants with REVEL between 0.15 and 0.290 now receive BP4 (previously they did not). Some VUS variants may shift to Likely Benign if they also carry strong benign evidence (BS1).

To see the impact: re-run your existing analyses with v0.14.0 and diff the classification column.

### New criteria (additive)

PS1 and PM5 only fire when a protein index is provided. Without it, classification behavior is identical to v0.13.0. You can upgrade the package and only enable PS1/PM5 when you're ready to prepare the protein index.

## Evidence strength map

Complete mapping of tags to strength tiers in vartriage:

| Tag | Strength | Direction | Status |
| --- | -------- | --------- | ------ |
| PVS1 | Very Strong | Pathogenic | Evaluated |
| PVS1_Strong | Strong | Pathogenic | Evaluated (v0.17.0) |
| PS1 | Strong | Pathogenic | Evaluated |
| PM1 | Moderate | Pathogenic | Evaluated (v0.17.0) |
| PM2 | Moderate | Pathogenic | Evaluated |
| PM4 | Moderate | Pathogenic | Evaluated (v0.17.0) |
| PM5 | Moderate | Pathogenic | Evaluated |
| PP3 | Supporting | Pathogenic | Evaluated |
| PP3_MODERATE | Moderate | Pathogenic | Evaluated |
| PP5 | Supporting | Pathogenic | Evaluated |
| BA1 | Standalone | Benign | Evaluated |
| BS1 | Strong | Benign | Evaluated |
| BS2 | Strong | Benign | Placeholder (not evaluated) |
| BP4 | Supporting | Benign | Evaluated |
| BP4_MODERATE | Moderate | Benign | Evaluated |
| BP7 | Supporting | Benign | Evaluated |

### Criteria not implemented

The following ACMG/AMP 2015 criteria are not currently evaluated:

| Criterion | Reason |
| --------- | ------ |
| PM3 | Requires phase data (detected in trans with pathogenic variant) |
| PP1 | Requires pedigree cosegregation data |
| PP2 | Requires per-gene benign missense rate data |
| BS2 | Requires gnomAD homozygote count data |
| BS3 | Requires functional assay data |

## Missing data handling

The classifier tracks which data sources were unavailable for each variant. This is useful for interpreting VUS results: a VUS with three missing sources may warrant re-analysis once those sources become available.

Common missing data patterns:

| Missing source | Criteria affected | How to resolve |
| -------------- | ----------------- | -------------- |
| gnomAD | PM2, BA1, BS1 | Provide gnomAD frequencies (local or API) |
| gnomAD_constraint | PVS1 strength | Provide gene knowledge data (--knowledge-dir or bundled) |
| functional_domain | PM1 | Provide gene knowledge data with constraint metrics |
| REVEL | PP3, BP4 | Provide REVEL scores TSV |
| SpliceAI | PVS1 (splice), PP3 (splice), BP7 | Provide SpliceAI scores TSV |
| ClinVar | PP5 | Provide ClinVar annotation file |
| ClinVar_protein_index | PS1, PM5 | Generate and provide the protein index |
| codon_resolution | PS1, PM5 | Provide reference FASTA |

## Conflicting evidence

When both pathogenic and benign evidence tags are present for the same variant, the classifier returns VUS and sets `has_conflicting_evidence = True` on the output. This distinguishes "VUS because insufficient evidence" from "VUS because contradictory evidence." Clinical reports display a note when this flag is set.

## Limitations and assumptions

**Multiallelic input:** vartriage expects biallelic VCF input. Multiallelic sites must be decomposed upstream (e.g., via `bcftools norm -m-`). Passing multiallelic records may produce incorrect annotation lookups.

**Genome build:** Coordinates are build-agnostic internally. Reference files (GTF, gnomAD, ClinVar, CADD, SpliceAI) must match the VCF's genome build. The library ships GRCh38 presets and bundles; GRCh37 is supported with user-supplied references.

**Mitochondrial haplogroup context:** Haplogroup assignment is not performed. Variants defining rare haplogroups may appear rare in HelixMTdb and receive false-positive pathogenicity scores. Manual haplogroup confirmation is recommended for ambiguous mtDNA results.

**PP5 single-assertion model:** ClinVar entries frequently have multiple submissions with conflicting interpretations. The library stores a single assertion per variant. A multi-submitter ClinVar entry classified as "Pathogenic" by one lab and "Likely Benign" by another will show whichever assertion the reference file provides. Future releases may incorporate review_status filtering.

**Disease-specific thresholds:** PM2, BA1, and BS1 use fixed population frequency thresholds regardless of disease inheritance mode. Disease-specific carrier frequency calibration is not implemented.

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. Genet Med. 2015;17(5):405-424.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. Am J Hum Genet. 2022;109(12):2163-2177.
- Riggs ER, et al. Technical standards for the interpretation and reporting of constitutional copy-number variants. Genet Med. 2020;22(2):245-257.
- Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. Hum Mutat. 2018;39(11):1517-1524.
