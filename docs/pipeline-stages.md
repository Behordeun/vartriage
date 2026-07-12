# Pipeline Stages

## Data flow

```text
┌──────────┐    ┌───────────────┐    ┌──────────────────┐    ┌──────────────────────┐    ┌────────────────┐    ┌─────────────────┐
│ VCFParser │───▶│ QualityFilter │───▶│ AnnotationEngine │───▶│ PrioritizationEngine │───▶│ ACMGClassifier │───▶│ ReportGenerator │
└──────────┘    └───────────────┘    └──────────────────┘    └──────────────────────┘    └────────────────┘    └─────────────────┘
    Variant          Variant           AnnotatedVariant           ScoredVariant          ClassifiedVariant         JSON/CSV/PDF
```

Each stage consumes an iterator and yields an iterator. Only one batch lives in memory at a time.

## VCF Parsing

**Class:** `VCFParser`

Streams `Variant` records from `.vcf` or `.vcf.gz` files using pysam (htslib wrapper). Compressed files require a `.tbi` tabix index.

**Input:** File path to a VCF file.

**Output:** Iterator of `Variant` dataclass instances.

**Behavior:**

- Validates the VCF header at construction (fail-fast on missing `##fileformat` or `#CHROM`)
- Yields one record at a time, no buffering
- Multiallelic sites: takes the first ALT allele (split upstream)
- Missing QUAL becomes `None`
- FILTER field preserved as-is ("PASS", ".", or semicolon-joined filter names)

**Errors:**

- `FileNotFoundError` if VCF or `.tbi` index missing
- `ParseError` on malformed headers or data lines

## Quality Filtering

**Class:** `QualityFilter`

Drops variants that fail machine-level quality controls.

**Input:** Iterator of `Variant`.

**Output:** Iterator of `Variant` (subset).

**Rules:**

1. FILTER field must be `"PASS"` or `"."` (missing/not-applied)
2. QUAL score must be present and at least `min_qual` (default 20.0)

Variants with missing QUAL emit a `MissingDataWarning` before being dropped.

**Configuration:** `QualityFilterConfig(min_qual=20.0)`

## Inheritance Pattern Classification (optional)

**Class:** `InheritanceFilter`

When a proband-mother-father trio is configured, replaces SampleExtractor and classifies each variant into Mendelian inheritance patterns based on family genotypes.

**Input:** Iterator of `Variant` (with `_pysam_samples` in info dict).

**Output:** Iterator of `Variant` (with `inheritance_pattern`, `sample_gt`, `sample_name` in info dict).

**Patterns classified:**

| Pattern | Condition |
|---------|-----------|
| de_novo | Proband has alt, both parents hom-ref |
| dominant | Proband het, exactly one parent het |
| recessive | Proband hom-alt, both parents het |
| compound_het | Two+ het variants in same gene from different parents (trans) |
| x_linked | ChrX variant, proband has alt, mother is het carrier |

A single variant can carry multiple labels when it satisfies more than one rule. Variants where the proband is hom-ref or has missing genotype are skipped.

Compound het requires gene annotation, so the pipeline positions InheritanceFilter after AnnotationEngine when compound_het is in the patterns list.

**Configuration:** `InheritanceConfig(proband="CHILD", mother="MOM", father="DAD")`

**CLI:** `--proband CHILD --mother MOM --father DAD [--inheritance-pattern de_novo ...]`

Mutually exclusive with `--sample`.

## Annotation

**Class:** `AnnotationEngine`

Enriches each variant with three annotations, processed in configurable batches (default 10,000).

**Input:** Iterator of `Variant`.

**Output:** Iterator of `AnnotatedVariant`.

**Annotations added:**

| Annotation | Source | Lookup method |
|------------|--------|---------------|
| Functional consequence | GTF/GFF gene models | Coordinate overlap |
| Population frequency | gnomAD | Exact match (chrom, pos, ref, alt) |
| ClinVar assertion | ClinVar | Exact match (chrom, pos, ref, alt) |

**Backend auto-detection:**

The engine picks the fastest available backend at construction:

- Consequence: PyRanges (if installed) or pure-Python sorted interval tree
- Frequency: Polars (if installed) or dict-based lookup
- ClinVar: Polars (if installed) or dict-based lookup

Both produce the same results. The accelerated backends run faster on large reference files.

**Missing data handling:**

- Variant not in gnomAD: `allele_frequency=None`, `frequency_unknown=True`, emits `MissingDataWarning`
- Variant not in ClinVar: `clinvar_assertion=None`, `clinvar_unknown=True`, emits `MissingDataWarning`

**Configuration:** `AnnotationConfig(gene_annotation_path, gnomad_path, clinvar_path=None, batch_size=10_000)`

## Prioritization

**Class:** `PrioritizationEngine`

Filters by allele frequency and computes composite pathogenicity scores.

**Input:** Iterator of `AnnotatedVariant`.

**Output:** Iterator of `ScoredVariant` (sorted descending by composite rank within each batch).

**Two phases:**

1. **Frequency gate:** Drops variants with `allele_frequency > max_allele_frequency` (default 0.01). Variants marked `frequency_unknown` always pass.

2. **Composite scoring:** Normalizes and combines CADD + REVEL scores:
   - CADD normalized: `min(cadd_phred / 99.0, 1.0)`
   - Composite: `(REVEL * 0.6) + (CADD_normalized * 0.4)`
   - Falls back to the single available score when only one source exists
   - Variants without any scores get `composite_rank=None` and sort last

**Memory safety:** On `MemoryError`, automatically retries with smaller chunk sizes (capped at 500,000 per chunk).

**Configuration:** `PrioritizationConfig(max_allele_frequency=0.01, cadd_scores_path=None, revel_scores_path=None, batch_size=10_000)`

## ACMG Classification

**Class:** `ACMGClassifier`

Assigns ACMG/AMP 2015 evidence tags and determines final classification.

**Input:** Iterator of `ScoredVariant`.

**Output:** Iterator of `ClassifiedVariant`.

**Evidence criteria evaluated:**

| Tag | Strength | Condition |
|-----|----------|-----------|
| PVS1 | Very Strong | Consequence is Nonsense or Frameshift |
| PM2 | Moderate | gnomAD AF < 0.0001 |
| PP3 | Supporting | REVEL score > 0.7 |
| PP5 | Supporting | ClinVar Pathogenic, no conflicting Benign/Likely_Benign |

**Combining rules:**

Tags combine into a final classification:

- **Pathogenic:** 1 Very Strong + 1 Moderate, or 1 Very Strong + 2 Supporting
- **Likely_Pathogenic:** 1 Very Strong + 1 Supporting, or 1 Moderate + 2 Supporting
- **VUS:** Everything else

When a required data source is unavailable for a criterion, that tag is omitted and the source name is recorded in `missing_data_sources`.

## Report Generation

**Class:** `ReportGenerator`

Serializes classified variants to JSON, CSV, or PDF.

**Input:** Sequence of `ClassifiedVariant` + output path.

**Output:** File at the specified path.

**Formats:**

- **JSON:** Array of variant objects with all fields
- **CSV:** One row per variant, header row included
- **PDF:** Formatted clinical report (requires `reportlab` or uses fallback renderer)

**Atomicity:** Writes to a temp file first, then performs an atomic rename. If the write fails, no partial output exists at the target path.

**Configuration:** `ReportConfig(output_format="json")`
