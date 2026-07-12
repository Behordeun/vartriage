# Pipeline Stages

## Data flow

```text
VCFParser → QualityFilter → AnnotationEngine → [GeneFilter] → PrioritizationEngine → ACMGClassifier → ReportGenerator
```

Stages in brackets are optional. GeneFilter activates when `--gene-list` is provided.

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

## Gene Filtering (optional)

**Class:** `GeneFilter`

Restricts the annotated variant stream to only those variants whose gene symbol appears in a user-supplied text file.

**Input:** Iterator of `AnnotatedVariant`.

**Output:** Iterator of `AnnotatedVariant` (subset).

**Behavior:**

- Loads a plain text file at construction: one gene symbol per line, comment lines (`#`) and blank lines skipped
- Normalizes all symbols to uppercase for case-insensitive matching
- Yields only variants whose `gene_name` (from annotation) matches a gene in the set
- Intergenic variants (`gene_name=None`) are silently excluded
- After the stream is consumed, logs a WARNING listing any panel genes with zero matching variants

**Errors:**

- `FileNotFoundError` if the gene list file does not exist
- `ValueError` if the file contains zero valid gene symbols

**Configuration:** `GeneFilterConfig(gene_list_path=Path("my_panel.txt"))`

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
