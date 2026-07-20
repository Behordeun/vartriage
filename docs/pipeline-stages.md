# Pipeline Stages

## Data flow

```text
VCFParser → [SampleExtractor|InheritanceFilter] → [RegionFilter] → QualityFilter → AnnotationEngine → [GeneFilter] → [SecondaryFindingsFilter] → PrioritizationEngine → ACMGClassifier → ReportGenerator
```

Stages in brackets are optional. They activate based on configuration:

- `SampleExtractor` activates with `--sample` (mutually exclusive with InheritanceFilter)
- `InheritanceFilter` activates with `--proband/--mother/--father`
- `RegionFilter` activates with `--regions`
- `GeneFilter` activates with `--gene-list`
- `SecondaryFindingsFilter` activates with `--secondary-findings`

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

## Sample Extraction (optional, v0.3.0+)

**Class:** `SampleExtractor`

Extracts a single sample from multi-sample VCFs, keeping only variants where the named sample carries an alternate allele.

**Input:** Iterator of `Variant` (with `_pysam_samples` in info dict).

**Output:** Iterator of `Variant` (subset, with `sample_gt` and `sample_name` in info dict).

**Behavior:**

- Looks up the named sample in the VCF header
- Filters to variants where the sample's genotype contains at least one alt allele
- Optionally applies a genotype quality (GQ) threshold
- Variants below GQ threshold are dropped
- Populates `zygosity` and `quality_metrics` on downstream `AnnotatedVariant` (v0.10.0+)

**Errors:**

- `ValueError` if the sample name does not exist in the VCF header

**Configuration:** `SampleConfig(sample_name="PROBAND_01", min_gq=20)`

**CLI:** `--sample PROBAND_01 --min-gq 20`

Mutually exclusive with trio arguments (`--proband/--mother/--father`).

## Region Filtering (optional, v0.3.0+)

**Class:** `RegionFilter`

Restricts analysis to variants overlapping genomic intervals from a BED file.

**Input:** Iterator of `Variant`.

**Output:** Iterator of `Variant` (subset overlapping at least one BED interval).

**Behavior:**

- Parses BED file at construction (0-based half-open intervals)
- Builds an interval index per chromosome
- Each variant is checked for overlap with any target interval
- Variants outside all intervals are dropped

**Errors:**

- `FileNotFoundError` if BED file does not exist
- `ValueError` if BED file is malformed (fewer than 3 columns, non-numeric coordinates)

**Configuration:** `RegionFilterConfig(bed_path=Path("target_regions.bed"))`

**CLI:** `--regions target_regions.bed`

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
| --------- | ----------- |
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
| ------------ | -------- | --------------- |
| Functional consequence | GTF/GFF gene models | Coordinate overlap |
| Population frequency | gnomAD | Exact match (chrom, pos, ref, alt) |
| ClinVar assertion | ClinVar | Exact match (chrom, pos, ref, alt) |

**Backend auto-detection:**

The engine picks the fastest available backend at construction:

- Consequence: PyRanges (if installed) or pure-Python sorted interval tree
- Frequency: Polars (if installed) or dict-based lookup; tabix VCF (zero-memory) for `.vcf.bgz`/`.vcf.gz` gnomAD files
- ClinVar: Polars (if installed) or dict-based lookup

Both produce the same results. The accelerated backends run faster on large reference files.

**Codon-level consequence calling (v0.8.0+):**

When `reference_fasta_path` is configured, the engine uses actual amino acid comparison for SNVs in CDS regions instead of the positional heuristic. `TranscriptCDSIndex` maps genomic coordinates to CDS positions via the GTF, then `CodonResolver` extracts the reference codon from FASTA, substitutes the variant base, translates both, and compares. Correctly distinguishes synonymous, missense, and nonsense. Without FASTA, the fallback heuristic ("SNV in CDS = Missense") is used.

**Variant normalization (v0.8.0+):**

When a reference FASTA is available, indels are left-aligned and trimmed before database lookups. This reduces silent lookup failures from representation differences between callers and databases. The algorithm follows Tan et al. 2015 with a 1000-iteration safety cap.

**Population-specific frequencies (v0.9.0+):**

When gnomAD data includes per-population columns (AFR, AMR, ASJ, EAS, FIN, NFE, SAS), these are stored in `PopulationFrequencies` on `AnnotatedVariant`. PM2 checks ALL populations below threshold. BA1/BS1 check if ANY population exceeds threshold.

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

## Secondary Findings Screening (optional, v0.10.0+)

**Class:** `SecondaryFindingsFilter`

Flags variants in the 71 ACMG Secondary Findings (SF v3.2) medically actionable genes, regardless of whether they pass the primary gene panel filter.

**Input:** Iterator of `AnnotatedVariant`.

**Output:** Iterator of `AnnotatedVariant` (unchanged, but flagged variants are tracked internally).

**Behavior:**

- Loads the built-in gene list from `vartriage/data/acmg_sf_v3.2.txt` (shipped as package data)
- Checks each variant's `gene_name` against the 71-gene set
- Does not remove any variants from the stream
- `is_secondary_finding(gene_name)` method for point queries
- `split_stream()` method separates the variant stream into primary and secondary finding subsets

Clinical reports include a dedicated secondary findings section when this filter is active.

**CLI:** `--secondary-findings`

## Prioritization

**Class:** `PrioritizationEngine`

Filters by allele frequency and computes composite pathogenicity scores.

**Input:** Iterator of `AnnotatedVariant`.

**Output:** Iterator of `ScoredVariant` (sorted descending by composite rank within each batch).

**Two phases:**

1. **Frequency gate:** Drops variants with `allele_frequency > max_allele_frequency` (default 0.01). Variants marked `frequency_unknown` always pass.

2. **Composite scoring:** Normalizes and combines CADD + REVEL + SpliceAI scores:
   - CADD normalized: `min(cadd_phred / 99.0, 1.0)`
   - Three scores present: `(REVEL * 0.5) + (CADD_normalized * 0.3) + (SpliceAI * 0.2)`
   - Two scores present: weights redistribute proportionally among available scores
   - Single score present: used as-is
   - Falls back to legacy 0.6/0.4 formula when SpliceAI is not configured
   - Variants without any scores get `composite_rank=None` and sort last

**Memory safety:** On `MemoryError`, automatically retries with smaller chunk sizes (capped at 500,000 per chunk).

**Prioritization score (v0.8.0+):**

A separate `prioritization_score` field uses literature-validated scoring:
- Missense variants: REVEL score directly (threshold 0.7 validated against ClinGen)
- Splice-adjacent variants: SpliceAI delta score
- Non-missense: CADD Phred / 60

This is the recommended ranking method. `composite_rank` (legacy weighted average) is kept for backward compatibility but deprecated. Both are present in output; `composite_rank` will be removed in v1.0.0.

**Configuration:** `PrioritizationConfig(max_allele_frequency=0.01, cadd_scores_path=None, revel_scores_path=None, spliceai_scores_path=None, batch_size=10_000)`

## ACMG Classification

**Class:** `ACMGClassifier`

Assigns ACMG/AMP 2015 evidence tags and determines final classification.

**Input:** Iterator of `ScoredVariant`.

**Output:** Iterator of `ClassifiedVariant`.

**Evidence criteria evaluated:**

Pathogenic criteria:

| Tag | Strength | Condition |
| ----- | ---------- | ----------- |
| PVS1 | Very Strong | Consequence is Nonsense, Frameshift, or Splice_Site + SpliceAI > 0.8 |
| PM2 | Moderate | All population AFs < 0.0001 (population-specific when available) |
| PP3 | Supporting | REVEL score > 0.7, or SpliceAI > 0.5 on splice-adjacent variant |
| PP5 | Supporting | ClinVar Pathogenic, no conflicting Benign/Likely_Benign |

Benign criteria (v0.9.0+):

| Tag | Strength | Condition |
| ----- | ---------- | ----------- |
| BA1 | Standalone | Any population AF > 5% |
| BS1 | Strong | Any population AF > 1% (only when BA1 not already assigned) |
| BP4 | Supporting | Missense with REVEL < 0.15, or non-missense with CADD Phred < 10 |
| BP7 | Supporting | Synonymous with SpliceAI < 0.1 |

**Combining rules:**

Tags combine into a final classification across all five ACMG tiers:

- **Pathogenic:** 1 Very Strong + 1 Moderate, or 1 Very Strong + 2 Supporting
- **Likely_Pathogenic:** 1 Very Strong + 1 Supporting, or 1 Moderate + 2 Supporting
- **VUS:** Insufficient evidence for either direction, or conflicting pathogenic + benign evidence
- **Likely_Benign:** 1 Strong benign + 1 Supporting benign
- **Benign:** BA1 alone (standalone), or 2 Strong benign

When a required data source is unavailable for a criterion, that tag is omitted and the source name is recorded in `missing_data_sources`.

## Report Generation

**Class:** `ReportGenerator`

Serializes classified variants to JSON, CSV, PDF, or clinical report formats.

**Input:** Sequence of `ClassifiedVariant` + output path.

**Output:** File at the specified path.

**Formats:**

- **JSON:** Array of variant objects with all fields
- **CSV:** One row per variant, header row included
- **PDF:** Formatted clinical report (requires `reportlab` or uses fallback renderer)
- **VCF:** Bgzipped VCF with tabix index, injecting VARTRIAGE_* INFO fields for classified variants
- **Clinical HTML/PDF/DOCX:** Structured clinical reports with evidence narratives (see below)

**Atomicity:** Writes to a temp file first, then performs an atomic rename. If the write fails, no partial output exists at the target path.

**Configuration:** `ReportConfig(output_format="json")`

## Clinical Report Generation

**Class:** `ClinicalReportGenerator`

Produces structured, sign-off-ready clinical variant reports from ClassifiedVariant data. Activated when `--output-format` is `clinical-html`, `clinical-pdf`, or `clinical-docx`.

**Input:** Iterator of `ClassifiedVariant` + `ClinicalReportConfig` + output path.

**Output:** Clinical report file + `.audit.json` sidecar.

**Report sections (in order):**

1. **Header:** patient ID, panel name, analysis date (ISO 8601), pipeline version.
2. **Executive Summary:** total variants analyzed, count per classification tier (Pathogenic, Likely_Pathogenic, VUS).
3. **Findings Table:** variants ranked by tier (Pathogenic first), then by composite rank descending within tier.
4. **Evidence Cards:** one per variant. Each card contains gene name, consequence, population frequency, predictor scores, ClinVar data, inheritance pattern, and ACMG criteria with plain-language explanations.
5. **Limitations:** lists any data sources that were unavailable during classification.
6. **Methodology:** pipeline version, reference file paths, classification parameters, analysis timestamp.
7. **Sign-off:** placeholder fields for reviewer name, review date, and digital signature.

**Evidence narratives:**

The `EvidenceNarrativeBuilder` transforms raw ClassifiedVariant fields into human-readable text using hardcoded string templates. No LLM or generative AI is invoked. Each narrative includes:

- Gene name and protein consequence
- Allele frequency with denominator context (e.g., "0.000008 (1 in 125,000)")
- Predictor scores with scale context (e.g., "REVEL: 0.95 (scale 0-1, threshold 0.7)")
- ACMG evidence tags with one-sentence explanations of why each tag fired
- Inheritance pattern when available
- Notes on missing data sources

**Audit trail:**

A `.audit.json` sidecar is written alongside the report. Contains:

- Run manifest: all config parameters, reference file paths with SHA-256 checksums, pipeline version, execution timestamp, Python version, platform
- Decision log: one entry per variant with evidence tags assigned, tags skipped (with reasons), intermediate scores, and final classification

**Output formats:**

- `clinical-html`: Self-contained HTML with inlined CSS. No external stylesheets, no JavaScript. Viewable offline in any browser.
- `clinical-pdf`: Renders the HTML to PDF via WeasyPrint. Raises `ImportError` if weasyprint is not installed. Text remains selectable.
- `clinical-docx`: Word document via python-docx. Uses Heading 1, Heading 2, Normal, and Table Grid styles. Raises `ImportError` if python-docx is not installed.

**Empty variant handling:**

When zero variants pass triage, the report includes all sections with zero counts and a statement: "No clinically significant variants meeting triage criteria were identified within the requested panel target areas."

**Configuration:** `ClinicalReportConfig(patient_id="...", panel_name="...", output_format="clinical-html")`
