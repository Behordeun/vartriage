# Architecture

Package structure, the Protocol-based backend system, and extension points.

## Package structure

```text
vartriage/
├── __init__.py              # Public API exports
├── pipeline.py              # Top-level orchestrator
├── cli.py                   # Command-line interface (single-sample + cohort subcommand)
├── protocols.py             # Protocol interfaces for swappable backends
├── exceptions.py            # VarTriageWarning base class
├── py.typed                 # PEP 561 marker
├── io/
│   ├── vcf_parser.py       # pysam-based VCF streaming with sample extraction
│   └── exceptions.py       # Error hierarchy (ParseError, ConfigurationError, etc.)
├── filter/
│   ├── quality_filter.py   # FILTER + QUAL gate
│   ├── gene_filter.py      # Gene-list-based variant filtering (v0.4.0)
│   ├── region_filter.py    # BED-based genomic interval filtering (v0.3.0)
│   ├── sample_extractor.py # Single-sample extraction from multi-sample VCFs (v0.3.0)
│   ├── inheritance_filter.py  # Trio-based inheritance patterns (v0.4.0)
│   └── secondary_findings.py # ACMG SF v3.2 screening (v0.10.0)
├── annotation/
│   ├── engine.py            # Orchestrator with backend auto-detection
│   ├── consequence.py       # Pure-Python interval tree consequence lookup
│   ├── consequence_pyranges.py  # PyRanges backend (vectorized batch join)
│   ├── codon_resolver.py   # FASTA-backed codon extraction + translation (v0.8.0)
│   ├── transcript_index.py # Per-transcript CDS exon map from GTF (v0.8.0)
│   ├── frequency.py         # Dict-based gnomAD lookup
│   ├── frequency_polars.py  # Polars backend (optional)
│   ├── frequency_tabix.py   # Tabix VCF backend (pysam, zero-memory)
│   ├── clinvar.py           # Dict-based ClinVar lookup
│   └── clinvar_polars.py    # Polars backend (optional)
├── prioritization/
│   ├── engine.py            # AF gating + scoring orchestrator
│   ├── frequency_filter.py  # Allele frequency threshold filter
│   ├── score_loader.py     # CADD/REVEL/SpliceAI TSV file loading
│   └── scoring.py           # Score normalization + composite + prioritization_score
├── classification/
│   ├── acmg.py              # Evidence tag assignment (pathogenic + benign criteria)
│   └── combining.py         # ACMG combining rules (all 5 tiers)
├── cohort/                  # Multi-sample cohort analysis (v0.11.0)
│   ├── __init__.py          # Public exports
│   ├── aggregator.py        # Cross-sample variant merging by coordinate
│   ├── statistics.py        # Gene burden + recurrence stats
│   ├── report.py            # JSON/CSV cohort report writer
│   └── pipeline.py          # Multi-sample orchestrator
├── reporting/
│   ├── generator.py         # Format routing + atomic write
│   ├── json_writer.py       # JSON serialization (streaming)
│   ├── csv_writer.py        # CSV serialization (streaming)
│   ├── vcf_writer.py        # Annotated VCF output with tabix index (v0.4.0)
│   ├── pdf_writer.py        # ReportLab PDF (optional)
│   ├── pdf_fallback.py      # Fallback when reportlab not installed
│   └── clinical/            # Clinical report generation (v0.5.0)
│       ├── generator.py     # Clinical report orchestrator
│       ├── models.py        # Clinical report data structures
│       ├── narrative.py     # Per-variant evidence narrative builder
│       ├── template_engine.py  # Template rendering engine
│       ├── templates.py     # HTML/DOCX template definitions
│       └── audit.py         # JSON audit trail sidecar writer
├── bundle/                  # Score bundle downloader (v0.6.0)
│   ├── cli.py              # Bundle subcommand (download, list, verify, status)
│   ├── config.py           # TOML configuration loader
│   ├── downloader.py       # HTTP download with resume + retry
│   ├── transformer.py      # Post-download format conversion
│   ├── registry.py         # Bundle registry (available bundles + versions)
│   ├── registry.json       # Static registry data
│   ├── storage.py          # Path resolution + bundle layout management
│   ├── manifest.py         # Per-bundle manifest tracking
│   ├── _checksums.py       # SHA-256 verification
│   ├── _disk.py            # Disk space pre-flight checks
│   └── _progress.py        # Download progress display
├── api/                     # API annotation backend (v0.7.0)
│   ├── __init__.py          # Lazy exports
│   ├── _base.py             # BaseAPIClient (retry, rate limit, circuit breaker)
│   ├── _rate_limiter.py     # Token bucket with daily caps
│   ├── _circuit_breaker.py  # CLOSED/OPEN/HALF_OPEN state machine
│   ├── _cache.py            # SQLite response cache with TTL
│   ├── _notation.py         # VCF-to-VEP coordinate converter
│   ├── _consequence_map.py  # SO term to FunctionalConsequence mapping
│   ├── config.py            # APIConfig dataclass
│   ├── vep_client.py        # Ensembl VEP batch POST client
│   ├── clinvar_client.py    # NCBI ClinVar E-utilities client
│   ├── cadd_client.py       # CADD REST score lookups
│   ├── spliceai_client.py   # SpliceAI Lookup with smart filtering
│   ├── annotation_engine.py # Composes VEP + ClinVar into annotate() interface
│   └── score_provider.py    # CADD hierarchy + SpliceAI
├── models/
│   ├── config.py            # All config dataclasses
│   ├── variant.py           # Variant, AnnotatedVariant, ScoredVariant, ClassifiedVariant, enums
│   ├── cohort.py            # CohortConfig, CohortVariant, GeneBurden, CohortSummary
│   └── warnings.py          # MissingDataWarning
├── data/
│   └── acmg_sf_v3.2.txt    # 71-gene ACMG Secondary Findings list (v0.10.0)
└── _internal/
    ├── batch.py             # Batch iteration utilities
    ├── cache.py             # Pickle caching with mtime invalidation
    ├── genetic_code.py      # Standard genetic code + translate_codon() (v0.8.0)
    ├── normalizer.py        # Left-align + trim indel normalization (v0.8.0)
    ├── interval_tree.py     # Sorted-array interval tree
    ├── vectorized.py        # NumPy vectorized operations
    └── warning_accumulator.py  # Warning collection + threshold
```

## Filtering stages (v0.3.0+)

The pipeline supports several optional filtering stages activated by configuration:

```text
VCFParser → [SampleExtractor] → [InheritanceFilter] → [RegionFilter] → QualityFilter → AnnotationEngine → [GeneFilter] → [SecondaryFindingsFilter] → PrioritizationEngine → ACMGClassifier → ReportGenerator
```

Stages in brackets activate conditionally:

- **SampleExtractor** (`--sample`): pulls a single sample from multi-sample VCFs, filtering to variants where the named sample carries an alt allele. Applies optional GQ threshold.
- **InheritanceFilter** (`--proband/--mother/--father`): classifies variants into Mendelian inheritance patterns (de_novo, dominant, recessive, compound_het, x_linked) based on trio genotypes. Mutually exclusive with SampleExtractor.
- **RegionFilter** (`--regions`): restricts to variants overlapping BED intervals.
- **GeneFilter** (`--gene-list`): post-annotation filter keeping only variants in specified genes.
- **SecondaryFindingsFilter** (`--secondary-findings`): flags variants in ACMG SF v3.2 genes regardless of primary gene panel.

## Codon-level consequence calling (v0.8.0+)

When `--reference-fasta` is provided, the annotation engine uses a two-step process for SNVs in CDS regions:

1. `TranscriptCDSIndex` builds a per-transcript exon map from GTF CDS features, mapping genomic coordinates to CDS positions with strand awareness.
2. `CodonResolver` extracts the reference codon from the FASTA (handling split codons at exon junctions), substitutes the variant base, translates both codons via `genetic_code.translate_codon()`, and compares amino acids.

This correctly distinguishes synonymous, missense, and nonsense changes. Without the FASTA, the positional heuristic ("SNV in CDS = Missense") is used as a fallback.

## Variant normalization (v0.8.0+)

`_internal/normalizer.py` implements left-alignment and trimming for indels before database lookups. The algorithm (Tan et al. 2015):

1. Right-trim matching suffix bases
2. Left-trim matching prefix bases
3. Left-align by shifting the variant leftward while the last base of ref matches the last base of alt

A 1000-iteration safety cap prevents infinite loops on pathological inputs. Normalization reduces silent lookup failures caused by representation differences between VCF callers and reference databases.

## Bundle downloader (v0.6.0+)

```text
vartriage/bundle/
├── cli.py              # Subcommands: download, list, verify, status, update-registry
├── config.py           # TOML config (~/.vartriage/config.toml)
├── downloader.py       # HTTP GET with Range resume, exponential backoff, progress bar
├── transformer.py      # Post-download conversion (VCF-to-TSV, ClinVar normalize, etc.)
├── registry.py         # Available bundles: clinvar, gnomad-exomes-chr22, revel, gencode, spliceai
├── storage.py          # Layout: ~/.vartriage/bundles/{build}/{bundle_name}/
├── manifest.py         # Per-bundle tracking (version, checksum, timestamp)
├── _checksums.py       # SHA-256 file verification
├── _disk.py            # Free space pre-flight
└── _progress.py        # stderr progress bar
```

Bundles are reference files packaged for automatic download and transformation. `--use-bundles` in the CLI auto-resolves paths from installed bundles. Explicit path flags always take precedence.

Storage layout: `~/.vartriage/bundles/grch38/clinvar/clinvar.tsv` (configurable via `VARTRIAGE_BUNDLE_STORAGE` env var).

## Clinical report generation (v0.5.0+)

```text
vartriage/reporting/clinical/
├── generator.py         # Orchestrator: collects classified variants, routes to template
├── models.py            # ClinicalReportData, FindingEntry, EvidenceCard structures
├── narrative.py         # Plain-language ACMG evidence descriptions per variant
├── template_engine.py   # Renders HTML (self-contained, inlined CSS) or DOCX
├── templates.py         # Section templates: header, summary, findings, limitations, methodology
└── audit.py             # JSON audit trail sidecar (.audit.json) with SHA-256 checksums
```

Three output formats: `clinical-html` (self-contained, no JS), `clinical-pdf` (via WeasyPrint), `clinical-docx` (via python-docx). Each produces an audit trail sidecar recording the run manifest, reference file checksums, and per-variant decision log.

Reports include: computational-only disclaimer (citing ACMG/AMP 2015), executive summary, findings table sorted by classification tier, per-variant evidence cards with criterion explanations, limitations section, methodology description, and sign-off block.

## API annotation backend (v0.7.0+)

The API package provides an alternative annotation backend that queries remote services instead of local files. All components use deferred httpx imports so that `import vartriage` works without httpx installed.

The `BaseAPIClient` composes rate limiting, circuit breaking, and retry into a single request path shared by all service-specific clients:

- **Rate limiter**: token bucket with configurable per-service limits and daily caps
- **Circuit breaker**: CLOSED/OPEN/HALF_OPEN state machine, opens after 5 failures in 60s, half-open recovery after 30s
- **Retry**: exponential backoff (1s, 2s, 4s) with Retry-After header parsing
- **Cache**: SQLite-backed with configurable TTL (7 days default, 30 days for ClinVar, pinned mode for clinical reproducibility)

See [API Mode Guide](api-mode.md) for usage and configuration.

## Cohort analysis module (v0.11.0+)

The cohort module builds on top of the standard pipeline. `CohortPipeline` runs each sample VCF through the full `Pipeline` (parse, filter, annotate, score, classify), collects `ClassifiedVariant` results, then feeds them to `CohortAggregator`. The aggregator groups variants by `(chrom, pos, ref, alt)`, resolves cross-sample severity, and produces `CohortVariant` records. `CohortStatistics` computes derived metrics (gene burden, recurrence distribution) and `CohortReportGenerator` serializes everything to disk.

Parallel processing uses `concurrent.futures.ThreadPoolExecutor` since per-sample work is I/O-bound (VCF parsing via pysam, reference file reads). The GIL releases during pysam's C-level I/O operations.

See [Cohort Analysis Guide](cohort-analysis.md) for usage.

## Protocol-based backend system

The library uses Python `Protocol` classes (structural subtyping) to define contracts between the orchestrator and backend implementations. These are declared in `protocols.py`.

**Key protocols:**

- `IntervalIndex` - genomic interval overlap queries
- `FrequencyDatabase` - population frequency lookups
- `ClinVarDatabase` - clinical significance lookups
- `PDFRenderer` - PDF output rendering

Each protocol has two implementations:

1. A pure-Python fallback that works with only `pysam` + `numpy` installed
2. An optimized backend activated when optional extras (`polars`, `pyranges`, `reportlab`) are present

At construction time, `AnnotationEngine` auto-detects available backends. It attempts the optimized backend first; if the import fails or initialization raises, it falls back to the pure-Python implementation.

**Why protocols instead of ABC:**

- No inheritance required; backends only need the correct method signatures.
- Testable with plain mock objects.
- Compatible with `isinstance` checks via `runtime_checkable`.

## Adding a new annotation source

To add a new annotation source (e.g., OMIM):

1. Define the Protocol in `protocols.py`:

    ```python
    class OMIMDatabase(Protocol):
        def load(self, reference_path: Path) -> None: ...
        def lookup_batch(self, variants: list[tuple[str, int, str, str]]) -> list[Optional[str]]: ...
    ```

2. Write the pure-Python implementation in `annotation/omim.py`.

3. Optionally write an accelerated implementation in `annotation/omim_polars.py`.

4. Update `AnnotationEngine.__init__` to build the backend (with fallback logic).

5. Update `AnnotationEngine._annotate_batch` to call the new backend and incorporate results into `AnnotatedVariant`.

6. If the annotation adds a new field, update the `AnnotatedVariant` dataclass in `models/variant.py`.

## Adding a new evidence tag

To add a new ACMG evidence criterion:

1. Add the tag to the `EvidenceTag` enum in `models/variant.py`.
2. Add its strength to `EVIDENCE_STRENGTH_MAP`.
3. Add an `_evaluate_*` method to `ACMGClassifier` in `classification/acmg.py`.
4. Call it from `_assign_tags`.
5. Update combining rules in `classification/combining.py` if the new strength tier changes thresholds.

## Cache layer

`_internal/cache.py` provides transparent pickle caching for expensive reference file parsing. The flow:

1. Caller asks for cached data via `load_cached(source_path, parser_fn, version)`.
2. Cache checks for a `.vartriage.cache` file adjacent to the source.
3. If the cache exists, has a matching version stamp, and the source mtime hasn't changed, it deserializes and returns.
4. Otherwise, it calls `parser_fn`, serializes the result, writes atomically (temp file + `os.rename`), and returns.

Currently cached: GTF interval trees, CADD score dicts, REVEL score dicts, SpliceAI score dicts.

Cache invalidation triggers:

- Source file mtime changes (you downloaded a new GENCODE release, etc.)
- vartriage version changes (internal data structures may have changed)
- Cache file deleted manually

The atomic write ensures a crash during serialization won't leave a corrupt cache.

## Frequency backend routing

`AnnotationEngine` selects the frequency backend based on the file extension of `gnomad_path`:

| Extension | Backend | Memory profile |
| --- | --- | --- |
| `.vcf.bgz`, `.vcf.gz` | `TabixFrequencyDatabase` | Near-zero (queries on the fly) |
| `.tsv`, `.tsv.gz` | `PolarsFrequencyDatabase` or `DictFrequencyDatabase` | Full dict in RAM |

`TabixFrequencyDatabase` uses pysam's `TabixFile` to query the bgzipped VCF by region. Each batch of variants is looked up individually by `chrom:pos`. This is slower per-query than the dict backend but uses no memory for the reference file, making it practical for full gnomAD (17 GB compressed).

The dict/polars backends remain the better choice when the file fits in RAM and you're annotating large batches.

## Scoring model (v0.8.0+)

Two scoring paths coexist:

- **Composite rank** (legacy, v0.1.0): weighted average of normalized CADD and REVEL. When SpliceAI is available (v0.4.0+), weights are 0.5/0.3/0.2 (REVEL/CADD/SpliceAI). Falls back to proportional redistribution when scores are missing.
- **Prioritization score** (v0.8.0+): literature-validated scoring. Uses REVEL directly for missense (0.7 threshold validated against ClinGen data), SpliceAI for splice-adjacent variants, CADD Phred/60 for non-missense. The recommended ranking method.

Both are present in output for backward compatibility. `composite_rank` will be removed in v1.0.0.

## Memory model

The pipeline processes data in batches (default 10,000 variants). At any time, only one batch is fully materialized in memory. The final `ClassifiedVariant` list is accumulated for report generation, but this is bounded by the number of variants that pass filtering (typically <1% of a WGS VCF).

Peak RSS for a 4M-variant WGS file stays under 2GB with default settings. The chunked fallback in `PrioritizationEngine` handles memory pressure by reducing batch sizes on `MemoryError` (cap: 500K per chunk).
