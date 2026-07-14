# Changelog

All notable changes to vartriage are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.7.0] - 2026-07-14

### Added

- **API annotation mode** (`--mode api|hybrid`): query Ensembl VEP, ClinVar E-utilities, CADD, and SpliceAI via HTTP instead of local reference files. Zero-config variant triage for gene panels and exploratory analysis.
- **Ensembl VEP client**: batch POST annotation (200 variants/request) with consequence, gene symbol, gnomAD frequency, and CADD Phred extraction. VEP Sequence Ontology terms mapped to FunctionalConsequence enum across 14 severity ranks. GRCh37 and GRCh38 support.
- **ClinVar E-utilities client**: esearch + esummary lookups with clinical significance mapping, review status ranking for conflicting interpretations, and NCBI API key support for higher rate limits.
- **CADD REST API client**: position-based Phred score lookups with ref/alt allele filtering. Used as fallback when VEP's CADD plugin has no pre-computed score.
- **SpliceAI Lookup client**: max delta score extraction with consequence-based smart filtering (only queries splice-relevant variants to conserve rate limit).
- **API annotation engine**: composes VEP + ClinVar into the same `annotate()` interface as the local engine. Concurrent ClinVar lookups via ThreadPoolExecutor.
- **API score provider**: CADD score hierarchy (VEP plugin first, standalone API fallback) and SpliceAI lookups. Documents REVEL limitation (no API exists).
- **Response caching**: SQLite-backed at `~/.vartriage/api_cache.db` with configurable TTL (7 days default, 30 days for ClinVar). Pinned mode (`cache_ttl_days = -1`) for clinical reproducibility.
- **Resilience stack**: token bucket rate limiter (per-service with daily caps), circuit breaker (CLOSED/OPEN/HALF_OPEN with 60s recovery), exponential backoff retry (1s/2s/4s, max 3), Retry-After header parsing.
- **VCF-to-VEP notation converter**: handles SNVs, deletions, insertions, MNVs, complex indels, and chr prefix stripping.
- CLI flags: `--mode` (local/api/hybrid), `--api-key` (NCBI), `--no-confirm` (skip large-run prompt).
- HTTP proxy support via `HTTPS_PROXY`/`HTTP_PROXY` env vars and `[api.proxy]` TOML config.
- User-Agent header on all requests identifying vartriage version and project URL.
- `docs/api-mode.md` user guide.
- `pip install vartriage[api]` optional dependency group (httpx).

### Changed

- `PipelineConfig` accepts an optional `api` field for API backend configuration.
- `Pipeline.run()` routes to `APIAnnotationEngine` when API mode is active.
- `pyproject.toml`: added `api` and updated `all` optional dependency groups.

## [0.6.0] - 2026-07-13

### Added

- **Score bundle downloader** (`vartriage bundle`): automated downloading, transformation, and management of reference files. Supported bundles: ClinVar, gnomAD exomes (chr22), REVEL, GENCODE, SpliceAI. Subcommands: `download`, `list`, `verify`, `status`, `update-registry`.
- **Bundle-aware pipeline** (`--use-bundles`, `--genome-build`): auto-resolve reference file paths from installed bundles. Explicit CLI paths take precedence. No configuration needed after `vartriage bundle download`.
- HTTP download engine with resume support (Range header), atomic `.partial` files, exponential backoff retry (1s, 2s, 4s) for 429/5xx errors, and streaming progress bar on stderr.
- Post-download transformers: VCF-to-TSV (bcftools with pysam fallback), ClinVar normalizer (adds `chr` prefix + CLNSIG mapping), CSV-to-TSV (REVEL), SpliceAI max-delta extractor, and passthrough (gz decompression).
- Per-bundle manifests tracking version, checksums, download timestamps, and source URLs.
- Configurable storage layout at `~/.vartriage/bundles/` with env var override (`VARTRIAGE_BUNDLE_STORAGE`).
- TOML configuration file support (`~/.vartriage/config.toml`) for default build, concurrency, and proxy settings.
- Disk space pre-flight validation before downloads.
- SHA-256 checksum verification for downloaded and transformed files.
- `use_bundles` and `genome_build` fields on `PipelineConfig`.
- `docs/bundles.md` user guide with command reference and storage layout documentation.

### Changed

- `pyproject.toml`: added `tomli` and `tomllib` to mypy ignore list for Python 3.10 compatibility.

## [0.5.0] - 2026-07-13

### Added

- **Clinical report generation** (`--output-format clinical-html|clinical-pdf|clinical-docx`): produces structured, sign-off-ready clinical variant reports from ClassifiedVariant data. Per-variant evidence narratives with plain-language ACMG criteria explanations. Structured sections: Header, Executive Summary, Findings Table, Evidence Cards, Limitations, Methodology, Sign-off. JSON audit trail sidecar (`.audit.json`) with run manifest and per-variant decision log.
- CLI flags `--patient-id` and `--panel-name` for clinical report metadata. Required when using any `clinical-*` output format.
- `ClinicalReportConfig` dataclass: `patient_id`, `panel_name`, `output_format`, `report_template`.
- Self-contained HTML output (all CSS inlined, no JavaScript, no external dependencies).
- PDF rendering via WeasyPrint (optional dependency). DOCX rendering via python-docx (optional dependency).
- Atomic file writing for clinical reports (temp file + rename).
- Findings table sorted by classification tier (Pathogenic > Likely_Pathogenic > VUS), then by composite rank descending within each tier.
- Empty variant set produces a valid report with zero counts and a negative-finding statement.

### Dependencies

- `weasyprint` (optional, for `clinical-pdf`)
- `python-docx` (optional, for `clinical-docx`)

## [0.4.0] - 2026-07-12

### Added

- **Trio-based inheritance pattern classification** (`--proband`, `--mother`, `--father`): classifies variants into Mendelian inheritance patterns (de novo, dominant, recessive, compound heterozygous, X-linked) based on family genotypes. Multi-label: a single variant can carry multiple patterns when criteria overlap. Compound het uses gene-aware buffering to detect trans inheritance. Replaces SampleExtractor when trio mode is active. Mutually exclusive with `--sample`.
- **SpliceAI score integration** (`--spliceai-scores`): third pathogenicity predictor alongside CADD and REVEL. Dynamic weight redistribution in the composite formula (0.5/0.3/0.2 when all three present, proportional fallback for any two, single-score identity). PP3 now also fires when SpliceAI > 0.5 on splice-site or missense variants. PVS1 fires for SPLICE_SITE + SpliceAI > 0.8. Fully backward-compatible: existing two-score behavior unchanged when SpliceAI is not configured.
- **VCF output format** (`--output-format vcf`): produces an annotated bgzipped VCF (.vcf.gz) with a tabix index (.tbi). Re-reads the source VCF, injects VARTRIAGE_CONSEQUENCE, VARTRIAGE_AF, VARTRIAGE_RANK, VARTRIAGE_ACMG, and VARTRIAGE_TAGS INFO fields for classified variants, and writes all records (matched or not) to output. Directly loadable in IGV and queryable with bcftools.
- **Gene list filtering** (`--gene-list`): restrict analysis to variants in a user-supplied gene panel file. One gene symbol per line, case-insensitive matching. Genes in the list with zero matching variants produce a logged WARNING so you can catch typos or outdated nomenclature.
- `gene_name` field on `AnnotatedVariant` populated during consequence annotation.

## [0.3.0] - 2026-07-12

### Added

- **BED-based region filtering** (`--regions`): target analysis to specific genomic intervals from a BED file.
- **Multi-sample VCF support** (`--sample`, `--min-gq`): extract a single sample from multi-sample VCFs with optional genotype quality filtering.
- `RegionFilterConfig`, `SampleConfig` configuration dataclasses.
- `GeneFilterConfig`, `RegionFilterConfig`, `SampleConfig` configuration dataclasses.
- **VCF output format** (`--output-format vcf`): produces an annotated bgzipped VCF (.vcf.gz) with a tabix index (.tbi). Re-reads the source VCF, injects VARTRIAGE_CONSEQUENCE, VARTRIAGE_AF, VARTRIAGE_RANK, VARTRIAGE_ACMG, and VARTRIAGE_TAGS INFO fields for classified variants, and writes all records (matched or not) to output. Directly loadable in IGV and queryable with bcftools.
- **SpliceAI score integration** (`--spliceai-scores`): third pathogenicity predictor alongside CADD and REVEL. Dynamic weight redistribution in the composite formula (0.5/0.3/0.2 when all three present, proportional fallback for any two, single-score identity). PP3 now also fires when SpliceAI > 0.5 on splice-site or missense variants. PVS1 fires for SPLICE_SITE + SpliceAI > 0.8. Fully backward-compatible: existing two-score behavior unchanged when SpliceAI is not configured.

## [0.2.0] - 2026-07-10

### Added

- **Vectorized batch annotation** in `PyRangesConsequenceAnnotator`: `assign_batch` now performs a single PyRanges join instead of per-variant loops. 680x speedup (7 variants/sec to 4,749 variants/sec).
- **Reference file caching** (`_internal/cache.py`): GTF interval trees, CADD score dicts, and REVEL score dicts are serialized to pickle after first parse. Subsequent runs skip parsing entirely. Cache invalidation uses source file mtime and a version stamp. Writes are atomic (write-to-temp then rename). Cache files sit next to their source with a `.vartriage.cache` suffix.
- **TabixFrequencyDatabase** (`annotation/frequency_tabix.py`): queries bgzipped+tabix-indexed gnomAD VCFs on the fly via pysam. Zero memory footprint for the reference file. Auto-selected when `gnomad_path` ends with `.vcf.bgz` or `.vcf.gz`.
- **Extension-based backend routing** in `AnnotationEngine`: `.vcf.bgz`/`.vcf.gz` paths route to the tabix backend; `.tsv`/`.tsv.gz` paths use the existing polars/dict backends.

### Performance

- chr22 full annotation benchmark (130K variants + GENCODE + 4.8M gnomAD entries): 36.3s wall time, ~2 GB peak RSS
- With 100K gnomAD subset: 19.5s wall time, 453 MB peak RSS

## [0.1.1] - 2026-07-10

### Fixed

- Frequency loader now treats `'.'` in the af column as null instead of raising a parse error (gnomAD files use `'.'` for missing values)

### Added

- CLI entry point: `vartriage --vcf ... --output ...` with full argument parsing
- Streaming report generation: JSON and CSV write directly from iterators, no full-variant buffering
- `ScoreLoader` class for loading CADD/REVEL TSV files into coordinate-keyed dicts
- `VarTriageWarning` base class for the warning hierarchy; `ScoreValidationWarning` and `MissingDataSummaryWarning` inherit from it
- `py.typed` marker for PEP 561 compliance
- Typed Protocol return values in annotation engine (removed all `Any` from return annotations)
- `assign_batch` method on the `IntervalIndex` Protocol
- Upper bounds on all dependency version specs (pysam <1.0, numpy <3.0, polars <2.0, etc.)
- GitHub Actions CI workflow (Python 3.10, 3.11, 3.12)
- Automated PyPI publishing via trusted publisher workflow
- `CONTRIBUTING.md`
- `LICENSE` file (MIT)

## [0.1.0] - 2026-07-09

Initial release on PyPI.

### Added

- VCF parsing via pysam with streaming iteration
- Quality filtering on FILTER field and QUAL score threshold
- Annotation engine with auto-detected backends:
  - Functional consequence via GTF/GFF gene models
  - Population frequency via gnomAD
  - Clinical significance via ClinVar
- Prioritization: allele frequency gate + composite CADD/REVEL scoring
- ACMG/AMP 2015 evidence classification (PVS1, PM2, PP3, PP5)
- Report generation in JSON, CSV, and PDF
- Protocol-based backend system with pure-Python fallbacks
- Memory-bounded processing for whole-genome scale data (4M+ variants under 2 GB RSS)
- Optional extras: `[accelerated]` (polars, pyranges), `[pdf]` (reportlab), `[all]`
