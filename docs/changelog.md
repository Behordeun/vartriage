# Changelog

All notable changes to vartriage are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.10.0] - 2026-07-18

### Added

- **Zygosity model**: `Zygosity` enum (HETEROZYGOUS, HOMOZYGOUS_ALT, HEMIZYGOUS, UNKNOWN) and `zygosity` field on `AnnotatedVariant`. Ready for population from VCF FORMAT GT field.
- **Variant quality metrics**: `VariantQualityMetrics` dataclass (depth, genotype_quality, allele_balance, is_low_confidence) and `quality_metrics` field on `AnnotatedVariant`.
- **ACMG Secondary Findings (SF v3.2)**: shipped 71-gene list as package data (`vartriage/data/acmg_sf_v3.2.txt`). `SecondaryFindingsFilter` class with `is_secondary_finding()` and `split_stream()`. CLI flag `--secondary-findings`.
- **Computational-only disclaimer**: clinical reports now display a prominent banner after the header citing ACMG/AMP 2015 (Richards et al.) and stating findings require clinical geneticist review.

### Notes

- Zygosity extraction from VCF FORMAT fields and quality metrics population planned for v0.10.1 (requires VCFParser changes for FORMAT field access).
- HGVS nomenclature generation planned for v0.10.1 (requires CodonContext pipeline integration).

## [0.9.0] - 2026-07-17

### Added

- **Benign evidence criteria**: the pipeline can now classify variants as Benign or Likely Benign for the first time.
  - BA1: any gnomAD population AF > 5% = standalone Benign
  - BS1: any population AF > 1% = strong benign evidence
  - BP4: low computational predictor scores (REVEL < 0.15 for missense, CADD Phred < 10 for non-missense)
  - BP7: synonymous variant with SpliceAI < 0.1 (no splice impact)
- **Population-specific frequency model** (`PopulationFrequencies` dataclass): stores per-population gnomAD AFs (AFR, AMR, ASJ, EAS, FIN, NFE, SAS) with `max_population_af`, `any_exceeds()`, and `all_below()` helpers.
- **Updated combining rules** (ACMG Table 5): BA1 standalone = Benign, 2 BS = Benign, 1 BS + 1 BP = Likely Benign, conflicting pathogenic + benign = VUS.
- `STANDALONE` strength tier in `EvidenceStrength` enum for BA1.
- 9 new `EvidenceTag` enum members: PS1, PM1, PM4, PM5 (pathogenic, evaluators planned for v0.9.1), BA1, BS1, BS2, BP4, BP7 (benign, evaluators active).
- `has_conflicting_evidence()` helper in combining module.
- `population_frequencies` field on `AnnotatedVariant`.
- 24 unit tests for benign criteria covering positive/negative/edge cases.

### Changed

- PM2 now checks ALL population-specific frequencies below threshold (not just global AF). Falls back to global AF when per-population data is absent.
- `EVIDENCE_STRENGTH_MAP` expanded from 4 to 13 entries covering all new tags.
- `combine_evidence()` rewritten: separates pathogenic and benign evidence, detects conflicts, applies benign combining rules.
- README links converted to absolute GitHub URLs (fixes dead links on PyPI).

### Notes

- PS1, PM1, PM4, PM5 evaluators are not yet implemented (tags exist in the enum, evaluators planned for v0.9.1 when ClinVar amino acid index and functional domain data are available).
- BS2 tag exists but evaluator requires gnomAD homozygote count data not currently parsed (planned for v0.9.1).

## [0.8.0] - 2026-07-16

### Added

- **Codon-level consequence calling** (`--reference-fasta`): SNVs in CDS regions now use actual amino acid comparison instead of the positional heuristic. Requires an indexed reference FASTA. Correctly distinguishes synonymous, missense, and nonsense changes by extracting the reference codon, substituting the variant base, and translating both codons.
- **Variant normalization**: left-align and trim indels before gnomAD/ClinVar/score lookups using the reference FASTA (Tan et al. 2015 algorithm). Reduces silent lookup failures caused by representation differences between VCF callers and reference databases.
- **Prioritization score** (`prioritization_score` output field): literature-backed scoring using REVEL directly for missense (validated 0.7 threshold), SpliceAI for splice-adjacent, CADD Phred/60 for non-missense. Replaces the unvalidated 0.4/0.6 weighted average as the recommended ranking method.
- `TranscriptCDSIndex`: per-transcript CDS exon map built from GTF, enabling genomic-to-CDS position mapping for forward and reverse strand genes.
- `CodonResolver`: FASTA-backed codon extraction with split-codon-at-exon-junction handling and negative-strand reverse complement.
- `VariantNormalizer`: three-step normalization (right-trim, left-trim, left-align) with 1000-iteration safety cap.
- `compute_prioritization_score()` function in scoring module.
- Standard genetic code table (`_internal/genetic_code.py`) with `translate_codon()` and `reverse_complement()`.
- 36 unit tests covering genetic code, transcript index, normalizer, and prioritization score.
- 3 live VEP concordance integration tests (marked `@pytest.mark.slow`).

### Changed

- `ScoredVariant` now carries both `composite_rank` (legacy) and `prioritization_score` (new). Both synced via `__post_init__` when only one is set.
- JSON output includes `prioritization_score` field alongside `composite_rank`.
- `_determine_consequence()` accepts optional `codon_resolver` parameter; uses it for SNVs in CDS when available, falls back to positional heuristic without FASTA.
- `AnnotationConfig` gains `reference_fasta_path: Optional[Path]` field.
- GTF parsing now populates a `TranscriptCDSIndex` from CDS features (frame column extracted).

### Notes

- Without `--reference-fasta`, all behavior is unchanged from v0.7.0 (backward compatible).
- The positional heuristic ("SNV in CDS = Missense") remains as a fallback. To get correct consequence calling, provide the reference FASTA.
- `composite_rank` is deprecated in favor of `prioritization_score`. Both are present in output for backward compatibility. `composite_rank` will be removed in v1.0.0.

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
