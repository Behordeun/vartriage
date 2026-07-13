# Changelog

All notable changes to vartriage are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.5.0] - Unreleased

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

## [0.4.0] - 2025-07-12

### Added

- **Trio-based inheritance pattern classification** (`--proband`, `--mother`, `--father`): classifies variants into Mendelian inheritance patterns (de novo, dominant, recessive, compound heterozygous, X-linked) based on family genotypes. Multi-label: a single variant can carry multiple patterns when criteria overlap. Compound het uses gene-aware buffering to detect trans inheritance. Replaces SampleExtractor when trio mode is active. Mutually exclusive with `--sample`.
- **SpliceAI score integration** (`--spliceai-scores`): third pathogenicity predictor alongside CADD and REVEL. Dynamic weight redistribution in the composite formula (0.5/0.3/0.2 when all three present, proportional fallback for any two, single-score identity). PP3 now also fires when SpliceAI > 0.5 on splice-site or missense variants. PVS1 fires for SPLICE_SITE + SpliceAI > 0.8. Fully backward-compatible: existing two-score behavior unchanged when SpliceAI is not configured.
- **VCF output format** (`--output-format vcf`): produces an annotated bgzipped VCF (.vcf.gz) with a tabix index (.tbi). Re-reads the source VCF, injects VARTRIAGE_CONSEQUENCE, VARTRIAGE_AF, VARTRIAGE_RANK, VARTRIAGE_ACMG, and VARTRIAGE_TAGS INFO fields for classified variants, and writes all records (matched or not) to output. Directly loadable in IGV and queryable with bcftools.
- **Gene list filtering** (`--gene-list`): restrict analysis to variants in a user-supplied gene panel file. One gene symbol per line, case-insensitive matching. Genes in the list with zero matching variants produce a logged WARNING so you can catch typos or outdated nomenclature.
- `gene_name` field on `AnnotatedVariant` populated during consequence annotation.

## [0.3.0] - 2025-07-12

### Added

- **BED-based region filtering** (`--regions`): target analysis to specific genomic intervals from a BED file.
- **Multi-sample VCF support** (`--sample`, `--min-gq`): extract a single sample from multi-sample VCFs with optional genotype quality filtering.
- `RegionFilterConfig`, `SampleConfig` configuration dataclasses.
- `GeneFilterConfig`, `RegionFilterConfig`, `SampleConfig` configuration dataclasses.
- **VCF output format** (`--output-format vcf`): produces an annotated bgzipped VCF (.vcf.gz) with a tabix index (.tbi). Re-reads the source VCF, injects VARTRIAGE_CONSEQUENCE, VARTRIAGE_AF, VARTRIAGE_RANK, VARTRIAGE_ACMG, and VARTRIAGE_TAGS INFO fields for classified variants, and writes all records (matched or not) to output. Directly loadable in IGV and queryable with bcftools.
- **SpliceAI score integration** (`--spliceai-scores`): third pathogenicity predictor alongside CADD and REVEL. Dynamic weight redistribution in the composite formula (0.5/0.3/0.2 when all three present, proportional fallback for any two, single-score identity). PP3 now also fires when SpliceAI > 0.5 on splice-site or missense variants. PVS1 fires for SPLICE_SITE + SpliceAI > 0.8. Fully backward-compatible: existing two-score behavior unchanged when SpliceAI is not configured.

## [0.2.0] - 2025-07-10

### Added

- **Vectorized batch annotation** in `PyRangesConsequenceAnnotator`: `assign_batch` now performs a single PyRanges join instead of per-variant loops. 680x speedup (7 variants/sec to 4,749 variants/sec).
- **Reference file caching** (`_internal/cache.py`): GTF interval trees, CADD score dicts, and REVEL score dicts are serialized to pickle after first parse. Subsequent runs skip parsing entirely. Cache invalidation uses source file mtime and a version stamp. Writes are atomic (write-to-temp then rename). Cache files sit next to their source with a `.vartriage.cache` suffix.
- **TabixFrequencyDatabase** (`annotation/frequency_tabix.py`): queries bgzipped+tabix-indexed gnomAD VCFs on the fly via pysam. Zero memory footprint for the reference file. Auto-selected when `gnomad_path` ends with `.vcf.bgz` or `.vcf.gz`.
- **Extension-based backend routing** in `AnnotationEngine`: `.vcf.bgz`/`.vcf.gz` paths route to the tabix backend; `.tsv`/`.tsv.gz` paths use the existing polars/dict backends.

### Performance

- chr22 full annotation benchmark (130K variants + GENCODE + 4.8M gnomAD entries): 36.3s wall time, ~2 GB peak RSS
- With 100K gnomAD subset: 19.5s wall time, 453 MB peak RSS

## [0.1.1] - 2025-07-10

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

## [0.1.0] - 2025-07-09

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
