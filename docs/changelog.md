# Changelog

All notable changes to vartriage are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **Trio-based inheritance pattern classification** (`--proband`, `--mother`, `--father`): classifies variants into Mendelian inheritance patterns (de novo, dominant, recessive, compound heterozygous, X-linked) based on family genotypes. Multi-label: a single variant can carry multiple patterns when criteria overlap. Compound het uses gene-aware buffering to detect trans inheritance. Replaces SampleExtractor when trio mode is active. Mutually exclusive with `--sample`.

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
