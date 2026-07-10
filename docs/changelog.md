# Changelog

All notable changes to vartriage are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.1] - 2025-07-10

### Fixed

- Frequency loader now treats `'.'` in the af column as null instead of raising a parse error (gnomAD files use `'.'` for missing values)

### Added

- CLI entry point: `vartriage --vcf ... --output ...` with full argument parsing
- Streaming report generation — JSON and CSV write directly from iterators, no full-variant buffering
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
