# Changelog

## v0.1.0

Initial release.

- VCF parsing via pysam with streaming iteration
- Quality filtering (FILTER field + QUAL score threshold)
- Annotation engine with auto-detected backends:
    - Functional consequence via GTF/GFF gene models
    - Population frequency via gnomAD
    - Clinical significance via ClinVar
- Prioritization with allele frequency gate and composite CADD/REVEL scoring
- ACMG/AMP 2015 evidence classification (PVS1, PM2, PP3, PP5)
- Report generation in JSON, CSV, and PDF formats
- Protocol-based backend system with pure-Python fallbacks
- Memory-bounded processing for whole-genome scale data
