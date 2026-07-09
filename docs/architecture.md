# Architecture

Package structure, the Protocol-based backend system, and extension points.

## Package structure

```text
vartriage/
├── __init__.py              # Public API exports
├── pipeline.py              # Top-level orchestrator
├── protocols.py             # Protocol interfaces for swappable backends
├── io/
│   ├── vcf_parser.py       # pysam-based VCF streaming
│   └── exceptions.py       # Error hierarchy
├── filter/
│   └── quality_filter.py   # FILTER + QUAL gate
├── annotation/
│   ├── engine.py            # Orchestrator with backend auto-detection
│   ├── consequence.py       # Pure-Python interval tree consequence lookup
│   ├── consequence_pyranges.py  # PyRanges backend (optional)
│   ├── frequency.py         # Dict-based gnomAD lookup
│   ├── frequency_polars.py  # Polars backend (optional)
│   ├── clinvar.py           # Dict-based ClinVar lookup
│   └── clinvar_polars.py    # Polars backend (optional)
├── prioritization/
│   ├── engine.py            # AF gating + scoring orchestrator
│   ├── frequency_filter.py  # Allele frequency threshold filter
│   └── scoring.py           # CADD/REVEL normalization + composite
├── classification/
│   ├── acmg.py              # Evidence tag assignment
│   └── combining.py         # ACMG combining rules
├── reporting/
│   ├── generator.py         # Format routing + atomic write
│   ├── json_writer.py       # JSON serialization
│   ├── csv_writer.py        # CSV serialization
│   ├── pdf_writer.py        # ReportLab PDF (optional)
│   └── pdf_fallback.py      # Fallback when reportlab not installed
├── models/
│   ├── config.py            # All config dataclasses
│   ├── variant.py           # Variant, AnnotatedVariant, ScoredVariant, ClassifiedVariant, enums
│   └── warnings.py          # MissingDataWarning
└── _internal/
    ├── batch.py             # Batch iteration utilities
    ├── interval_tree.py     # Sorted-array interval tree
    ├── vectorized.py        # NumPy vectorized operations
    └── warning_accumulator.py  # Warning collection + threshold
```

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

To add a new annotation source (e.g., OMIM, SpliceAI):

1. Define the Protocol in `protocols.py`:

```python
class SpliceAIDatabase(Protocol):
    def load(self, reference_path: Path) -> None: ...
    def lookup_batch(self, variants: list[tuple[str, int, str, str]]) -> list[Optional[float]]: ...
```

2. Write the pure-Python implementation in `annotation/spliceai.py`.

3. Optionally write an accelerated implementation in `annotation/spliceai_polars.py`.

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

## Memory model

The pipeline processes data in batches (default 10,000 variants). At any time, only one batch is fully materialized in memory. The final `ClassifiedVariant` list is accumulated for report generation, but this is bounded by the number of variants that pass filtering (typically <1% of a WGS VCF).

Peak RSS for a 4M-variant WGS file stays under 2GB with default settings.
