"""Hypothesis tests for missing data handling.

Covers MissingDataWarning completeness, summary warning emission when
thresholds are exceeded, and partial multi-source frequency resolution.
"""

from __future__ import annotations

import tempfile
import warnings as python_warnings
from pathlib import Path

from hypothesis import given, settings, assume
from hypothesis import strategies as st

from vartriage.annotation.frequency import DictFrequencyDatabase
from vartriage.annotation.clinvar import DictClinVarDatabase
from vartriage.annotation.engine import AnnotationEngine
from vartriage._internal.warning_accumulator import (
    MissingDataSummaryWarning,
    WarningAccumulator,
)
from vartriage.models.config import (
    AnnotationConfig,
    MissingDataConfig,
)
from vartriage.models.variant import (
    ClinVarAssertion,
    Variant,
)
from vartriage.models.warnings import MissingDataWarning

from tests.generators.variants import (
    chromosome,
    genomic_position,
    snv_allele,
)

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

@st.composite
def variant_key(draw: st.DrawFn) -> tuple[str, int, str, str]:
    """Generate a (chrom, pos, ref, alt) lookup key."""
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())
    return (chrom, pos, ref, alt)

@st.composite
def missing_data_scenario(
    draw: st.DrawFn,
) -> tuple[
    dict[tuple[str, int, str, str], float],
    list[tuple[str, int, str, str]],
]:
    """Generate a frequency database and query keys with guaranteed misses.

    Returns a dict of known entries and a list of query keys where at least
    one key is guaranteed to be absent from the database.
    """
    num_entries = draw(st.integers(min_value=1, max_value=15))
    db_entries: dict[tuple[str, int, str, str], float] = {}

    for _ in range(num_entries):
        key = draw(variant_key())
        freq = draw(
            st.floats(
                min_value=0.0,
                max_value=1.0,
                allow_nan=False,
                allow_infinity=False,
            )
        )
        db_entries[key] = freq

    # Generate keys that are NOT in the database
    absent_keys: list[tuple[str, int, str, str]] = []
    for _ in range(draw(st.integers(min_value=1, max_value=5))):
        key = draw(variant_key())
        if key not in db_entries:
            absent_keys.append(key)

    assume(len(absent_keys) >= 1)
    return db_entries, absent_keys

@st.composite
def multi_source_scenario(
    draw: st.DrawFn,
) -> tuple[
    dict[tuple[str, int, str, str], float],
    dict[tuple[str, int, str, str], ClinVarAssertion],
    list[tuple[str, int, str, str]],
]:
    """Generate gnomAD and ClinVar databases with overlapping keys.

    Produces a set of variant keys where gnomAD has some data but ClinVar
    does not (or vice versa), simulating partial multi-source resolution.
    """
    num_variants = draw(st.integers(min_value=2, max_value=10))
    all_keys: list[tuple[str, int, str, str]] = []
    for _ in range(num_variants):
        all_keys.append(draw(variant_key()))

    # Deduplicate
    all_keys = list(set(all_keys))
    assume(len(all_keys) >= 2)

    # gnomAD has some keys but not all
    gnomad_count = draw(
        st.integers(min_value=1, max_value=max(1, len(all_keys) - 1))
    )
    gnomad_keys = all_keys[:gnomad_count]
    gnomad_db: dict[tuple[str, int, str, str], float] = {}
    for key in gnomad_keys:
        gnomad_db[key] = draw(
            st.floats(
                min_value=0.0,
                max_value=1.0,
                allow_nan=False,
                allow_infinity=False,
            )
        )

    # ClinVar has some keys, possibly overlapping with gnomAD
    clinvar_count = draw(
        st.integers(min_value=1, max_value=max(1, len(all_keys) - 1))
    )
    clinvar_keys = all_keys[-clinvar_count:]
    clinvar_db: dict[tuple[str, int, str, str], ClinVarAssertion] = {}
    for key in clinvar_keys:
        clinvar_db[key] = draw(
            st.sampled_from(list(ClinVarAssertion))
        )

    return gnomad_db, clinvar_db, all_keys

# ---------------------------------------------------------------------------

@given(scenario=missing_data_scenario())
@settings(max_examples=100)
def test_missing_data_warning_has_all_required_fields(
    scenario: tuple[
        dict[tuple[str, int, str, str], float],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """Every MissingDataWarning contains chrom, pos, ref, alt, and source.

    When a variant is not found in the gnomAD reference database, the
    emitted MissingDataWarning must have all five required fields populated
    with non-None values that match the queried variant coordinates.
    """
    db_entries, absent_keys = scenario

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as f:
        f.write("chrom\tpos\tref\talt\taf\n")
        for (chrom, pos, ref, alt), freq in db_entries.items():
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq}\n")
        tmp_path = Path(f.name)

    try:
        db = DictFrequencyDatabase()
        db.load(tmp_path)
        db.warnings.clear()

        db.lookup_batch(absent_keys)

        assert len(db.warnings) == len(absent_keys), (
            f"Expected {len(absent_keys)} warnings, "
            f"got {len(db.warnings)}"
        )

        for warning, key in zip(db.warnings, absent_keys):
            # Every warning must have all required fields populated
            assert warning.chrom is not None and warning.chrom != ""
            assert warning.pos is not None and warning.pos > 0
            assert warning.ref is not None and warning.ref != ""
            assert warning.alt is not None and warning.alt != ""
            assert warning.source is not None and warning.source != ""

            # Fields must match the queried variant
            assert warning.chrom == key[0], (
                f"Warning chrom={warning.chrom} != queried {key[0]}"
            )
            assert warning.pos == key[1], (
                f"Warning pos={warning.pos} != queried {key[1]}"
            )
            assert warning.ref == key[2], (
                f"Warning ref={warning.ref} != queried {key[2]}"
            )
            assert warning.alt == key[3], (
                f"Warning alt={warning.alt} != queried {key[3]}"
            )
            assert warning.source == "gnomAD", (
                f"Warning source={warning.source}, expected 'gnomAD'"
            )
    finally:
        tmp_path.unlink(missing_ok=True)

@given(scenario=missing_data_scenario())
@settings(max_examples=100)
def test_missing_data_warning_from_annotation_engine(
    scenario: tuple[
        dict[tuple[str, int, str, str], float],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """AnnotationEngine emits complete MissingDataWarnings for missing data.

    When the annotation engine processes variants not found in gnomAD,
    each emitted warning has the full set of required fields: chrom,
    pos, ref, alt, and source.
    """
    db_entries, absent_keys = scenario

    # Build a minimal GTF file for the annotation engine
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".gtf", delete=False
    ) as gtf_f:
        gtf_f.write("# GTF file\n")
        gtf_path = Path(gtf_f.name)

    # Build gnomAD TSV
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as gnomad_f:
        gnomad_f.write("chrom\tpos\tref\talt\taf\n")
        for (chrom, pos, ref, alt), freq in db_entries.items():
            gnomad_f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq}\n")
        gnomad_path = Path(gnomad_f.name)

    try:
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=None,
            batch_size=1000,
        )
        engine = AnnotationEngine(config)

        # Create Variant objects from absent keys
        variants = [
            Variant(
                chrom=k[0],
                pos=k[1],
                id=None,
                ref=k[2],
                alt=k[3],
                qual=30.0,
                filter_status="PASS",
                info={},
            )
            for k in absent_keys
        ]

        # Run annotation
        list(engine.annotate(iter(variants)))

        # Each absent variant should have produced a gnomAD warning
        gnomad_warnings = [
            w for w in engine.warnings if w.source == "gnomAD"
        ]
        assert len(gnomad_warnings) == len(absent_keys)

        for warning in gnomad_warnings:
            assert warning.chrom is not None and warning.chrom != ""
            assert warning.pos is not None and warning.pos > 0
            assert warning.ref is not None and warning.ref != ""
            assert warning.alt is not None and warning.alt != ""
            assert warning.source == "gnomAD"
    finally:
        gtf_path.unlink(missing_ok=True)
        gnomad_path.unlink(missing_ok=True)

# ---------------------------------------------------------------------------

@given(
    threshold=st.integers(min_value=1, max_value=20),
    num_warnings=st.integers(min_value=1, max_value=50),
)
@settings(max_examples=100)
def test_summary_warning_emitted_when_threshold_exceeded(
    threshold: int,
    num_warnings: int,
) -> None:
    """Summary warning emitted when MissingDataWarning count exceeds threshold.

    When the number of missing data warnings in a pipeline run exceeds the
    configured threshold, the WarningAccumulator emits a summary containing
    the total count and the names of contributing reference sources.
    """
    config = MissingDataConfig(warning_threshold=threshold)
    accumulator = WarningAccumulator(config)

    with python_warnings.catch_warnings(record=True) as caught:
        python_warnings.simplefilter("always")

        for i in range(num_warnings):
            source = "gnomAD" if i % 2 == 0 else "ClinVar"
            accumulator.add(
                MissingDataWarning(
                    chrom=f"chr{(i % 22) + 1}",
                    pos=1000 + i,
                    ref="A",
                    alt="T",
                    source=source,
                    reason="not_found",
                )
            )

    exceeds_threshold = num_warnings > threshold

    if exceeds_threshold:
        assert accumulator.threshold_exceeded is True
        assert accumulator.summary_emitted is True

        # A summary warning should have been emitted via Python warnings
        summary_warnings = [
            w for w in caught if issubclass(w.category, UserWarning)
        ]
        assert len(summary_warnings) >= 1

        summary = accumulator.build_summary()
        assert summary.total_count == num_warnings
        assert len(summary.sources) > 0
        assert "gnomAD" in summary.sources
    else:
        assert accumulator.threshold_exceeded is False
        assert accumulator.summary_emitted is False

@given(
    threshold=st.integers(min_value=1, max_value=10),
    extra_count=st.integers(min_value=1, max_value=20),
)
@settings(max_examples=100)
def test_summary_warning_contains_total_count_and_sources(
    threshold: int,
    extra_count: int,
) -> None:
    """Summary warning includes the total count and all contributing sources.

    When a summary is triggered, the total count equals the actual number
    of warnings emitted, and the source list includes every distinct
    reference source that contributed.
    """
    num_warnings = threshold + extra_count  # always exceeds
    config = MissingDataConfig(warning_threshold=threshold)
    accumulator = WarningAccumulator(config)

    sources = ["gnomAD", "ClinVar"]

    with python_warnings.catch_warnings(record=True):
        python_warnings.simplefilter("always")

        for i in range(num_warnings):
            source = sources[i % len(sources)]
            accumulator.add(
                MissingDataWarning(
                    chrom=f"chr{(i % 22) + 1}",
                    pos=1000 + i,
                    ref="A",
                    alt="T",
                    source=source,
                    reason="not_found",
                )
            )

    summary = accumulator.build_summary()
    assert summary.total_count == num_warnings
    assert summary.total_count > threshold
    assert "gnomAD" in summary.sources
    assert "ClinVar" in summary.sources

@given(
    threshold=st.integers(min_value=5, max_value=50),
    num_warnings=st.integers(min_value=1, max_value=4),
)
@settings(max_examples=100)
def test_no_summary_when_below_threshold(
    threshold: int,
    num_warnings: int,
) -> None:
    """No summary warning when count is at or below the threshold.
    """
    assume(num_warnings <= threshold)

    config = MissingDataConfig(warning_threshold=threshold)
    accumulator = WarningAccumulator(config)

    with python_warnings.catch_warnings(record=True) as caught:
        python_warnings.simplefilter("always")

        for i in range(num_warnings):
            accumulator.add(
                MissingDataWarning(
                    chrom=f"chr{(i % 22) + 1}",
                    pos=1000 + i,
                    ref="A",
                    alt="T",
                    source="gnomAD",
                    reason="not_found",
                )
            )

    # Below threshold: no summary should be emitted
    assert accumulator.threshold_exceeded is False
    assert accumulator.summary_emitted is False

    summary_warnings = [
        w for w in caught if issubclass(w.category, UserWarning)
    ]
    assert len(summary_warnings) == 0

# ---------------------------------------------------------------------------

@given(scenario=multi_source_scenario())
@settings(max_examples=100)
def test_partial_multi_source_uses_available_frequency(
    scenario: tuple[
        dict[tuple[str, int, str, str], float],
        dict[tuple[str, int, str, str], ClinVarAssertion],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """Uses available frequency value when one source has data and another doesn't.

    When a variant is queried against multiple reference sources and at
    least one source returns a valid frequency while others return no data,
    the pipeline uses the available frequency value for the variant and
    emits MissingDataWarnings only for the sources that returned no data.
    """
    gnomad_db, clinvar_db, all_keys = scenario

    # Build gnomAD TSV
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as gnomad_f:
        gnomad_f.write("chrom\tpos\tref\talt\taf\n")
        for (chrom, pos, ref, alt), freq in gnomad_db.items():
            gnomad_f.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq}\n"
            )
        gnomad_path = Path(gnomad_f.name)

    # Build ClinVar TSV
    assertion_to_str = {
        ClinVarAssertion.PATHOGENIC: "Pathogenic",
        ClinVarAssertion.LIKELY_PATHOGENIC: "Likely pathogenic",
        ClinVarAssertion.VUS: "Uncertain significance",
        ClinVarAssertion.LIKELY_BENIGN: "Likely benign",
        ClinVarAssertion.BENIGN: "Benign",
    }

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as clinvar_f:
        clinvar_f.write(
            "chrom\tpos\tref\talt\tclinical_significance\n"
        )
        for (chrom, pos, ref, alt), assertion in clinvar_db.items():
            sig_str = assertion_to_str[assertion]
            clinvar_f.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{sig_str}\n"
            )
        clinvar_path = Path(clinvar_f.name)

    # Build minimal GTF
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".gtf", delete=False
    ) as gtf_f:
        gtf_f.write("# GTF file\n")
        gtf_path = Path(gtf_f.name)

    try:
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=clinvar_path,
            batch_size=1000,
        )
        engine = AnnotationEngine(config)

        # Create Variant objects for all keys
        variants = [
            Variant(
                chrom=k[0],
                pos=k[1],
                id=None,
                ref=k[2],
                alt=k[3],
                qual=30.0,
                filter_status="PASS",
                info={},
            )
            for k in all_keys
        ]

        results = list(engine.annotate(iter(variants)))

        # Verify: for keys in gnomAD, frequency is assigned
        for i, key in enumerate(all_keys):
            result = results[i]
            if key in gnomad_db:
                assert result.allele_frequency is not None, (
                    f"Key {key} is in gnomAD but got None frequency"
                )
                assert abs(
                    result.allele_frequency - gnomad_db[key]
                ) < 1e-7
                assert result.frequency_unknown is False
            else:
                assert result.allele_frequency is None
                assert result.frequency_unknown is True

        # Verify: warnings only emitted for sources that returned no data
        gnomad_warnings = [
            w for w in engine.warnings if w.source == "gnomAD"
        ]
        clinvar_warnings = [
            w for w in engine.warnings if w.source == "ClinVar"
        ]

        # gnomAD warnings: only for keys NOT in gnomad_db
        expected_gnomad_misses = [
            k for k in all_keys if k not in gnomad_db
        ]
        assert len(gnomad_warnings) == len(expected_gnomad_misses), (
            f"Expected {len(expected_gnomad_misses)} gnomAD warnings, "
            f"got {len(gnomad_warnings)}"
        )

        # ClinVar warnings: only for keys NOT in clinvar_db
        expected_clinvar_misses = [
            k for k in all_keys if k not in clinvar_db
        ]
        assert len(clinvar_warnings) == len(expected_clinvar_misses), (
            f"Expected {len(expected_clinvar_misses)} ClinVar warnings, "
            f"got {len(clinvar_warnings)}"
        )

        # Every gnomAD warning matches an absent key
        for warning in gnomad_warnings:
            w_key = (warning.chrom, warning.pos, warning.ref, warning.alt)
            assert w_key not in gnomad_db, (
                f"Got gnomAD warning for key {w_key} that IS in the DB"
            )

        # Every ClinVar warning matches an absent key
        for warning in clinvar_warnings:
            w_key = (warning.chrom, warning.pos, warning.ref, warning.alt)
            assert w_key not in clinvar_db, (
                f"Got ClinVar warning for key {w_key} that IS in the DB"
            )
    finally:
        gnomad_path.unlink(missing_ok=True)
        clinvar_path.unlink(missing_ok=True)
        gtf_path.unlink(missing_ok=True)

@given(scenario=multi_source_scenario())
@settings(max_examples=100)
def test_partial_resolution_no_warning_for_present_sources(
    scenario: tuple[
        dict[tuple[str, int, str, str], float],
        dict[tuple[str, int, str, str], ClinVarAssertion],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """No warning emitted for reference sources that return valid data.

    When gnomAD returns a frequency for a variant, no gnomAD
    MissingDataWarning should exist for that variant. Same principle
    applies to ClinVar.
    """
    gnomad_db, clinvar_db, all_keys = scenario

    # Build gnomAD TSV
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as gnomad_f:
        gnomad_f.write("chrom\tpos\tref\talt\taf\n")
        for (chrom, pos, ref, alt), freq in gnomad_db.items():
            gnomad_f.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq}\n"
            )
        gnomad_path = Path(gnomad_f.name)

    # Build ClinVar TSV
    assertion_to_str = {
        ClinVarAssertion.PATHOGENIC: "Pathogenic",
        ClinVarAssertion.LIKELY_PATHOGENIC: "Likely pathogenic",
        ClinVarAssertion.VUS: "Uncertain significance",
        ClinVarAssertion.LIKELY_BENIGN: "Likely benign",
        ClinVarAssertion.BENIGN: "Benign",
    }

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as clinvar_f:
        clinvar_f.write(
            "chrom\tpos\tref\talt\tclinical_significance\n"
        )
        for (chrom, pos, ref, alt), assertion in clinvar_db.items():
            sig_str = assertion_to_str[assertion]
            clinvar_f.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{sig_str}\n"
            )
        clinvar_path = Path(clinvar_f.name)

    # Build minimal GTF
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".gtf", delete=False
    ) as gtf_f:
        gtf_f.write("# GTF file\n")
        gtf_path = Path(gtf_f.name)

    try:
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=clinvar_path,
            batch_size=1000,
        )
        engine = AnnotationEngine(config)

        variants = [
            Variant(
                chrom=k[0],
                pos=k[1],
                id=None,
                ref=k[2],
                alt=k[3],
                qual=30.0,
                filter_status="PASS",
                info={},
            )
            for k in all_keys
        ]

        list(engine.annotate(iter(variants)))

        # For keys present in gnomAD, no gnomAD warning should exist
        for key in all_keys:
            if key in gnomad_db:
                gnomad_warnings_for_key = [
                    w
                    for w in engine.warnings
                    if w.source == "gnomAD"
                    and (w.chrom, w.pos, w.ref, w.alt) == key
                ]
                assert len(gnomad_warnings_for_key) == 0, (
                    f"Unexpected gnomAD warning for present key {key}"
                )

            if key in clinvar_db:
                clinvar_warnings_for_key = [
                    w
                    for w in engine.warnings
                    if w.source == "ClinVar"
                    and (w.chrom, w.pos, w.ref, w.alt) == key
                ]
                assert len(clinvar_warnings_for_key) == 0, (
                    f"Unexpected ClinVar warning for present key {key}"
                )
    finally:
        gnomad_path.unlink(missing_ok=True)
        clinvar_path.unlink(missing_ok=True)
        gtf_path.unlink(missing_ok=True)
