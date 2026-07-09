"""Hypothesis tests for report serialization (JSON and CSV writers).

Verifies round-trip fidelity for JSON output and structural validity
for CSV output across generated variant data.
"""

from __future__ import annotations

import json
import csv
import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from vartriage.models.variant import (
    ACMGClassification,
    ClassifiedVariant,
    EvidenceTag,
)
from vartriage.reporting.json_writer import write_json
from vartriage.reporting.csv_writer import write_csv, CSV_FIELDS

from tests.generators.variants import scored_variant, evidence_tag_set

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

@st.composite
def classified_variant(draw: st.DrawFn) -> ClassifiedVariant:
    """Generate a ClassifiedVariant from a scored variant with random tags."""
    sv = draw(scored_variant())
    tags = draw(evidence_tag_set())
    classification = draw(st.sampled_from(list(ACMGClassification)))
    missing_sources = draw(
        st.frozensets(
            st.sampled_from(["gnomAD", "REVEL", "ClinVar", "CADD"]),
            min_size=0,
            max_size=3,
        )
    )
    return ClassifiedVariant(
        scored=sv,
        evidence_tags=tags,
        classification=classification,
        missing_data_sources=missing_sources,
    )

@st.composite
def classified_variant_list(
    draw: st.DrawFn, min_size: int = 0, max_size: int = 20
) -> list[ClassifiedVariant]:
    """Generate a list of ClassifiedVariants for serialization testing."""
    return draw(
        st.lists(classified_variant(), min_size=min_size, max_size=max_size)
    )

# ---------------------------------------------------------------------------

@given(variants=classified_variant_list(min_size=0, max_size=15))
@settings(max_examples=100)
def test_json_round_trip_preserves_field_values(
    variants: list[ClassifiedVariant],
) -> None:
    """Serialize to JSON then deserialize produces identical field values, types, and order.

    For every variant in the list, all output fields in the deserialized record
    must match the original ClassifiedVariant field values exactly in value and
    type. The list order is preserved.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "report.json"
        write_json(variants, output_path)

        with open(output_path, "r", encoding="utf-8") as f:
            deserialized = json.load(f)

    assert len(deserialized) == len(variants), (
        f"Expected {len(variants)} records, got {len(deserialized)}"
    )

    for i, (record, original) in enumerate(zip(deserialized, variants)):
        scored = original.scored
        annotated = scored.annotated
        raw = annotated.variant

        assert record["chromosome"] == raw.chrom, (
            f"Record {i}: chromosome mismatch"
        )
        assert record["position"] == raw.pos, (
            f"Record {i}: position mismatch"
        )
        assert isinstance(record["position"], int), (
            f"Record {i}: position should be int, got {type(record['position'])}"
        )
        assert record["ref_allele"] == raw.ref, (
            f"Record {i}: ref_allele mismatch"
        )
        assert record["alt_allele"] == raw.alt, (
            f"Record {i}: alt_allele mismatch"
        )

        expected_consequence = (
            annotated.consequence.value
            if annotated.consequence is not None
            else None
        )
        assert record["functional_consequence"] == expected_consequence, (
            f"Record {i}: functional_consequence mismatch"
        )

        assert record["allele_frequency"] == annotated.allele_frequency, (
            f"Record {i}: allele_frequency mismatch"
        )
        if annotated.allele_frequency is not None:
            assert isinstance(record["allele_frequency"], float), (
                f"Record {i}: allele_frequency should be float"
            )
        else:
            assert record["allele_frequency"] is None

        assert record["composite_rank"] == scored.composite_rank, (
            f"Record {i}: composite_rank mismatch"
        )
        if scored.composite_rank is not None:
            assert isinstance(record["composite_rank"], float), (
                f"Record {i}: composite_rank should be float"
            )
        else:
            assert record["composite_rank"] is None

        expected_clinvar = (
            annotated.clinvar_assertion.value
            if annotated.clinvar_assertion is not None
            else None
        )
        assert record["clinvar_assertion"] == expected_clinvar, (
            f"Record {i}: clinvar_assertion mismatch"
        )

        expected_classification = (
            original.classification.value
            if original.classification is not None
            else None
        )
        assert record["acmg_classification"] == expected_classification, (
            f"Record {i}: acmg_classification mismatch"
        )

        expected_tags = sorted(tag.value for tag in original.evidence_tags)
        assert record["evidence_tags"] == expected_tags, (
            f"Record {i}: evidence_tags mismatch"
        )

@given(variants=classified_variant_list(min_size=0, max_size=15))
@settings(max_examples=100)
def test_json_output_field_order(
    variants: list[ClassifiedVariant],
) -> None:
    """JSON output fields appear in the specified order for every record.

    Order: chromosome, position, ref_allele, alt_allele, functional_consequence,
    allele_frequency, composite_rank, clinvar_assertion, acmg_classification,
    evidence_tags.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "report.json"
        write_json(variants, output_path)

        with open(output_path, "r", encoding="utf-8") as f:
            content = f.read()

    deserialized = json.loads(content, object_pairs_hook=list)

    expected_field_order = [
        "chromosome",
        "position",
        "ref_allele",
        "alt_allele",
        "functional_consequence",
        "allele_frequency",
        "composite_rank",
        "clinvar_assertion",
        "acmg_classification",
        "evidence_tags",
    ]

    for i, record_pairs in enumerate(deserialized):
        keys = [pair[0] for pair in record_pairs]
        assert keys == expected_field_order, (
            f"Record {i}: field order {keys} does not match expected {expected_field_order}"
        )

# ---------------------------------------------------------------------------

@given(variants=classified_variant_list(min_size=0, max_size=15))
@settings(max_examples=100)
def test_csv_has_correct_header_and_row_count(
    variants: list[ClassifiedVariant],
) -> None:
    """CSV output has one header row with correct field order and one data row per variant.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "report.csv"
        write_csv(variants, output_path)

        with open(output_path, "r", encoding="utf-8", newline="") as f:
            reader = csv.reader(f)
            rows = list(reader)

    assert len(rows) >= 1, "CSV must have at least a header row"

    header = rows[0]
    assert header == CSV_FIELDS, (
        f"CSV header {header} does not match expected {CSV_FIELDS}"
    )

    data_rows = rows[1:]
    assert len(data_rows) == len(variants), (
        f"Expected {len(variants)} data rows, got {len(data_rows)}"
    )

@given(variants=classified_variant_list(min_size=1, max_size=15))
@settings(max_examples=100)
def test_csv_rows_have_consistent_field_count(
    variants: list[ClassifiedVariant],
) -> None:
    """Every row in the CSV (header and data) has the same number of fields.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "report.csv"
        write_csv(variants, output_path)

        with open(output_path, "r", encoding="utf-8", newline="") as f:
            reader = csv.reader(f)
            rows = list(reader)

    expected_field_count = len(CSV_FIELDS)

    for i, row in enumerate(rows):
        assert len(row) == expected_field_count, (
            f"Row {i} has {len(row)} fields, expected {expected_field_count}"
        )
