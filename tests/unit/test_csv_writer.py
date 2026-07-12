"""Unit tests for the CSV report writer."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from vartriage.models.variant import (ACMGClassification, AnnotatedVariant,
                                      ClassifiedVariant, ClinVarAssertion,
                                      EvidenceTag, FunctionalConsequence,
                                      ScoredVariant, Variant)
from vartriage.reporting.csv_writer import CSV_FIELDS, write_csv


def _make_classified_variant(
    chrom: str = "chr1",
    pos: int = 12345,
    ref: str = "A",
    alt: str = "T",
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    allele_frequency: float | None = 0.002,
    composite_rank: float | None = 0.78,
    clinvar_assertion: ClinVarAssertion | None = ClinVarAssertion.PATHOGENIC,
    classification: ACMGClassification = ACMGClassification.LIKELY_PATHOGENIC,
    evidence_tags: frozenset[EvidenceTag] = frozenset(
        {EvidenceTag.PVS1, EvidenceTag.PM2}
    ),
) -> ClassifiedVariant:
    """Helper to build a ClassifiedVariant for testing."""
    base = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=base,
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
    )
    scored = ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
    )
    return ClassifiedVariant(
        scored=scored,
        evidence_tags=evidence_tags,
        classification=classification,
    )


class TestCSVWriter:
    """Tests for the write_csv function."""

    def test_writes_header_row(self, tmp_path: Path) -> None:
        """Header row matches the specified field order."""
        output = tmp_path / "output.csv"
        write_csv([], output)

        with open(output, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            header = next(reader)

        assert header == CSV_FIELDS

    def test_empty_variant_list_produces_header_only(self, tmp_path: Path) -> None:
        """An empty variant list produces a valid CSV with just the header."""
        output = tmp_path / "output.csv"
        write_csv([], output)

        with open(output, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert rows[0] == CSV_FIELDS

    def test_single_variant_row_values(self, tmp_path: Path) -> None:
        """A single fully-populated variant produces the correct field values."""
        variant = _make_classified_variant()
        output = tmp_path / "output.csv"
        write_csv([variant], output)

        with open(output, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        assert len(rows) == 2
        data_row = rows[1]
        assert data_row[0] == "chr1"
        assert data_row[1] == "12345"
        assert data_row[2] == "A"
        assert data_row[3] == "T"
        assert data_row[4] == "Missense"
        assert data_row[5] == "0.002"
        assert data_row[6] == "0.78"
        assert data_row[7] == "Pathogenic"
        assert data_row[8] == "Likely_Pathogenic"
        # Evidence tags sorted alphabetically by value
        assert data_row[9] == "PM2;PVS1"

    def test_absent_values_become_empty_fields(self, tmp_path: Path) -> None:
        """None values are represented as empty strings in CSV output."""
        variant = _make_classified_variant(
            allele_frequency=None,
            composite_rank=None,
            clinvar_assertion=None,
            evidence_tags=frozenset(),
        )
        output = tmp_path / "output.csv"
        write_csv([variant], output)

        with open(output, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        data_row = rows[1]
        assert data_row[5] == ""  # allele_frequency
        assert data_row[6] == ""  # composite_rank
        assert data_row[7] == ""  # clinvar_assertion
        assert data_row[9] == ""  # evidence_tags (empty frozenset)

    def test_multiple_variants_preserve_order(self, tmp_path: Path) -> None:
        """Multiple variants appear in the output in the same order as input."""
        variants = [
            _make_classified_variant(chrom="chr1", pos=100),
            _make_classified_variant(chrom="chr2", pos=200),
            _make_classified_variant(chrom="chr3", pos=300),
        ]
        output = tmp_path / "output.csv"
        write_csv(variants, output)

        with open(output, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        assert len(rows) == 4  # header + 3 data rows
        assert rows[1][0] == "chr1"
        assert rows[2][0] == "chr2"
        assert rows[3][0] == "chr3"

    def test_consistent_field_count_per_row(self, tmp_path: Path) -> None:
        """Every row has the same number of fields as the header."""
        variants = [
            _make_classified_variant(),
            _make_classified_variant(
                allele_frequency=None,
                composite_rank=None,
                clinvar_assertion=None,
                evidence_tags=frozenset(),
            ),
        ]
        output = tmp_path / "output.csv"
        write_csv(variants, output)

        with open(output, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        expected_count = len(CSV_FIELDS)
        for row in rows:
            assert len(row) == expected_count

    def test_utf8_encoding(self, tmp_path: Path) -> None:
        """Output file is UTF-8 encoded."""
        output = tmp_path / "output.csv"
        write_csv([], output)

        raw_bytes = output.read_bytes()
        # Should be decodable as UTF-8 without errors
        raw_bytes.decode("utf-8")

    def test_returns_output_path(self, tmp_path: Path) -> None:
        """write_csv returns the output path."""
        output = tmp_path / "output.csv"
        result = write_csv([], output)
        assert result == output

    def test_rfc4180_comma_delimiter(self, tmp_path: Path) -> None:
        """Fields are comma-delimited per RFC 4180."""
        variant = _make_classified_variant()
        output = tmp_path / "output.csv"
        write_csv([variant], output)

        lines = output.read_text(encoding="utf-8").splitlines()
        # Header line uses commas
        assert "," in lines[0]
        # Each header field appears as expected
        assert lines[0].startswith("chromosome,position,")
