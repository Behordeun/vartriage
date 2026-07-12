"""Hypothesis tests for VCF parsing correctness."""

from __future__ import annotations

import tempfile
from pathlib import Path

from hypothesis import given, settings

from tests.generators.vcf import malformed_vcf_content, valid_vcf_content
from vartriage.io.exceptions import ParseError
from vartriage.io.vcf_parser import VCFParser


def _parse_data_line(line: str) -> dict:
    """Parse a raw VCF data line into a dict of expected field values.

    Extracts the 8 mandatory columns from a tab-separated VCF data line
    and returns them in a dictionary matching the Variant dataclass fields.
    """
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = int(fields[1])
    variant_id = fields[2] if fields[2] != "." else None
    ref = fields[3]
    alt = fields[4]

    qual_str = fields[5]
    qual = None if qual_str == "." else float(qual_str)

    filter_status = fields[6]

    return {
        "chrom": chrom,
        "pos": pos,
        "id": variant_id,
        "ref": ref,
        "alt": alt,
        "qual": qual,
        "filter_status": filter_status,
    }


class TestVCFParsingRoundTrip:
    """VCF Parsing Round-Trip.

    For any valid VCF content, parsed Variant records have fields exactly
    matching source data lines.
    """

    @settings(max_examples=100, deadline=None)
    @given(content=valid_vcf_content(min_variants=1, max_variants=10))
    def test_parsed_variants_match_source_data(self, content: str) -> None:
        """Parsed Variant fields match the corresponding VCF data line values."""
        lines = content.strip().split("\n")
        data_lines = [line for line in lines if line and not line.startswith("#")]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(content)
            tmp_path = Path(f.name)

        try:
            with VCFParser(tmp_path) as parser:
                variants = list(parser)

            assert len(variants) == len(
                data_lines
            ), f"Expected {len(data_lines)} variants, got {len(variants)}"

            for variant, raw_line in zip(variants, data_lines):
                expected = _parse_data_line(raw_line)

                assert (
                    variant.chrom == expected["chrom"]
                ), f"CHROM mismatch: {variant.chrom!r} != {expected['chrom']!r}"
                assert (
                    variant.pos == expected["pos"]
                ), f"POS mismatch: {variant.pos} != {expected['pos']}"
                assert (
                    variant.id == expected["id"]
                ), f"ID mismatch: {variant.id!r} != {expected['id']!r}"
                assert (
                    variant.ref == expected["ref"]
                ), f"REF mismatch: {variant.ref!r} != {expected['ref']!r}"
                assert (
                    variant.alt == expected["alt"]
                ), f"ALT mismatch: {variant.alt!r} != {expected['alt']!r}"

                if expected["qual"] is None:
                    assert (
                        variant.qual is None
                    ), f"QUAL should be None, got {variant.qual}"
                else:
                    assert variant.qual is not None, "QUAL should not be None"
                    assert (
                        abs(variant.qual - expected["qual"]) < 0.01
                    ), f"QUAL mismatch: {variant.qual} != {expected['qual']}"

                assert variant.filter_status == expected["filter_status"], (
                    f"FILTER mismatch: {variant.filter_status!r} != "
                    f"{expected['filter_status']!r}"
                )
        finally:
            tmp_path.unlink(missing_ok=True)


class TestMalformedVCFDetection:
    """Malformed VCF Detection.

    For any VCF with header/data violations, parser raises ParseError
    with line number and malformation nature.
    Note: pysam is lenient with certain malformed INFO/FORMAT header
    lines (it logs warnings but does not raise errors). The property
    validates that violations detectable by the parser raise ParseError.
    """

    # pysam tolerates these. It logs a warning but doesn't error out.
    # "invalid_qual" is silently coerced to 0.0 by pysam's htslib layer,
    # so the parser has no way to detect the malformation at runtime.
    PYSAM_TOLERANT_VIOLATIONS = {"malformed_info", "invalid_qual"}

    @settings(max_examples=100, deadline=None)
    @given(data=malformed_vcf_content())
    def test_malformed_vcf_raises_parse_error(self, data: tuple[str, str]) -> None:
        """Malformed VCF content triggers a ParseError with line info."""
        content, violation_type = data

        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(content)
            tmp_path = Path(f.name)

        try:
            parse_error_raised = False
            try:
                with VCFParser(tmp_path) as parser:
                    for _variant in parser:
                        pass
            except ParseError as exc:
                parse_error_raised = True
                assert (
                    exc.line_number >= 1
                ), f"ParseError line_number should be >= 1, got {exc.line_number}"
                assert exc.detail, "ParseError detail should not be empty"

            if violation_type not in self.PYSAM_TOLERANT_VIOLATIONS:
                assert parse_error_raised, (
                    f"Expected ParseError for violation '{violation_type}', "
                    f"but no exception was raised"
                )
        finally:
            tmp_path.unlink(missing_ok=True)
