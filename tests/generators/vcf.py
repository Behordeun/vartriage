"""Hypothesis strategies for generating VCF file content.

Produces VCF 4.2 format content with proper headers and data lines,
as well as intentionally malformed content for error-handling tests.
"""

from __future__ import annotations

from hypothesis import strategies as st
from hypothesis.strategies import SearchStrategy

from tests.generators.variants import (
    CHROMOSOMES,
    FILTER_FAIL_VALUES,
    FILTER_PASS_VALUES,
    NUCLEOTIDES,
)


VCF_FILEFORMAT_LINE = "##fileformat=VCFv4.2"

VCF_INFO_LINES = [
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">',
]

VCF_FORMAT_LINES = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
]

VCF_COLUMN_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def _build_vcf_header(
    include_fileformat: bool = True,
    include_info: bool = True,
    include_format: bool = True,
    include_column_header: bool = True,
) -> str:
    """Construct a VCF header from components."""
    lines: list[str] = []
    if include_fileformat:
        lines.append(VCF_FILEFORMAT_LINE)
    if include_info:
        lines.extend(VCF_INFO_LINES)
    if include_format:
        lines.extend(VCF_FORMAT_LINES)
    if include_column_header:
        lines.append(VCF_COLUMN_HEADER)
    return "\n".join(lines)


@st.composite
def vcf_data_line(draw: st.DrawFn) -> str:
    """Generate a single valid VCF data line.

    Produces a tab-separated line with the 8 mandatory VCF columns:
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
    """
    chrom = draw(st.sampled_from(CHROMOSOMES))
    pos = draw(st.integers(min_value=1, max_value=250_000_000))

    variant_id = draw(
        st.one_of(
            st.just("."),
            st.from_regex(r"rs[0-9]{1,12}", fullmatch=True),
        )
    )

    ref_len = draw(st.integers(min_value=1, max_value=5))
    ref = "".join(draw(st.sampled_from(NUCLEOTIDES)) for _ in range(ref_len))

    alt_len = draw(st.integers(min_value=1, max_value=5))
    alt = "".join(draw(st.sampled_from(NUCLEOTIDES)) for _ in range(alt_len))

    qual = draw(
        st.one_of(
            st.just("."),
            st.floats(min_value=0.0, max_value=10000.0, allow_nan=False).map(
                lambda x: f"{x:.2f}"
            ),
        )
    )

    filter_status = draw(
        st.sampled_from(FILTER_PASS_VALUES + FILTER_FAIL_VALUES + ["."])
    )

    info = draw(
        st.one_of(
            st.just("."),
            st.just(f"DP={draw(st.integers(min_value=1, max_value=1000))}"),
            st.just(
                f"DP={draw(st.integers(min_value=1, max_value=1000))};"
                f"MQ={draw(st.floats(min_value=0, max_value=60, allow_nan=False)):.1f}"
            ),
        )
    )

    fields = [chrom, str(pos), variant_id, ref, alt, qual, filter_status, info]
    return "\t".join(fields)


@st.composite
def valid_vcf_content(
    draw: st.DrawFn,
    min_variants: int = 1,
    max_variants: int = 20,
) -> str:
    """Generate complete, valid VCF 4.2 file content.

    Produces a well-formed VCF with:
    - ``##fileformat=VCFv4.2`` declaration
    - One or more ``##INFO`` meta-information lines
    - The mandatory ``#CHROM`` column header line
    - Between ``min_variants`` and ``max_variants`` data lines

    Parameters
    ----------
    min_variants : int
        Minimum number of variant data lines.
    max_variants : int
        Maximum number of variant data lines.
    """
    header = _build_vcf_header()
    num_variants = draw(st.integers(min_value=min_variants, max_value=max_variants))
    data_lines = [draw(vcf_data_line()) for _ in range(num_variants)]
    return header + "\n" + "\n".join(data_lines) + "\n"


@st.composite
def malformed_vcf_content(draw: st.DrawFn) -> tuple[str, str]:
    """Generate VCF content with a specific header or data violation.

    Returns a tuple of (malformed_content, violation_type) where
    violation_type describes the nature of the malformation for
    assertion purposes.

    Violation types:
    - "missing_fileformat": Missing ##fileformat declaration
    - "missing_column_header": Missing #CHROM column header line
    - "malformed_info": Invalid ##INFO meta-information format
    - "missing_columns": Data line with fewer than 8 tab-separated fields
    - "invalid_pos": Data line with a non-integer POS field
    - "invalid_qual": Data line with non-numeric, non-dot QUAL field
    """
    violation = draw(
        st.sampled_from(
            [
                "missing_fileformat",
                "missing_column_header",
                "malformed_info",
                "missing_columns",
                "invalid_pos",
                "invalid_qual",
            ]
        )
    )

    if violation == "missing_fileformat":
        lines = [
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
            VCF_COLUMN_HEADER,
            "chr1\t100\t.\tA\tT\t30.0\tPASS\t.",
        ]
        content = "\n".join(lines) + "\n"

    elif violation == "missing_column_header":
        lines = [
            VCF_FILEFORMAT_LINE,
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
            "chr1\t100\t.\tA\tT\t30.0\tPASS\t.",
        ]
        content = "\n".join(lines) + "\n"

    elif violation == "malformed_info":
        lines = [
            VCF_FILEFORMAT_LINE,
            "##INFO=<BROKEN LINE WITHOUT PROPER FORMAT>",
            VCF_COLUMN_HEADER,
            "chr1\t100\t.\tA\tT\t30.0\tPASS\t.",
        ]
        content = "\n".join(lines) + "\n"

    elif violation == "missing_columns":
        lines = [
            VCF_FILEFORMAT_LINE,
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
            VCF_COLUMN_HEADER,
            "chr1\t100\t.\tA",
        ]
        content = "\n".join(lines) + "\n"

    elif violation == "invalid_pos":
        lines = [
            VCF_FILEFORMAT_LINE,
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
            VCF_COLUMN_HEADER,
            "chr1\tNOT_A_NUMBER\t.\tA\tT\t30.0\tPASS\t.",
        ]
        content = "\n".join(lines) + "\n"

    elif violation == "invalid_qual":
        lines = [
            VCF_FILEFORMAT_LINE,
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
            VCF_COLUMN_HEADER,
            "chr1\t100\t.\tA\tT\tINVALID_QUAL\tPASS\t.",
        ]
        content = "\n".join(lines) + "\n"

    else:
        content = ""

    return content, violation
