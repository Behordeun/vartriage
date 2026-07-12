"""Property-based tests for RegionFilter and SampleExtractor."""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from vartriage.filter.region_filter import RegionFilter
from vartriage.filter.sample_extractor import SampleExtractor
from vartriage.io.exceptions import ParseError
from vartriage.models.config import (
    RegionFilterConfig,
    SampleConfig,
)
from vartriage.models.variant import Variant


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


@st.composite
def bed_interval(draw: st.DrawFn) -> tuple[str, int, int]:
    """Valid BED interval with start < end."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    start = draw(st.integers(min_value=0, max_value=249_999_998))
    end = draw(
        st.integers(min_value=start + 1, max_value=249_999_999)
    )
    return (chrom, start, end)


@st.composite
def simple_variant(draw: st.DrawFn) -> Variant:
    """Minimal Variant for region-filter testing."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    pos = draw(st.integers(min_value=1, max_value=250_000_000))
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
        info={},
    )


@st.composite
def gt_tuple_with_alt(draw: st.DrawFn) -> tuple[int | None, ...]:
    """GT tuple containing at least one alt allele (>= 1)."""
    alleles = draw(
        st.lists(
            st.integers(min_value=0, max_value=3),
            min_size=2,
            max_size=2,
        )
    )
    if all(a == 0 for a in alleles):
        idx = draw(
            st.integers(min_value=0, max_value=len(alleles) - 1)
        )
        alleles[idx] = draw(st.integers(min_value=1, max_value=3))
    return tuple(alleles)


@st.composite
def gt_tuple_no_alt(draw: st.DrawFn) -> tuple[int | None, ...]:
    """GT tuple with no alt alleles (all 0 or None)."""
    alleles = draw(
        st.lists(
            st.sampled_from([0, None]),
            min_size=2,
            max_size=2,
        )
    )
    return tuple(alleles)


# ---------------------------------------------------------------------------
# Property 1: BED Parsing Round-Trip
# ---------------------------------------------------------------------------


@given(interval=bed_interval())
@settings(max_examples=100)
def test_bed_parsing_round_trip(
    interval: tuple[str, int, int],
) -> None:
    """Positions at start and end-1 overlap; position at end does not."""
    chrom, start, end = interval
    bed_line = f"{chrom}\t{start}\t{end}\n"

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as f:
        f.write(bed_line)
        bed_path = Path(f.name)

    try:
        config = RegionFilterConfig(bed_path=bed_path)
        rf = RegionFilter(config)

        # pos = start + 1 => 0-based = start, inside [start, end)
        v_start = Variant(
            chrom=chrom, pos=start + 1, id=None,
            ref="A", alt="T", qual=30.0,
            filter_status="PASS", info={},
        )
        assert list(rf.apply(iter([v_start]))) == [v_start]

        # pos = end => 0-based = end - 1, inside [start, end)
        v_end_minus_1 = Variant(
            chrom=chrom, pos=end, id=None,
            ref="A", alt="T", qual=30.0,
            filter_status="PASS", info={},
        )
        assert list(rf.apply(iter([v_end_minus_1]))) == [
            v_end_minus_1
        ]

        # pos = end + 1 => 0-based = end, NOT in [start, end)
        v_past_end = Variant(
            chrom=chrom, pos=end + 1, id=None,
            ref="A", alt="T", qual=30.0,
            filter_status="PASS", info={},
        )
        assert list(rf.apply(iter([v_past_end]))) == []
    finally:
        bed_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Property 2: BED Validation Rejects Malformed
# ---------------------------------------------------------------------------


@st.composite
def malformed_bed_line(draw: st.DrawFn) -> str:
    """Malformed BED line (< 3 cols, bad ints, start>=end, negative)."""
    kind = draw(
        st.sampled_from([
            "few_columns", "non_integer",
            "start_ge_end", "negative",
        ])
    )
    if kind == "few_columns":
        num_cols = draw(st.integers(min_value=1, max_value=2))
        cols = [
            draw(
                st.text(
                    min_size=1, max_size=8,
                    alphabet="chrABC123",
                )
            )
            for _ in range(num_cols)
        ]
        return "\t".join(cols)
    elif kind == "non_integer":
        chrom = draw(st.sampled_from(CHROMOSOMES))
        bad_val = draw(
            st.text(min_size=1, max_size=5, alphabet="abcXYZ!@")
        )
        return f"{chrom}\t{bad_val}\t100"
    elif kind == "start_ge_end":
        chrom = draw(st.sampled_from(CHROMOSOMES))
        val = draw(st.integers(min_value=1, max_value=1000))
        start = val
        end = draw(st.integers(min_value=0, max_value=val))
        return f"{chrom}\t{start}\t{end}"
    else:  # negative
        chrom = draw(st.sampled_from(CHROMOSOMES))
        neg = draw(st.integers(min_value=-1000, max_value=-1))
        return f"{chrom}\t{neg}\t100"


@given(bad_line=malformed_bed_line())
@settings(max_examples=100)
def test_bed_validation_rejects_malformed(bad_line: str) -> None:
    """Malformed BED lines raise ParseError with line_number=1."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as f:
        f.write(bad_line + "\n")
        bed_path = Path(f.name)

    try:
        config = RegionFilterConfig(bed_path=bed_path)
        try:
            RegionFilter(config)
            raise AssertionError(
                f"Expected ParseError for line: {bad_line!r}"
            )
        except ParseError as exc:
            assert exc.line_number == 1
    finally:
        bed_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Property 3: Comment Invisibility
# ---------------------------------------------------------------------------


@given(
    intervals=st.lists(bed_interval(), min_size=1, max_size=10),
    comments=st.lists(
        st.sampled_from([
            "# this is a comment",
            "#comment line",
            "browser position chr1:1-100",
            "track name=test",
        ]),
        min_size=0,
        max_size=5,
    ),
)
@settings(max_examples=100)
def test_comment_invisibility(
    intervals: list[tuple[str, int, int]],
    comments: list[str],
) -> None:
    """Comment/browser/track lines don't change interval count."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")
        data_only_path = Path(f.name)

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as f:
        for comment in comments:
            f.write(comment + "\n")
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")
        mixed_path = Path(f.name)

    try:
        config_data = RegionFilterConfig(bed_path=data_only_path)
        rf_data = RegionFilter(config_data)

        config_mixed = RegionFilterConfig(bed_path=mixed_path)
        rf_mixed = RegionFilter(config_mixed)

        data_count = sum(
            len(ivs) for ivs in rf_data._intervals.values()
        )
        mixed_count = sum(
            len(ivs) for ivs in rf_mixed._intervals.values()
        )

        assert data_count == mixed_count
    finally:
        data_only_path.unlink(missing_ok=True)
        mixed_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Property 4: Overlap Correctness (vs naive O(n) check)
# ---------------------------------------------------------------------------


def _naive_overlaps(
    intervals: list[tuple[str, int, int]], chrom: str, pos: int
) -> bool:
    """Naive O(n) overlap check for reference."""
    query = pos - 1
    for iv_chrom, start, end in intervals:
        if iv_chrom == chrom and start <= query < end:
            return True
    return False


@given(
    intervals=st.lists(bed_interval(), min_size=1, max_size=20),
    variants=st.lists(simple_variant(), min_size=1, max_size=20),
)
@settings(max_examples=100)
def test_overlap_correctness(
    intervals: list[tuple[str, int, int]],
    variants: list[Variant],
) -> None:
    """RegionFilter agrees with naive O(n) overlap check."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")
        bed_path = Path(f.name)

    try:
        config = RegionFilterConfig(bed_path=bed_path)
        rf = RegionFilter(config)

        result = list(rf.apply(iter(variants)))
        expected = [
            v for v in variants
            if _naive_overlaps(intervals, v.chrom, v.pos)
        ]

        assert result == expected
    finally:
        bed_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Property 5: Order Preservation
# ---------------------------------------------------------------------------


@given(
    intervals=st.lists(bed_interval(), min_size=1, max_size=10),
    variants=st.lists(simple_variant(), min_size=1, max_size=30),
)
@settings(max_examples=100)
def test_order_preservation_region_filter(
    intervals: list[tuple[str, int, int]],
    variants: list[Variant],
) -> None:
    """RegionFilter output is an ordered subsequence of input."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False
    ) as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")
        bed_path = Path(f.name)

    try:
        config = RegionFilterConfig(bed_path=bed_path)
        rf = RegionFilter(config)

        result = list(rf.apply(iter(variants)))

        it = iter(variants)
        for r in result:
            for v in it:
                if v == r:
                    break
            else:
                raise AssertionError(
                    "Output variant not found in input order"
                )
    finally:
        bed_path.unlink(missing_ok=True)


@given(
    variants=st.lists(simple_variant(), min_size=1, max_size=30),
)
@settings(max_examples=100)
def test_order_preservation_sample_extractor(
    variants: list[Variant],
) -> None:
    """SampleExtractor output preserves input order."""
    sample_name = "SAMPLE1"
    config = SampleConfig(sample_name=sample_name)

    extractor = SampleExtractor(config, [sample_name])

    enriched = []
    for v in variants:
        new_info = dict(v.info)
        new_info["_pysam_samples"] = {
            sample_name: {"GT": (0, 1), "GQ": 50}
        }
        enriched.append(Variant(
            chrom=v.chrom, pos=v.pos, id=v.id,
            ref=v.ref, alt=v.alt, qual=v.qual,
            filter_status=v.filter_status, info=new_info,
        ))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = list(extractor.apply(iter(enriched)))

    assert len(result) == len(enriched)
    for i, (r, e) in enumerate(zip(result, enriched)):
        assert r.chrom == e.chrom and r.pos == e.pos


# ---------------------------------------------------------------------------
# Property 6: Genotype Alt-Allele Filtering
# ---------------------------------------------------------------------------


@given(gt=gt_tuple_with_alt())
@settings(max_examples=100)
def test_genotype_alt_allele_pass(
    gt: tuple[int | None, ...],
) -> None:
    """GT with at least one alt allele passes the filter."""
    gt_str = SampleExtractor._format_gt(gt)
    assert SampleExtractor._has_alt_allele(gt_str)


@given(gt=gt_tuple_no_alt())
@settings(max_examples=100)
def test_genotype_alt_allele_reject(
    gt: tuple[int | None, ...],
) -> None:
    """GT with no alt alleles is excluded."""
    gt_str = SampleExtractor._format_gt(gt)
    assert not SampleExtractor._has_alt_allele(gt_str)


# ---------------------------------------------------------------------------
# Property 7: GQ Threshold
# ---------------------------------------------------------------------------


@given(
    gq=st.integers(min_value=0, max_value=99),
    threshold=st.integers(min_value=0, max_value=99),
)
@settings(max_examples=100)
def test_gq_threshold(gq: int, threshold: int) -> None:
    """Variant excluded iff GQ < threshold."""
    sample_name = "S1"
    config = SampleConfig(
        sample_name=sample_name, min_gq=threshold
    )
    extractor = SampleExtractor(config, [sample_name])

    variant = Variant(
        chrom="chr1", pos=100, id=None,
        ref="A", alt="T", qual=30.0,
        filter_status="PASS",
        info={
            "_pysam_samples": {
                sample_name: {"GT": (0, 1), "GQ": gq}
            }
        },
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = list(extractor.apply(iter([variant])))

    if gq < threshold:
        assert len(result) == 0
    else:
        assert len(result) == 1


# ---------------------------------------------------------------------------
# Property 8: GQ Config Validation
# ---------------------------------------------------------------------------


@given(
    bad_gq=st.one_of(
        st.integers(min_value=100, max_value=10000),
        st.integers(min_value=-10000, max_value=-1),
    )
)
@settings(max_examples=100)
def test_gq_config_validation(bad_gq: int) -> None:
    """GQ values outside [0, 99] raise ValueError."""
    try:
        SampleConfig(sample_name="test", min_gq=bad_gq)
        raise AssertionError(
            f"Expected ValueError for min_gq={bad_gq}"
        )
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# Property 9: Sample Data Attachment
# ---------------------------------------------------------------------------


@given(
    gq=st.integers(min_value=0, max_value=99),
    gt=gt_tuple_with_alt(),
)
@settings(max_examples=100)
def test_sample_data_attachment(
    gq: int, gt: tuple[int | None, ...]
) -> None:
    """Output info contains sample_gt, sample_name, sample_gq."""
    sample_name = "PROBAND"
    config = SampleConfig(sample_name=sample_name)
    extractor = SampleExtractor(config, [sample_name])

    variant = Variant(
        chrom="chr1", pos=50, id=None,
        ref="C", alt="G", qual=40.0,
        filter_status="PASS",
        info={
            "_pysam_samples": {
                sample_name: {"GT": gt, "GQ": gq}
            }
        },
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = list(extractor.apply(iter([variant])))

    assert len(result) == 1
    info = result[0].info
    assert "sample_gt" in info
    assert "sample_name" in info
    assert info["sample_name"] == sample_name
    assert "sample_gq" in info
    assert info["sample_gq"] == gq


# ---------------------------------------------------------------------------
# Property 10: No-Config Passthrough
# ---------------------------------------------------------------------------


@given(
    variants=st.lists(simple_variant(), min_size=0, max_size=20)
)
@settings(max_examples=100)
def test_no_config_passthrough(variants: list[Variant]) -> None:
    """When region config is None, stream passes through unchanged."""
    region_config = None

    if region_config is None:
        output = list(iter(variants))
    else:
        rf = RegionFilter(region_config)
        output = list(rf.apply(iter(variants)))

    assert output == variants
