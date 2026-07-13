"""Hypothesis tests for quality filtering and configuration validation."""

from __future__ import annotations

import warnings

from hypothesis import given, settings
from hypothesis import strategies as st

from tests.generators.variants import (FILTER_FAIL_VALUES, FILTER_PASS_VALUES,
                                       chromosome, genomic_position,
                                       snv_allele)
from vartriage.filter.quality_filter import QualityFilter
from vartriage.models.config import (AnnotationConfig, PrioritizationConfig,
                                     QualityFilterConfig)
from vartriage.models.variant import Variant
from vartriage.models.warnings import MissingDataWarning

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

VALID_QUAL_THRESHOLD = st.floats(
    min_value=0.0, max_value=1_000_000.0, allow_nan=False, allow_infinity=False
)


@st.composite
def arbitrary_variant_for_filtering(draw: st.DrawFn) -> Variant:
    """Generate a Variant with random FILTER and QUAL values.

    Covers passing filters, failing filters, present QUAL, and missing QUAL.
    """
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())

    filter_status = draw(st.sampled_from(FILTER_PASS_VALUES + FILTER_FAIL_VALUES))

    qual: float | None = draw(
        st.one_of(
            st.none(),
            st.floats(
                min_value=0.0,
                max_value=10000.0,
                allow_nan=False,
                allow_infinity=False,
            ),
        )
    )

    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=qual,
        filter_status=filter_status,
        info={},
    )


# ---------------------------------------------------------------------------


@given(
    variant=arbitrary_variant_for_filtering(),
    threshold=VALID_QUAL_THRESHOLD,
)
@settings(max_examples=200)
def test_quality_filter_correctness(variant: Variant, threshold: float) -> None:
    """Variant passes iff FILTER in {PASS, .} AND QUAL is not None AND QUAL >= threshold."""
    config = QualityFilterConfig(min_qual=threshold)
    qf = QualityFilter(config)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = list(qf.apply(iter([variant])))

    should_pass = (
        variant.filter_status in {"PASS", "."}
        and variant.qual is not None
        and variant.qual >= threshold
    )

    if should_pass:
        assert result == [variant], (
            f"Expected variant to pass but it was excluded. "
            f"filter={variant.filter_status}, qual={variant.qual}, threshold={threshold}"
        )
    else:
        assert result == [], (
            f"Expected variant to be excluded but it passed. "
            f"filter={variant.filter_status}, qual={variant.qual}, threshold={threshold}"
        )

    # Variants excluded due to missing QUAL (with passing FILTER) trigger a warning
    if variant.filter_status in {"PASS", "."} and variant.qual is None:
        missing_qual_warnings = [
            w
            for w in caught
            if issubclass(w.category, UserWarning)
            and w.message.args
            and isinstance(w.message.args[0], MissingDataWarning)
        ]
        assert (
            len(missing_qual_warnings) == 1
        ), "Expected exactly one MissingDataWarning for missing QUAL"
        warning_obj: MissingDataWarning = missing_qual_warnings[0].message.args[0]  # type: ignore[union-attr]
        assert variant.chrom in warning_obj.reason  # type: ignore[operator]
        assert str(variant.pos) in warning_obj.reason  # type: ignore[operator]


# ---------------------------------------------------------------------------


@given(
    variants=st.lists(arbitrary_variant_for_filtering(), min_size=0, max_size=50),
    threshold=VALID_QUAL_THRESHOLD,
)
@settings(max_examples=100)
def test_filtering_preserves_ordering(
    variants: list[Variant], threshold: float
) -> None:
    """Relative ordering of passing variants matches input ordering."""
    config = QualityFilterConfig(min_qual=threshold)
    qf = QualityFilter(config)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = list(qf.apply(iter(variants)))

    # Compute which variants should pass using the same predicate
    expected = [
        v
        for v in variants
        if v.filter_status in {"PASS", "."}
        and v.qual is not None
        and v.qual >= threshold
    ]

    assert (
        result == expected
    ), "Filtered output does not preserve the relative input ordering"


# ---------------------------------------------------------------------------


@given(
    value=st.one_of(
        st.floats(
            max_value=-0.01,
            allow_nan=False,
            allow_infinity=False,
        ),
        st.floats(
            min_value=1_000_000.01,
            max_value=1e12,
            allow_nan=False,
            allow_infinity=False,
        ),
    )
)
@settings(max_examples=100)
def test_quality_filter_config_rejects_out_of_range(value: float) -> None:
    """QualityFilterConfig rejects min_qual outside [0, 1_000_000]."""
    try:
        QualityFilterConfig(min_qual=value)
        assert False, f"Expected ValueError for min_qual={value}"
    except ValueError as exc:
        assert "0" in str(exc) and "1000000" in str(
            exc
        ), f"Error message should specify valid range, got: {exc}"


@given(
    value=st.one_of(
        st.floats(
            max_value=-0.001,
            allow_nan=False,
            allow_infinity=False,
        ),
        st.floats(
            min_value=1.001,
            max_value=100.0,
            allow_nan=False,
            allow_infinity=False,
        ),
    )
)
@settings(max_examples=100)
def test_prioritization_config_rejects_out_of_range_af(value: float) -> None:
    """PrioritizationConfig rejects max_allele_frequency outside [0.0, 1.0]."""
    try:
        PrioritizationConfig(max_allele_frequency=value)
        assert False, f"Expected ValueError for max_allele_frequency={value}"
    except ValueError as exc:
        assert "0.0" in str(exc) and "1.0" in str(
            exc
        ), f"Error message should specify valid range, got: {exc}"


@given(
    value=st.one_of(
        st.integers(min_value=-10000, max_value=999),
        st.integers(min_value=100_001, max_value=1_000_000),
    )
)
@settings(max_examples=100)
def test_annotation_config_rejects_out_of_range_batch_size(value: int) -> None:
    """AnnotationConfig rejects batch_size outside [1_000, 100_000]."""
    from pathlib import Path

    try:
        AnnotationConfig(
            gene_annotation_path=Path("/fake/genes.gtf"),
            gnomad_path=Path("/fake/gnomad.tsv"),
            batch_size=value,
        )
        assert False, f"Expected ValueError for batch_size={value}"
    except ValueError as exc:
        assert "1000" in str(exc) and "100000" in str(
            exc
        ), f"Error message should specify valid range, got: {exc}"


@given(
    value=st.one_of(
        st.integers(min_value=-10000, max_value=999),
        st.integers(min_value=100_001, max_value=1_000_000),
    )
)
@settings(max_examples=100)
def test_prioritization_config_rejects_out_of_range_batch_size(
    value: int,
) -> None:
    """PrioritizationConfig rejects batch_size outside [1_000, 100_000]."""
    try:
        PrioritizationConfig(batch_size=value)
        assert False, f"Expected ValueError for batch_size={value}"
    except ValueError as exc:
        assert "1000" in str(exc) and "100000" in str(
            exc
        ), f"Error message should specify valid range, got: {exc}"
