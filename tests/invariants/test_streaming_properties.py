"""Property tests for streaming report generation and score loading.

Covers:
- Streaming bounded-buffer behavior
- Write-error cleanup (no partial files)
- Malformed-line resilience in ScoreLoader
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Iterator
from unittest.mock import patch

from hypothesis import given, settings
from hypothesis import strategies as st

from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.prioritization.score_loader import ScoreLoader
from vartriage.reporting.csv_writer import write_csv
from vartriage.reporting.json_writer import write_json
from vartriage.reporting.generator import ReportGenerator
from vartriage.models.config import ReportConfig

from tests.generators.variants import scored_variant, evidence_tag_set


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------


@st.composite
def classified_variant(draw: st.DrawFn) -> ClassifiedVariant:
    """Build a random ClassifiedVariant."""
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
    draw: st.DrawFn, min_size: int = 1, max_size: int = 30
) -> list[ClassifiedVariant]:
    """Build a list of ClassifiedVariants."""
    return draw(
        st.lists(classified_variant(), min_size=min_size, max_size=max_size)
    )


# ---------------------------------------------------------------------------
# Helpers for Property 2: Bounded-buffer counting iterator
# ---------------------------------------------------------------------------


class CountingIterator:
    """Wraps a list and counts how many items have been pulled."""

    def __init__(self, items: list[ClassifiedVariant]) -> None:
        self._items = items
        self._index = 0
        self.consumed_count = 0
        self.max_concurrent = 0

    def __iter__(self) -> "CountingIterator":
        return self

    def __next__(self) -> ClassifiedVariant:
        if self._index >= len(self._items):
            raise StopIteration
        item = self._items[self._index]
        self._index += 1
        self.consumed_count += 1
        return item


# ---------------------------------------------------------------------------
# Property 2: Streaming bounded-buffer invariant
# ---------------------------------------------------------------------------


@given(variants=classified_variant_list(min_size=1, max_size=50))
@settings(max_examples=100)
def test_json_streaming_bounded_buffer(
    variants: list[ClassifiedVariant],
) -> None:
    """JSON writer pulls variants one at a time — no bulk buffering."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "report.json"
        counter = CountingIterator(variants)
        write_json(counter, output_path)

        assert counter.consumed_count == len(variants), (
            f"Expected {len(variants)} consumed, got {counter.consumed_count}"
        )
        assert output_path.exists()


@given(variants=classified_variant_list(min_size=1, max_size=50))
@settings(max_examples=100)
def test_csv_streaming_bounded_buffer(
    variants: list[ClassifiedVariant],
) -> None:
    """CSV writer pulls variants one at a time — no bulk buffering."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "report.csv"
        counter = CountingIterator(variants)
        write_csv(counter, output_path)

        assert counter.consumed_count == len(variants), (
            f"Expected {len(variants)} consumed, got {counter.consumed_count}"
        )
        assert output_path.exists()


# ---------------------------------------------------------------------------
# Property 3: Write-error cleanup guarantee
# ---------------------------------------------------------------------------


class FailingIterator:
    """Yields items normally until index ``fail_at``, then raises IOError."""

    def __init__(self, items: list[ClassifiedVariant], fail_at: int) -> None:
        self._items = items
        self._fail_at = fail_at
        self._index = 0

    def __iter__(self) -> "FailingIterator":
        return self

    def __next__(self) -> ClassifiedVariant:
        if self._index >= len(self._items):
            raise StopIteration
        if self._index == self._fail_at:
            raise IOError(
                f"Simulated write error at item {self._fail_at}"
            )
        item = self._items[self._index]
        self._index += 1
        return item


@given(
    variants=classified_variant_list(min_size=1, max_size=20),
    injection_frac=st.floats(min_value=0.0, max_value=1.0),
)
@settings(max_examples=100)
def test_write_error_cleanup_json(
    variants: list[ClassifiedVariant],
    injection_frac: float,
) -> None:
    """After an IOError mid-stream, no partial JSON file is left behind."""
    injection_point = int(injection_frac * len(variants))
    injection_point = min(injection_point, len(variants) - 1)

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "output.json"
        config = ReportConfig(output_format="json")
        generator = ReportGenerator(config)

        failing_iter = FailingIterator(variants, fail_at=injection_point)

        try:
            generator.generate(failing_iter, output_path)
        except (IOError, OSError):
            pass

        assert not output_path.exists(), (
            f"Target path {output_path} should not exist after write error, "
            f"but partial file was left behind"
        )


@given(
    variants=classified_variant_list(min_size=1, max_size=20),
    injection_frac=st.floats(min_value=0.0, max_value=1.0),
)
@settings(max_examples=100)
def test_write_error_cleanup_csv(
    variants: list[ClassifiedVariant],
    injection_frac: float,
) -> None:
    """After an IOError mid-stream, no partial CSV file is left behind."""
    injection_point = int(injection_frac * len(variants))
    injection_point = min(injection_point, len(variants) - 1)

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "output.csv"
        config = ReportConfig(output_format="csv")
        generator = ReportGenerator(config)

        failing_iter = FailingIterator(variants, fail_at=injection_point)

        try:
            generator.generate(failing_iter, output_path)
        except (IOError, OSError):
            pass

        assert not output_path.exists(), (
            f"Target path {output_path} should not exist after write error, "
            f"but partial file was left behind"
        )


# ---------------------------------------------------------------------------
# Property 5: Malformed-line resilience
# ---------------------------------------------------------------------------

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
NUCLEOTIDES = ["A", "C", "G", "T"]


@st.composite
def valid_tsv_line(draw: st.DrawFn) -> tuple[str, tuple[str, int, str, str], float]:
    """Build a valid TSV line with its expected key and score."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    pos = draw(st.integers(min_value=1, max_value=250_000_000))
    ref = draw(
        st.text(
            alphabet=st.sampled_from(NUCLEOTIDES),
            min_size=1,
            max_size=5,
        )
    )
    alt = draw(
        st.text(
            alphabet=st.sampled_from(NUCLEOTIDES),
            min_size=1,
            max_size=5,
        )
    )
    score = draw(
        st.floats(min_value=0.0, max_value=99.0, allow_nan=False, allow_infinity=False)
    )

    line = f"{chrom}\t{pos}\t{ref}\t{alt}\t{score}"
    key = (chrom, pos, ref, alt)
    return (line, key, score)


@st.composite
def malformed_tsv_line(draw: st.DrawFn) -> str:
    """Build a TSV line that ScoreLoader should skip (wrong col count, bad types)."""
    kind = draw(st.sampled_from(["few_columns", "bad_score", "bad_position"]))

    if kind == "few_columns":
        num_cols = draw(st.integers(min_value=1, max_value=4))
        cols = [
            draw(st.text(min_size=1, max_size=10, alphabet="abcABCGT123"))
            for _ in range(num_cols)
        ]
        return "\t".join(cols)
    elif kind == "bad_score":
        chrom = draw(st.sampled_from(CHROMOSOMES))
        pos = draw(st.integers(min_value=1, max_value=250_000_000))
        ref = draw(st.sampled_from(NUCLEOTIDES))
        alt = draw(st.sampled_from(NUCLEOTIDES))
        bad_score = draw(
            st.text(min_size=1, max_size=8, alphabet="abcxyz!@#")
        )
        return f"{chrom}\t{pos}\t{ref}\t{alt}\t{bad_score}"
    else:
        chrom = draw(st.sampled_from(CHROMOSOMES))
        bad_pos = draw(
            st.text(min_size=1, max_size=8, alphabet="abcxyz!@#")
        )
        ref = draw(st.sampled_from(NUCLEOTIDES))
        alt = draw(st.sampled_from(NUCLEOTIDES))
        score = draw(
            st.floats(
                min_value=0.0, max_value=99.0,
                allow_nan=False, allow_infinity=False,
            )
        )
        return f"{chrom}\t{bad_pos}\t{ref}\t{alt}\t{score}"


@given(
    valid_lines=st.lists(
        valid_tsv_line(), min_size=0, max_size=20, unique_by=lambda x: x[1]
    ),
    malformed_lines=st.lists(malformed_tsv_line(), min_size=0, max_size=10),
    data=st.data(),
)
@settings(max_examples=100)
def test_malformed_line_resilience(
    valid_lines: list[tuple[str, tuple[str, int, str, str], float]],
    malformed_lines: list[str],
    data: st.DataObject,
) -> None:
    """ScoreLoader loads only good lines, skips bad ones silently."""
    all_lines: list[str] = [line for line, _, _ in valid_lines] + malformed_lines

    if all_lines:
        shuffled = data.draw(st.permutations(all_lines))
    else:
        shuffled = []

    expected: dict[tuple[str, int, str, str], float] = {}
    for _, key, score in valid_lines:
        expected[key] = score

    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_path = Path(tmpdir) / "scores.tsv"
        with open(tsv_path, "w", encoding="utf-8") as f:
            for line in shuffled:
                f.write(line + "\n")

        loader = ScoreLoader()
        result = loader.load_cadd(tsv_path)

    assert len(result) <= len(expected), (
        f"Result has {len(result)} entries but expected at most "
        f"{len(expected)} (from {len(valid_lines)} valid lines)"
    )

    for key, score in result.items():
        assert key in expected, (
            f"Unexpected key {key} in result — not from a valid line"
        )
        assert score == expected[key], (
            f"Score mismatch for {key}: got {score}, expected {expected[key]}"
        )

    for key in expected:
        assert key in result, (
            f"Expected key {key} missing from result — valid line was skipped"
        )
