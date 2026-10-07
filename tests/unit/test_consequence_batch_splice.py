"""Batch splice-site detection in the pyranges consequence annotator.

Guards _find_splice_positions against both a correctness regression and the
per-variant exon rescan that made it scale with variants times exons. The
batch result must match the single-variant _check_splice_site decision for
every variant, and must stay fast when the exon model is large.
"""

from __future__ import annotations

import time

import pytest

pd = pytest.importorskip("pandas")
pr = pytest.importorskip("pyranges")

from vartriage.annotation.consequence_pyranges import (
    PyRangesConsequenceAnnotator,
    PyRangesIntervalIndex,
)
from vartriage.models.variant import Variant


def _annotator_with_exons(
    rows: list[tuple[str, int, int]],
) -> PyRangesConsequenceAnnotator:
    idx = PyRangesIntervalIndex()
    df = pd.DataFrame(
        {
            "Chromosome": [r[0] for r in rows],
            "Start": [r[1] for r in rows],
            "End": [r[2] for r in rows],
        }
    )
    idx._exon_gr = pr.PyRanges(df)
    # A gene frame is required for _loaded; reuse the exon frame, the splice
    # path only reads _exon_gr.
    idx._gr = pr.PyRanges(df.copy())
    idx._loaded = True

    annotator = PyRangesConsequenceAnnotator.__new__(PyRangesConsequenceAnnotator)
    annotator._index = idx
    annotator._annotation_path = None  # type: ignore[assignment]
    annotator._fallback = None
    return annotator


def _variant(chrom: str, pos: int) -> Variant:
    # 1-based pos; a 1bp SNV so Start=pos-1, End=pos in _build_query.
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=None,
        filter_status="PASS",
    )


def _query_df(annotator: PyRangesConsequenceAnnotator, variants: list[Variant]):
    query_df, _ = annotator._build_query(variants)
    return query_df


def test_batch_splice_matches_single_variant_decision() -> None:
    exons = [("chr1", 100, 200), ("chr1", 500, 600), ("chr2", 100, 200)]
    annotator = _annotator_with_exons(exons)

    # A mix of donor-edge, acceptor-edge, deep-exon, intronic, wrong-chrom.
    variants = [
        _variant("chr1", 200),  # donor edge of exon 100..200
        _variant("chr1", 101),  # acceptor edge of exon 100..200
        _variant("chr1", 150),  # deep in exon, not splice
        _variant("chr1", 600),  # donor edge of exon 500..600
        _variant("chr1", 350),  # intronic, not splice
        _variant("chr2", 100),  # acceptor edge on chr2
        _variant("chr3", 150),  # chromosome with no exons
    ]

    batch = annotator._find_splice_positions(_query_df(annotator, variants))

    # Oracle: the already-vectorized single-variant check, per variant.
    expected = {
        i
        for i, v in enumerate(variants)
        if annotator._index._check_splice_site(
            v.chrom, v.pos - 1, (v.pos - 1) + max(len(v.ref), len(v.alt))
        )
    }

    assert batch == expected, f"batch {sorted(batch)} != per-variant {sorted(expected)}"


def test_batch_splice_edge_sweep_matches_single_variant() -> None:
    # Sweep every position across both exon boundaries. This is the guard that
    # catches an off-by-one in the join window: the two outermost splice bases
    # (two before the exon start, two after the exon end) are easy to drop.
    annotator = _annotator_with_exons([("chr1", 100, 200)])
    variants = [_variant("chr1", pos) for pos in range(90, 211)]

    batch = annotator._find_splice_positions(_query_df(annotator, variants))
    expected = {
        i
        for i, v in enumerate(variants)
        if annotator._index._check_splice_site(
            v.chrom, v.pos - 1, (v.pos - 1) + max(len(v.ref), len(v.alt))
        )
    }

    assert batch == expected, (
        f"edge sweep mismatch: batch-only {sorted(batch - expected)}, "
        f"oracle-only {sorted(expected - batch)}"
    )


def test_batch_splice_empty_when_no_exons() -> None:
    annotator = _annotator_with_exons([])
    variants = [_variant("chr1", 150)]
    assert annotator._find_splice_positions(_query_df(annotator, variants)) == set()


def test_batch_splice_handles_empty_variant_set() -> None:
    annotator = _annotator_with_exons([("chr1", 100, 200)])
    empty = pd.DataFrame({"Chromosome": [], "Start": [], "End": [], "_idx": []})
    assert annotator._find_splice_positions(empty) == set()


def test_batch_splice_stays_fast_with_many_variants_and_exons() -> None:
    # The regression this guards: filtering a large exon frame once per variant
    # is O(variants x exons). At 200k exons by 5k variants the per-variant
    # rescan takes multiple seconds and climbs with the product; a single
    # interval join resolves the same batch in a fraction of a second. The
    # threshold sits below the measured per-variant cost so the guard fails if
    # the loop ever returns.
    exons = [("chr1", i * 300, i * 300 + 150) for i in range(200_000)]
    annotator = _annotator_with_exons(exons)
    variants = [_variant("chr1", i * 150) for i in range(5_000)]
    query_df = _query_df(annotator, variants)

    start = time.perf_counter()
    annotator._find_splice_positions(query_df)
    elapsed = time.perf_counter() - start

    assert elapsed < 1.0, (
        f"5000 variants over 200k exons took {elapsed:.2f}s; the per-variant "
        f"exon rescan has regressed"
    )
