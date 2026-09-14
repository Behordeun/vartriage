"""Regression guard: _ChromIndex.query overlap lookups stay sub-linear.

Before the max-end segment tree, query() scanned every interval whose start
preceded the query end, making a single lookup O(n) in the number of features
on the chromosome and the whole annotation pass O(variants x features). These
tests pin the sub-linear behavior so the quadratic cannot silently return.
"""

from __future__ import annotations

import time

from vartriage._internal.interval_tree import GenomicInterval, _ChromIndex


def _build_index(n: int) -> _ChromIndex:
    """A chromosome index of n short, non-overlapping intervals, 100 bp apart."""
    idx = _ChromIndex()
    for i in range(n):
        start = i * 100
        idx.add(
            GenomicInterval(
                chrom="chr1",
                start=start,
                end=start + 50,
                feature_type="exon",
                gene_name=f"G{i}",
                transcript_id=f"TX{i}",
                strand="+",
            )
        )
    idx.finalize()
    return idx


def _time_queries(idx: _ChromIndex, positions: list[int]) -> float:
    t0 = time.perf_counter()
    for p in positions:
        idx.query(p, p + 1)
    return time.perf_counter() - t0


def test_query_returns_only_true_overlaps_at_scale() -> None:
    idx = _build_index(50_000)
    # A point inside interval 40_000 ([4_000_000, 4_000_050)) hits exactly one.
    hits = idx.query(4_000_010, 4_000_011)
    assert len(hits) == 1
    assert hits[0].gene_name == "G40000"
    # A point in a 50 bp gap between intervals hits nothing.
    assert idx.query(4_000_075, 4_000_076) == []


def test_query_scaling_is_sublinear() -> None:
    """A 10x larger index must not make per-query cost grow ~10x.

    Queries hit the far end of each index (worst case for the old prefix scan,
    where right_idx ~= n). With the segment tree the cost is O(log n + k), so
    the 10x-larger index should stay well under a 4x time blow-up. The 4x bound
    is deliberately loose to absorb timing noise while still failing hard on a
    return to O(n) (which would be ~10x).
    """
    small_n, large_n = 20_000, 200_000
    small = _build_index(small_n)
    large = _build_index(large_n)

    # Query positions near the high end of each index (many starts precede them).
    small_pos = [(small_n - 1) * 100 + 10 for _ in range(2000)]
    large_pos = [(large_n - 1) * 100 + 10 for _ in range(2000)]

    # Warm once (touch code paths) before timing.
    _time_queries(small, small_pos[:50])
    _time_queries(large, large_pos[:50])

    t_small = _time_queries(small, small_pos)
    t_large = _time_queries(large, large_pos)

    # Each worst-case query returns exactly one interval regardless of size.
    assert len(small.query(small_pos[0], small_pos[0] + 1)) == 1
    assert len(large.query(large_pos[0], large_pos[0] + 1)) == 1

    # Sub-linear: 10x the intervals must cost far less than 10x per query.
    # Guard against timing flake when both are tiny.
    if t_small > 1e-3:
        assert t_large < 4.0 * t_small, (
            f"query scaling looks linear: {small_n} intervals took {t_small:.4f}s, "
            f"{large_n} took {t_large:.4f}s (ratio {t_large / t_small:.1f}x)"
        )
