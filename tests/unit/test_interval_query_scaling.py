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

    # Median of repeated samples resists scheduler/GC/CPU-frequency noise that a
    # single wall-clock reading would be at the mercy of.
    t_small = sorted(_time_queries(small, small_pos) for _ in range(5))[2]
    t_large = sorted(_time_queries(large, large_pos) for _ in range(5))[2]

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


def test_query_rebuilds_tree_for_legacy_finalized_index() -> None:
    """A _ChromIndex restored from a cache pickled before the max-end tree
    existed is _sorted=True but has no tree attribute; query() must rebuild it
    rather than raise AttributeError."""
    idx = _build_index(1000)
    # Simulate an object unpickled from an older build: sorted, tree attrs gone.
    del idx._max_end_tree
    del idx._tree_size
    assert idx._sorted is True

    hits = idx.query(500 * 100 + 10, 500 * 100 + 11)
    assert len(hits) == 1
    assert hits[0].gene_name == "G500"


def test_splice_cache_reset_across_reloads() -> None:
    """Reloading a different annotation into the same index must not serve
    splice windows derived from the first annotation."""
    import tempfile
    from pathlib import Path

    from vartriage._internal.interval_tree import SortedArrayIntervalIndex

    gtf_a = (
        "chr1\ttest\texon\t1000\t1200\t.\t+\t.\t"
        'gene_id "A"; transcript_id "TA"; gene_name "A";\n'
    )
    gtf_b = (
        "chr1\ttest\texon\t5000\t5200\t.\t+\t.\t"
        'gene_id "B"; transcript_id "TB"; gene_name "B";\n'
    )

    def _write(text: str) -> Path:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as tmp:
            tmp.write(text)
            return Path(tmp.name)

    path_a, path_b = _write(gtf_a), _write(gtf_b)
    try:
        index = SortedArrayIntervalIndex()
        index.load(path_a)
        # Prime the splice-window cache against annotation A (donor near 1200).
        assert index._is_splice_site("chr1", 1199, 1200) is True

        index.load(path_b)
        # After reload, A's boundaries are gone and B's are present.
        assert index._is_splice_site("chr1", 1199, 1200) is False
        assert index._is_splice_site("chr1", 4999, 5000) is True
    finally:
        path_a.unlink(missing_ok=True)
        path_b.unlink(missing_ok=True)
