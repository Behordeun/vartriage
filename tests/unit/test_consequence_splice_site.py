"""Splice-site detection in the pyranges consequence index.

Guards the vectorized _check_splice_site against both correctness regressions
and the row-by-row iteration that made it scale with exon count.
"""

from __future__ import annotations

import time

import pandas as pd
import pyranges as pr

from vartriage.annotation.consequence_pyranges import PyRangesIntervalIndex


def _index_with_exons(rows: list[tuple[str, int, int]]) -> PyRangesIntervalIndex:
    idx = PyRangesIntervalIndex()
    df = pd.DataFrame(
        {
            "Chromosome": [r[0] for r in rows],
            "Start": [r[1] for r in rows],
            "End": [r[2] for r in rows],
        }
    )
    idx._exon_gr = pr.PyRanges(df)
    idx._loaded = True
    return idx


def test_variant_at_donor_junction_is_a_splice_site() -> None:
    idx = _index_with_exons([("chr1", 100, 200)])
    # Donor window sits around the exon end (200); a variant at 199-200 overlaps it.
    assert idx._check_splice_site("chr1", 199, 200) is True


def test_variant_at_acceptor_junction_is_a_splice_site() -> None:
    idx = _index_with_exons([("chr1", 100, 200)])
    # Acceptor window sits around the exon start (100).
    assert idx._check_splice_site("chr1", 99, 100) is True


def test_variant_deep_in_exon_is_not_a_splice_site() -> None:
    idx = _index_with_exons([("chr1", 100, 200)])
    assert idx._check_splice_site("chr1", 150, 151) is False


def test_variant_on_other_chromosome_is_not_a_splice_site() -> None:
    idx = _index_with_exons([("chr1", 100, 200)])
    assert idx._check_splice_site("chr2", 199, 200) is False


def test_splice_check_stays_fast_with_many_exons() -> None:
    # A chromosome with thousands of exons must not be walked row by row per
    # query. Ten thousand exons times a thousand queries would take minutes
    # under iterrows; vectorized it is well under a second.
    exons = [("chr1", i * 300, i * 300 + 150) for i in range(10_000)]
    idx = _index_with_exons(exons)

    start = time.perf_counter()
    for pos in range(0, 1000):
        idx._check_splice_site("chr1", pos * 300, pos * 300 + 1)
    elapsed = time.perf_counter() - start

    assert elapsed < 5.0, f"1000 splice queries over 10k exons took {elapsed:.1f}s"
