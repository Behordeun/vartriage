"""Batch annotation hot paths iterate columns, not row Series.

Guards assign_batch, gene_names_batch and cds_overlaps_batch against the
pandas iterrows regression. The join between a variant batch and the gene
model produces roughly twenty hits per variant, so building a Series per
hit (iterrows) dominated the annotation stage. These tests pin correctness
on a small gene model and pin the scaling so a return to iterrows fails.
"""

from __future__ import annotations

import time

import pytest

pd = pytest.importorskip("pandas")
pr = pytest.importorskip("pyranges")

from vartriage.annotation.consequence_pyranges import PyRangesConsequenceAnnotator
from vartriage.models.variant import Variant


def _annotator(gtf_rows: list[dict[str, object]]) -> PyRangesConsequenceAnnotator:
    gr = pr.PyRanges(pd.DataFrame(gtf_rows))
    annotator = PyRangesConsequenceAnnotator.__new__(PyRangesConsequenceAnnotator)
    from vartriage.annotation.consequence_pyranges import PyRangesIntervalIndex

    idx = PyRangesIntervalIndex()
    idx._gr = gr
    idx._exon_gr = gr[gr.df["Feature"] == "exon"] if "Feature" in gr.df else gr
    idx._loaded = True
    annotator._index = idx
    annotator._annotation_path = None  # type: ignore[assignment]
    annotator._fallback = None
    return annotator


def _variant(chrom: str, pos: int, ref: str = "A", alt: str = "T") -> Variant:
    return Variant(
        chrom=chrom, pos=pos, id=None, ref=ref, alt=alt, qual=None, filter_status="PASS"
    )


def _gene_model() -> list[dict[str, object]]:
    # One gene with a transcript, two exons and a CDS, so a single variant
    # produces several overlapping hits (the condition that made iterrows slow).
    base = {"Chromosome": "chr1", "gene_name": "GENEA", "transcript_id": "T1"}
    return [
        {**base, "Feature": "gene", "Start": 100, "End": 400},
        {**base, "Feature": "transcript", "Start": 100, "End": 400},
        {**base, "Feature": "exon", "Start": 100, "End": 200},
        {**base, "Feature": "exon", "Start": 300, "End": 400},
        {**base, "Feature": "CDS", "Start": 120, "End": 180},
    ]


def test_assign_batch_runs_and_classifies_overlap() -> None:
    annotator = _annotator(_gene_model())
    variants = [_variant("chr1", 150), _variant("chr1", 350), _variant("chr2", 150)]
    out = annotator.assign_batch(variants)
    assert len(out) == 3
    # chr2 variant overlaps nothing -> intergenic.
    from vartriage.models.variant import FunctionalConsequence

    assert out[2] == FunctionalConsequence.INTERGENIC


def test_gene_names_batch_maps_overlapping_gene() -> None:
    annotator = _annotator(_gene_model())
    variants = [_variant("chr1", 150), _variant("chr2", 150)]
    names = annotator.gene_names_batch(variants)
    assert names[0] == "GENEA"
    assert names[1] is None


def test_cds_overlaps_batch_lists_cds_transcripts() -> None:
    annotator = _annotator(_gene_model())
    variants = [_variant("chr1", 150), _variant("chr1", 350)]
    cds = annotator.cds_overlaps_batch(variants)
    assert cds[0] == ["T1"]  # 150 is inside the CDS 120..180
    assert cds[1] == []  # 350 is exonic but not CDS


def test_batch_annotation_stays_fast_with_many_hits() -> None:
    # Many overlapping features per position force a high hit-to-variant
    # ratio, the exact condition iterrows handled badly. 300 overlapping
    # transcripts across 3000 variants is a fraction of a second with column
    # iteration and seconds with per-row Series construction.
    rows: list[dict[str, object]] = [
        {
            "Chromosome": "chr1",
            "gene_name": "GENEA",
            "transcript_id": f"T{t}",
            "Feature": "exon",
            "Start": 100,
            "End": 100_000,
        }
        for t in range(300)
    ]
    annotator = _annotator(rows)
    variants = [_variant("chr1", 100 + i) for i in range(3000)]

    start = time.perf_counter()
    annotator.assign_batch(variants)
    annotator.gene_names_batch(variants)
    elapsed = time.perf_counter() - start

    assert elapsed < 2.0, (
        f"batch annotation over a high hit ratio took {elapsed:.2f}s; a per-row "
        f"Series iteration has regressed"
    )
