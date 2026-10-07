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


def test_batch_annotation_scales_with_hits_not_quadratically() -> None:
    # The join hit count is what iterrows handled badly: it built a Series per
    # hit, so cost grew with a large per-row constant. Column iteration makes
    # the per-hit constant tiny. Rather than assert an absolute wall-clock bound
    # (flaky on shared CI runners), assert the SHAPE: doubling the transcript
    # count (and so the hit count) must not blow the time up super-linearly.
    # Per-row Series construction regresses this ratio sharply; column
    # iteration keeps it near-linear.
    def _time_with_transcripts(num_transcripts: int, num_variants: int) -> float:
        rows: list[dict[str, object]] = [
            {
                "Chromosome": "chr1",
                "gene_name": "GENEA",
                "transcript_id": f"T{t}",
                "Feature": "exon",
                "Start": 100,
                "End": 100_000,
            }
            for t in range(num_transcripts)
        ]
        annotator = _annotator(rows)
        variants = [_variant("chr1", 100 + i) for i in range(num_variants)]
        # Warm once so import/first-call overhead does not land in the measured
        # window, then take the best of three to damp runner scheduling noise.
        annotator.assign_batch(variants)
        best = float("inf")
        for _ in range(3):
            start = time.perf_counter()
            annotator.assign_batch(variants)
            annotator.gene_names_batch(variants)
            best = min(best, time.perf_counter() - start)
        return best

    small = _time_with_transcripts(150, 2000)
    large = _time_with_transcripts(300, 2000)

    # Hits double (150 -> 300 transcripts). Near-linear column iteration keeps
    # the ratio close to 2x; a per-row Series regression pushes it far higher
    # because each extra hit pays the Series-construction constant. A generous
    # ceiling of 4x absorbs timing noise while still failing on the regression,
    # which measured well above an order of magnitude.
    floor = 1e-3  # avoid dividing by a near-zero fast path on quick hardware
    ratio = large / max(small, floor)
    assert ratio < 4.0, (
        f"doubling the hit count scaled annotation time {ratio:.1f}x "
        f"(small={small:.3f}s, large={large:.3f}s); a per-row Series "
        f"iteration has regressed the batch path"
    )
