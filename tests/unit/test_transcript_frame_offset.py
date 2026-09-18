"""frame_offset must come from the translation-start exon, not GTF line order.

The start codon is in the lowest-coordinate CDS exon on the plus strand and the
highest-coordinate CDS exon on the minus strand. finalize() must pick the frame
of that exon regardless of the order exons were added, so a coordinate-sorted
GTF and a transcript-ordered GTF produce the same reading frame.
"""

from __future__ import annotations

from vartriage.annotation.transcript_index import TranscriptCDSIndex


def _build(strand: str, exons: list[tuple[int, int, int]]) -> int:
    """exons = list of (start, end, frame) added in the given order; return frame_offset."""
    idx = TranscriptCDSIndex()
    for start, end, frame in exons:
        idx.add_cds_exon(
            transcript_id="T1",
            gene_name="GENE",
            chrom="chr1",
            start=start,
            end=end,
            strand=strand,
            frame=frame,
        )
    idx.finalize()
    return idx.get_transcript("T1").frame_offset


def test_plus_strand_uses_lowest_coord_exon_frame() -> None:
    # Start codon in the lowest-coordinate exon (frame 0); added in coord order.
    assert _build("+", [(100, 130, 0), (200, 230, 2), (300, 330, 1)]) == 0


def test_plus_strand_independent_of_add_order() -> None:
    # Same transcript, exons added in reverse/shuffled order: still the low-coord frame.
    assert _build("+", [(300, 330, 1), (100, 130, 0), (200, 230, 2)]) == 0


def test_minus_strand_uses_highest_coord_exon_frame() -> None:
    # Negative strand: start codon in the highest-coordinate exon (frame 0),
    # listed first in GENCODE transcript order.
    assert _build("-", [(300, 330, 0), (200, 230, 2), (100, 130, 1)]) == 0


def test_minus_strand_independent_of_add_order() -> None:
    # A coordinate-sorted GTF would add the negative-strand exons ascending,
    # putting the translation-END exon first; frame_offset must still be the
    # highest-coordinate (start) exon's frame, not the first-added one.
    assert _build("-", [(100, 130, 1), (200, 230, 2), (300, 330, 0)]) == 0


def test_minus_strand_nonzero_start_frame_preserved() -> None:
    # Highest-coord exon carries frame 2; that must be the offset regardless of order.
    assert _build("-", [(100, 130, 0), (300, 330, 2), (200, 230, 1)]) == 2
