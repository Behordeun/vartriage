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


def test_finalize_preserves_offset_when_exon_frames_unknown() -> None:
    # A TranscriptCDS built with frame-less exons (legacy/direct construction)
    # but an explicit frame_offset must keep that offset through finalize.
    from vartriage.annotation.transcript_index import CDSExon, TranscriptCDS

    tc = TranscriptCDS(
        transcript_id="T1", gene_name="G", chrom="chr1", strand="-", frame_offset=2
    )
    tc.cds_exons.append(CDSExon(start=300, end=330))  # frame defaults to None
    tc.cds_exons.append(CDSExon(start=100, end=130))
    tc.finalize()
    assert tc.frame_offset == 2, "unknown exon frames must not overwrite frame_offset"


def test_legacy_cache_roundtrip_preserves_offset() -> None:
    # A legacy 2-field serialized index restores frame_offset directly; a later
    # finalize must not reset it to 0.
    from vartriage.annotation.transcript_index import TranscriptCDSIndex

    legacy = {
        "T1": {
            "gene_name": "G",
            "chrom": "chr1",
            "strand": "-",
            "frame_offset": 1,
            "cds_exons": [(300, 330), (100, 130)],  # 2-field, no frame
        }
    }
    idx = TranscriptCDSIndex.from_serializable(legacy)
    idx._transcripts["T1"].finalize()
    assert idx.get_transcript("T1").frame_offset == 1
