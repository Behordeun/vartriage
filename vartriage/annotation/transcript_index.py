"""Transcript CDS index for codon-level consequence resolution.

Builds a per-transcript map of CDS exon coordinates from the GTF
annotation, enabling genomic-to-CDS-position conversion needed for
determining which codon a variant falls in.

The index is built during GTF parsing (piggybacks on the existing
interval tree load) and stored alongside the interval data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class CDSExon:
    """A single CDS exon in genomic coordinates (0-based, half-open)."""

    start: int
    end: int


@dataclass
class TranscriptCDS:
    """CDS structure for a single transcript.

    Stores ordered CDS exons, strand, frame offset, and provides
    genomic-to-CDS position mapping for codon resolution.

    Parameters
    ----------
    transcript_id
        Ensembl or RefSeq transcript identifier.
    gene_name
        Gene symbol.
    chrom
        Chromosome name.
    strand
        "+" or "-".
    frame_offset
        Reading frame offset (0, 1, or 2) from the GTF frame column
        of the first CDS exon. Indicates how many bases to skip before
        the first complete codon.
    """

    transcript_id: str
    gene_name: str
    chrom: str
    strand: str
    frame_offset: int = 0
    cds_exons: list[CDSExon] = field(default_factory=list)

    @property
    def cds_length(self) -> int:
        """Total CDS length in bases."""
        return sum(e.end - e.start for e in self.cds_exons)

    def genomic_to_cds_position(self, genomic_pos: int) -> Optional[int]:
        """Map a 0-based genomic position to a 0-based CDS position.

        Returns None if the position doesn't fall within any CDS exon
        of this transcript.

        Parameters
        ----------
        genomic_pos
            0-based genomic coordinate.

        Returns
        -------
        Optional[int]
            0-based position within the concatenated CDS sequence.
            Accounts for strand (negative strand counts from the last exon).
        """
        if self.strand == "+":
            return self._forward_mapping(genomic_pos)
        return self._reverse_mapping(genomic_pos)

    def _forward_mapping(self, genomic_pos: int) -> Optional[int]:
        """Map position on the forward strand."""
        cds_offset = 0
        for exon in self.cds_exons:
            if exon.start <= genomic_pos < exon.end:
                return cds_offset + (genomic_pos - exon.start)
            cds_offset += exon.end - exon.start
        return None

    def _reverse_mapping(self, genomic_pos: int) -> Optional[int]:
        """Map position on the reverse strand.

        For negative-strand genes, the CDS is read from the last exon
        to the first (in genomic order). The CDS position counts from
        the 3' end of the last CDS exon.
        """
        cds_offset = 0
        # Reverse strand: iterate exons in reverse genomic order
        for exon in reversed(self.cds_exons):
            if exon.start <= genomic_pos < exon.end:
                # Distance from the end of this exon (since we read 3'->5')
                return cds_offset + (exon.end - 1 - genomic_pos)
            cds_offset += exon.end - exon.start
        return None

    def finalize(self) -> None:
        """Sort CDS exons by genomic position. Call after all exons are added."""
        self.cds_exons.sort(key=lambda e: e.start)


class TranscriptCDSIndex:
    """Index of CDS structures for all transcripts in a GTF file.

    Built during GTF parsing by collecting CDS feature lines per
    transcript. Provides lookup by transcript_id and by genomic
    coordinate (find which transcript's CDS overlaps a position).
    """

    def __init__(self) -> None:
        self._transcripts: dict[str, TranscriptCDS] = {}
        self._finalized: bool = False

    def add_cds_exon(
        self,
        transcript_id: str,
        gene_name: str,
        chrom: str,
        start: int,
        end: int,
        strand: str,
        frame: int,
    ) -> None:
        """Register a CDS exon for a transcript.

        Parameters
        ----------
        transcript_id
            Transcript identifier.
        gene_name
            Gene symbol.
        chrom
            Chromosome.
        start
            0-based start (inclusive).
        end
            0-based end (exclusive).
        strand
            "+" or "-".
        frame
            Reading frame (0, 1, 2) from GTF column 8.
        """
        if transcript_id not in self._transcripts:
            self._transcripts[transcript_id] = TranscriptCDS(
                transcript_id=transcript_id,
                gene_name=gene_name,
                chrom=chrom,
                strand=strand,
                frame_offset=frame,
            )
        self._transcripts[transcript_id].cds_exons.append(CDSExon(start=start, end=end))
        self._finalized = False

    def finalize(self) -> None:
        """Sort all transcript CDS exons. Call after GTF parsing is complete."""
        for transcript in self._transcripts.values():
            transcript.finalize()
        self._finalized = True
        logger.debug("TranscriptCDSIndex finalized: %d transcripts", len(self._transcripts))

    def get_transcript(self, transcript_id: str) -> Optional[TranscriptCDS]:
        """Look up a transcript by ID."""
        return self._transcripts.get(transcript_id)

    def find_overlapping(self, chrom: str, pos: int) -> list[TranscriptCDS]:
        """Find all transcripts whose CDS overlaps a genomic position.

        Parameters
        ----------
        chrom
            Chromosome name.
        pos
            0-based genomic position.

        Returns
        -------
        list[TranscriptCDS]
            Transcripts with a CDS exon covering this position.
            May be empty for intergenic/intronic positions.
        """
        results: list[TranscriptCDS] = []
        for transcript in self._transcripts.values():
            if transcript.chrom != chrom:
                continue
            for exon in transcript.cds_exons:
                if exon.start <= pos < exon.end:
                    results.append(transcript)
                    break
        return results

    @property
    def transcript_count(self) -> int:
        """Number of transcripts in the index."""
        return len(self._transcripts)

    def to_serializable(self) -> dict[str, object]:
        """Serialize for caching alongside the interval tree."""
        data: dict[str, object] = {}
        for tid, tc in self._transcripts.items():
            data[tid] = {
                "gene_name": tc.gene_name,
                "chrom": tc.chrom,
                "strand": tc.strand,
                "frame_offset": tc.frame_offset,
                "cds_exons": [(e.start, e.end) for e in tc.cds_exons],
            }
        return data

    @classmethod
    def from_serializable(cls, data: dict[str, object]) -> "TranscriptCDSIndex":
        """Reconstruct from cached data."""
        index = cls()
        for tid, info in data.items():
            if not isinstance(info, dict):
                continue
            tc = TranscriptCDS(
                transcript_id=tid,
                gene_name=str(info.get("gene_name", "")),
                chrom=str(info.get("chrom", "")),
                strand=str(info.get("strand", "+")),
                frame_offset=int(info.get("frame_offset", 0)),  # type: ignore[arg-type]
            )
            for exon_tuple in info.get("cds_exons", []):  # type: ignore[union-attr]
                if isinstance(exon_tuple, (list, tuple)) and len(exon_tuple) == 2:
                    tc.cds_exons.append(CDSExon(start=exon_tuple[0], end=exon_tuple[1]))
            index._transcripts[tid] = tc
        index._finalized = True
        return index
