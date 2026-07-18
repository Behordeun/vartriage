"""Codon-level consequence resolution using reference FASTA.

Given a SNV overlapping a CDS region, determines the exact amino acid
change by extracting the reference codon from the genome sequence,
substituting the variant base, and translating both codons. This
replaces the positional heuristic ("SNV in CDS = Missense") with
biologically correct consequence calling.

Handles:
- Forward and reverse strand genes
- Split codons at exon junctions (codon spans two exons)
- Frame offset from GTF (first CDS exon may not start at codon boundary)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from vartriage._internal.genetic_code import (reverse_complement,
                                              translate_codon)
from vartriage.annotation.transcript_index import (TranscriptCDS,
                                                   TranscriptCDSIndex)

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class CodonContext:
    """Result of resolving a variant to its codon-level impact.

    Attributes
    ----------
    transcript_id
        Transcript where the variant was resolved.
    gene_name
        Gene symbol.
    cds_position
        0-based position within the concatenated CDS sequence.
    codon_index
        Which codon the variant falls in (cds_position // 3).
    codon_position
        Position within the codon (0, 1, or 2).
    reference_codon
        3bp reference codon (sense strand).
    altered_codon
        3bp codon after variant substitution (sense strand).
    reference_aa
        Single-letter amino acid from the reference codon.
    altered_aa
        Single-letter amino acid from the altered codon ("*" = stop).
    is_synonymous
        True when reference_aa == altered_aa.
    is_nonsense
        True when altered_aa is a stop codon ("*") and reference is not.
    """

    transcript_id: str
    gene_name: str
    cds_position: int
    codon_index: int
    codon_position: int
    reference_codon: str
    altered_codon: str
    reference_aa: str
    altered_aa: str
    is_synonymous: bool
    is_nonsense: bool


class CodonResolver:
    """Resolves coding SNVs to amino acid changes using FASTA + transcript CDS.

    Parameters
    ----------
    fasta_path
        Path to an indexed reference genome FASTA (.fa with .fai).
    transcript_index
        Pre-built TranscriptCDSIndex from GTF parsing.
    """

    def __init__(self, fasta_path: Path, transcript_index: TranscriptCDSIndex) -> None:
        import pysam

        self._fasta = pysam.FastaFile(str(fasta_path))
        self._transcripts = transcript_index

    def resolve(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        transcript_id: Optional[str] = None,
    ) -> Optional[CodonContext]:
        """Resolve a CDS-overlapping SNV to its amino acid change.

        Parameters
        ----------
        chrom
            Chromosome name.
        pos
            1-based genomic position (VCF convention).
        ref
            Reference allele (single base for SNV).
        alt
            Alternate allele (single base for SNV).
        transcript_id
            Specific transcript to resolve against. If None, uses the
            first overlapping transcript found.

        Returns
        -------
        Optional[CodonContext]
            Codon resolution result, or None if:
            - Not a SNV (ref or alt length != 1)
            - Position doesn't overlap any CDS
            - Codon extraction fails
        """
        if len(ref) != 1 or len(alt) != 1:
            return None

        # Convert to 0-based for internal use
        genomic_pos_0 = pos - 1

        # Find the transcript to resolve against
        transcript = self._find_transcript(chrom, genomic_pos_0, transcript_id)
        if transcript is None:
            return None

        # Map genomic position to CDS position
        cds_pos = transcript.genomic_to_cds_position(genomic_pos_0)
        if cds_pos is None:
            return None

        # Account for frame offset
        adjusted_cds_pos = cds_pos - transcript.frame_offset
        if adjusted_cds_pos < 0:
            return None

        codon_index = adjusted_cds_pos // 3
        codon_position = adjusted_cds_pos % 3

        # Extract the 3bp reference codon
        ref_codon = self._extract_codon(transcript, codon_index)
        if ref_codon is None or len(ref_codon) != 3:
            return None

        # For negative strand, the variant base needs to be complemented
        # because the codon is already in sense (coding) orientation
        if transcript.strand == "-":
            variant_base = reverse_complement(alt)
        else:
            variant_base = alt.upper()

        # Build the altered codon by substituting at the correct position
        alt_codon = (
            ref_codon[:codon_position] + variant_base + ref_codon[codon_position + 1 :]
        )

        # Translate both
        ref_aa = translate_codon(ref_codon)
        alt_aa = translate_codon(alt_codon)

        if ref_aa == "?" or alt_aa == "?":
            # Ambiguous codon (contains N), can't resolve
            return None

        return CodonContext(
            transcript_id=transcript.transcript_id,
            gene_name=transcript.gene_name,
            cds_position=cds_pos,
            codon_index=codon_index,
            codon_position=codon_position,
            reference_codon=ref_codon,
            altered_codon=alt_codon,
            reference_aa=ref_aa,
            altered_aa=alt_aa,
            is_synonymous=(ref_aa == alt_aa),
            is_nonsense=(alt_aa == "*" and ref_aa != "*"),
        )

    def _find_transcript(
        self, chrom: str, genomic_pos_0: int, transcript_id: Optional[str]
    ) -> Optional[TranscriptCDS]:
        """Find the transcript to resolve against."""
        if transcript_id is not None:
            return self._transcripts.get_transcript(transcript_id)

        # Find overlapping transcripts, prefer canonical (longest CDS)
        overlapping = self._transcripts.find_overlapping(chrom, genomic_pos_0)
        if not overlapping:
            return None

        # Pick the transcript with the longest CDS (proxy for canonical)
        return max(overlapping, key=lambda t: t.cds_length)

    def _extract_codon(
        self, transcript: TranscriptCDS, codon_index: int
    ) -> Optional[str]:
        """Extract a 3bp codon from the reference genome.

        Handles split codons (codon spans an exon-intron junction) by
        fetching bases from multiple exon regions and concatenating.

        Parameters
        ----------
        transcript
            The transcript providing CDS exon coordinates.
        codon_index
            Which codon to extract (0-based).

        Returns
        -------
        Optional[str]
            3bp codon string in sense (coding) orientation, or None if
            extraction fails.
        """
        # The CDS-relative start of this codon
        codon_cds_start = (codon_index * 3) + transcript.frame_offset

        # Collect 3 genomic positions corresponding to CDS positions
        # codon_cds_start, codon_cds_start+1, codon_cds_start+2
        genomic_positions = []
        for offset in range(3):
            gpos = self._cds_to_genomic(transcript, codon_cds_start + offset)
            if gpos is None:
                return None
            genomic_positions.append(gpos)

        # Fetch each base from the FASTA
        bases = []
        for gpos in genomic_positions:
            try:
                base = self._fasta.fetch(transcript.chrom, gpos, gpos + 1).upper()
                bases.append(base)
            except (ValueError, IndexError):
                return None

        codon = "".join(bases)

        # For negative strand, reverse complement to get sense orientation
        if transcript.strand == "-":
            codon = reverse_complement(codon)

        return codon

    def _cds_to_genomic(self, transcript: TranscriptCDS, cds_pos: int) -> Optional[int]:
        """Map a CDS position back to a genomic position.

        Inverse of TranscriptCDS.genomic_to_cds_position().
        """
        if transcript.strand == "+":
            cumulative = 0
            for exon in transcript.cds_exons:
                exon_len = exon.end - exon.start
                if cumulative + exon_len > cds_pos:
                    return exon.start + (cds_pos - cumulative)
                cumulative += exon_len
        else:
            # Reverse strand: CDS positions count from 3' end
            cumulative = 0
            for exon in reversed(transcript.cds_exons):
                exon_len = exon.end - exon.start
                if cumulative + exon_len > cds_pos:
                    # Position within exon, counting from the end
                    return exon.end - 1 - (cds_pos - cumulative)
                cumulative += exon_len
        return None

    def close(self) -> None:
        """Close the FASTA file handle."""
        self._fasta.close()
