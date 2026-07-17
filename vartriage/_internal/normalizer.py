"""Variant normalization: left-align and trim using reference FASTA.

Normalizes VCF-style variants to their canonical left-aligned
representation before database lookups. Without normalization, the
same indel can be represented at different positions depending on
the VCF caller's alignment strategy, causing silent lookup failures
against gnomAD, ClinVar, and score files.

Algorithm based on Tan et al. 2015 (Bioinformatics):
1. Right-trim: remove shared suffix between REF and ALT
2. Left-trim: remove shared prefix (keep at least 1bp anchor for indels)
3. Left-align: shift indels left until the rightmost base no longer
   matches the preceding reference base
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class VariantNormalizer:
    """Left-align and trim variants using a reference FASTA.

    Thread-safe (pysam.FastaFile is read-only after open).

    Parameters
    ----------
    fasta_path
        Path to an indexed reference FASTA (.fa with .fai).
    """

    def __init__(self, fasta_path: Path) -> None:
        import pysam

        self._fasta = pysam.FastaFile(str(fasta_path))
        self._fasta_path = fasta_path

    def normalize(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> tuple[str, int, str, str]:
        """Normalize a variant to its canonical left-aligned form.

        Parameters
        ----------
        chrom
            Chromosome name.
        pos
            1-based position (VCF convention).
        ref
            Reference allele.
        alt
            Alternate allele.

        Returns
        -------
        tuple[str, int, str, str]
            Normalized (chrom, pos, ref, alt). The chrom is unchanged.
            For SNVs and already-normalized variants, returns the input
            unchanged.
        """
        # SNVs don't need normalization
        if len(ref) == 1 and len(alt) == 1:
            return (chrom, pos, ref, alt)

        ref_upper = ref.upper()
        alt_upper = alt.upper()

        # Step 1: Right-trim shared suffix
        ref_upper, alt_upper = self._right_trim(ref_upper, alt_upper)

        # Step 2: Left-trim shared prefix (keep 1bp for indels)
        ref_upper, alt_upper, pos = self._left_trim(ref_upper, alt_upper, pos)

        # Step 3: Left-align indels
        ref_upper, alt_upper, pos = self._left_align(
            chrom, pos, ref_upper, alt_upper
        )

        return (chrom, pos, ref_upper, alt_upper)

    def _right_trim(self, ref: str, alt: str) -> tuple[str, str]:
        """Remove shared suffix bases from both alleles."""
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
        return ref, alt

    def _left_trim(self, ref: str, alt: str, pos: int) -> tuple[str, str, int]:
        """Remove shared prefix bases, keeping at least 1bp for indels."""
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1
        return ref, alt, pos

    def _left_align(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> tuple[str, str, int]:
        """Shift indels left until they can't shift further.

        An indel can shift left when the rightmost base of the allele
        matches the reference base immediately preceding the variant.
        """
        # Only align if one allele is longer (insertion or deletion)
        if len(ref) == len(alt):
            return ref, alt, pos

        # Limit iterations to prevent infinite loops on pathological cases
        max_shifts = 1000
        shifts = 0

        while pos > 1 and shifts < max_shifts:
            # Fetch the base preceding the current position
            prev_pos_0 = pos - 2  # 0-based position of the preceding base
            try:
                prev_base = self._fasta.fetch(chrom, prev_pos_0, prev_pos_0 + 1).upper()
            except (ValueError, IndexError):
                break

            if not prev_base or len(prev_base) != 1:
                break

            # Check if rightmost bases match the preceding reference base
            if ref[-1] != prev_base or alt[-1] != prev_base:
                break

            # Shift left: prepend the preceding base and remove the last base
            ref = prev_base + ref[:-1]
            alt = prev_base + alt[:-1]
            pos -= 1
            shifts += 1

        return ref, alt, pos

    def close(self) -> None:
        """Close the FASTA file handle."""
        self._fasta.close()
