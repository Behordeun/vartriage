"""Vertebrate mitochondrial genetic code.

The mitochondrial genome uses a translation table that differs from
the standard nuclear code at 4 codons:
- TGA encodes Trp (W) instead of Stop (*)
- ATA encodes Met (M) instead of Ile (I)
- AGA encodes Stop (*) instead of Arg (R)
- AGG encodes Stop (*) instead of Arg (R)

This module provides translate_codon_mt() which checks these overrides
first, then falls through to the standard table for all other codons.
"""

from __future__ import annotations

from vartriage._internal.genetic_code import CODON_TABLE

MT_CODON_OVERRIDES: dict[str, str] = {
    "TGA": "W",  # Standard: Stop (*) -> Mito: Trp (W)
    "ATA": "M",  # Standard: Ile (I) -> Mito: Met (M)
    "AGA": "*",  # Standard: Arg (R) -> Mito: Stop (*)
    "AGG": "*",  # Standard: Arg (R) -> Mito: Stop (*)
}

# Chromosome names recognized as mitochondrial
_MT_CHROM_NAMES: frozenset[str] = frozenset({"CHRM", "MT", "M"})
# RefSeq accession for the revised Cambridge Reference Sequence (rCRS) mtDNA.
_MT_ACCESSIONS: frozenset[str] = frozenset({"NC_012920"})


def is_mitochondrial(chrom: str) -> bool:
    """Check whether a chromosome name refers to mitochondrial DNA.

    Recognizes the common naming conventions (chrM, chrMT, MT, M,
    case-insensitive) and the RefSeq mtDNA accession NC_012920. The prefix
    is matched as a literal "chr", not a character set, so an unrelated
    contig is never misrouted.

    Parameters
    ----------
    chrom
        Chromosome name from VCF (e.g., "chrM", "MT", "chr1").

    Returns
    -------
    bool
        True if the chromosome is mitochondrial.
    """
    upper = chrom.upper()
    if upper in _MT_ACCESSIONS or upper.split(".")[0] in _MT_ACCESSIONS:
        return True
    name = upper[3:] if upper.startswith("CHR") else upper
    return name in ("M", "MT")


def translate_codon_mt(codon: str) -> str:
    """Translate a 3bp DNA codon using the vertebrate mitochondrial code.

    Checks the 4 mitochondrial-specific overrides first, then falls
    through to the standard nuclear codon table for all other codons.

    Parameters
    ----------
    codon
        3-character uppercase DNA string (e.g., "TGA").

    Returns
    -------
    str
        Single-letter amino acid code, "*" for stop codons, or "?"
        for invalid/ambiguous codons containing N.
    """
    upper = codon.upper()
    override = MT_CODON_OVERRIDES.get(upper)
    if override is not None:
        return override
    return CODON_TABLE.get(upper, "?")
