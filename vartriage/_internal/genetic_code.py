"""Standard genetic code translation table.

Maps RNA codons (as DNA: T instead of U) to single-letter amino acid
codes. Stop codons map to "*". Provides translation and reverse
complement utilities for codon-level consequence resolution.
"""

from __future__ import annotations

CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def translate_codon(codon: str) -> str:
    """Translate a 3bp DNA codon to a single-letter amino acid.

    Parameters
    ----------
    codon
        3-character uppercase DNA string (e.g., "ATG").

    Returns
    -------
    str
        Single letter amino acid code, or "*" for stop codons.
        Returns "?" for invalid/ambiguous codons (containing N).
    """
    return CODON_TABLE.get(codon.upper(), "?")


def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence.

    Parameters
    ----------
    seq
        DNA string (uppercase A/T/C/G/N).

    Returns
    -------
    str
        Reverse complemented sequence.
    """
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(seq.upper()))
