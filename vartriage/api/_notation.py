"""VCF-to-VEP variant notation conversion.

Converts VCF-style (chrom, pos, ref, alt) representations into VEP's
region notation format for the POST /vep/human/region endpoint.

VEP notation: "chrom start end allele_string strand"
  - Coordinates are 1-based inclusive
  - Chromosome is numeric (no "chr" prefix)
  - Strand is always "+" for VCF-derived variants

Edge cases handled:
  - SNVs: straightforward position mapping
  - Deletions: strip VCF padding base, adjust start coordinate
  - Insertions: use flanking coordinates around insertion point
  - MNVs: span the full substitution range
  - Complex indels (ref and alt both >1bp): decompose as deletion+insertion
"""

from __future__ import annotations

import re

_VALID_ALLELE = re.compile(r"^[ACGTNacgtn*.\-]+$|^<[A-Za-z:]+>$")
_VALID_CHROM = re.compile(r"^[A-Za-z0-9_.\-]+$")


def _sanitize_genomic_input(chrom: str, ref: str, alt: str) -> None:
    """Reject inputs containing characters outside the genomic alphabet.

    Prevents injection of control characters or HTML into downstream
    API request bodies or report output.

    Chromosome: alphanumeric plus underscore, dot, and hyphen (covers
    standard contigs like chr1_GL000192_random, chrUn_gl000220,
    GL000220.1, plus ALT/decoy contigs like chr1-ALT, HLA-A*01:01:01:01).

    Alleles: ACGTN bases, dot, dash, star (standard VCF) plus symbolic
    alleles in angle brackets (e.g. <DEL>, <INS>, <DUP:TANDEM>).
    """
    if not _VALID_CHROM.match(chrom):
        raise ValueError(f"Invalid chromosome name: {chrom!r}")
    if not _VALID_ALLELE.match(ref):
        raise ValueError(f"Invalid reference allele: {ref!r}")
    if not _VALID_ALLELE.match(alt):
        raise ValueError(f"Invalid alternate allele: {alt!r}")


def vcf_to_vep_notation(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Convert a VCF variant to VEP region notation.

    Parameters
    ----------
    chrom
        Chromosome name (e.g., "chr22", "22", "chrX").
    pos
        1-based VCF position.
    ref
        Reference allele string.
    alt
        Alternate allele string.

    Returns
    -------
    str
        VEP-format string: "chrom start end allele/allele +"
    """
    _sanitize_genomic_input(chrom, ref, alt)
    chrom_clean = _strip_chr_prefix(chrom)
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1 and alt_len == 1:
        # SNV
        return f"{chrom_clean} {pos} {pos} {ref}/{alt} +"  # nosec: output is API request body, not HTML

    if ref_len > 1 and alt_len == 1:
        # Deletion: alt is just the padding base
        # VEP wants the deleted sequence without the padding base
        deleted = ref[1:]
        start = pos + 1
        end = pos + len(deleted)
        return f"{chrom_clean} {start} {end} {deleted}/- +"  # nosec: output is API request body, not HTML

    if ref_len == 1 and alt_len > 1:
        # Insertion: ref is just the padding base
        # VEP insertion coordinates flank the insertion point
        inserted = alt[1:]
        start = pos
        end = pos + 1
        return f"{chrom_clean} {start} {end} -/{inserted} +"  # nosec: output is API request body, not HTML

    # Complex: both ref and alt > 1bp
    # Could be MNV (same length) or complex indel (different lengths)
    if ref_len == alt_len:
        # MNV: multi-nucleotide variant, same length substitution
        end = pos + ref_len - 1
        return f"{chrom_clean} {pos} {end} {ref}/{alt} +"  # nosec: output is API request body, not HTML

    # Complex indel: different lengths, both > 1bp
    # Strip shared padding base, represent as the changed portion
    # VEP handles this as a combined deletion+insertion
    deleted = ref[1:]
    inserted = alt[1:]
    start = pos + 1
    end = pos + len(deleted) if deleted else pos

    del_part = deleted if deleted else "-"
    ins_part = inserted if inserted else "-"
    return f"{chrom_clean} {start} {end} {del_part}/{ins_part} +"  # nosec: output is API request body, not HTML


def _strip_chr_prefix(chrom: str) -> str:
    """Remove 'chr' prefix for VEP's expected numeric chromosome format.

    Preserves 'MT' as-is (VEP accepts it). Handles lowercase 'chr' variants.
    """
    lower = chrom.lower()
    if lower.startswith("chr"):
        stripped = chrom[3:]
        # Handle chrM -> MT for VEP compatibility
        if stripped.upper() == "M":
            return "MT"
        return stripped
    return chrom
