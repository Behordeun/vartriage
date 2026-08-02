"""ClinVar pathogenic variant index for PS1/PM5 evidence criteria.

Provides protein-level lookups against ClinVar Pathogenic assertions:
- PS1: same amino acid change at same position, different nucleotide change
- PM5: different missense change at same amino acid position

The index is keyed on (gene, amino_acid_position) and stores the set of
known pathogenic amino acid substitutions at each position. This enables
O(1) lookup for both criteria.

Input file format (TSV):
    gene    position    ref_aa    alt_aa    chrom    pos    ref    alt    significance

This file is generated from ClinVar VCF using the prepare_clinvar_protein_index.py
script, which filters for Pathogenic/Likely_Pathogenic missense variants and
annotates them with protein changes via VEP or codon resolution.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from vartriage._internal.path_safety import resolve_path
from vartriage.io.exceptions import ReferenceFileError

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class PathogenicMissense:
    """A known pathogenic missense variant from ClinVar.

    Stored in the protein index to support PS1/PM5 lookups.
    """

    gene: str
    position: int
    ref_aa: str
    alt_aa: str
    chrom: str
    genomic_pos: int
    ref_allele: str
    alt_allele: str


@dataclass
class ProteinPositionEntry:
    """All known pathogenic missense changes at a single amino acid position."""

    gene: str
    position: int
    variants: list[PathogenicMissense] = field(default_factory=list)

    def has_same_aa_change(self, ref_aa: str, alt_aa: str) -> bool:
        """True if any variant at this position produces the same substitution."""
        return any(v.ref_aa == ref_aa and v.alt_aa == alt_aa for v in self.variants)

    def has_different_aa_change(self, ref_aa: str, alt_aa: str) -> bool:
        """True if any variant at this position produces a different substitution.

        Requires at least one pathogenic missense that differs from the query.
        The reference AA must match (same position), but the alt must differ.
        """
        return any(
            v.ref_aa == ref_aa and v.alt_aa != alt_aa for v in self.variants
        )

    def is_different_nucleotide(
        self,
        chrom: str,
        genomic_pos: int,
        ref_allele: str,
        alt_allele: str,
        ref_aa: str,
        alt_aa: str,
    ) -> bool:
        """True if at least one entry has the same AA change via a different nucleotide.

        For PS1: same amino acid change but achieved via different codons.
        Only considers variants that produce the same ref_aa -> alt_aa substitution.
        """
        for v in self.variants:
            if v.ref_aa != ref_aa or v.alt_aa != alt_aa:
                continue
            # Same AA change — check if nucleotide differs
            same_nucleotide = (
                v.chrom == chrom
                and v.genomic_pos == genomic_pos
                and v.ref_allele == ref_allele
                and v.alt_allele == alt_allele
            )
            if not same_nucleotide:
                return True
        return False


class ClinVarProteinIndex:
    """Index of ClinVar pathogenic missense variants keyed by protein position.

    Supports PS1 and PM5 evidence assignment in the ACMG classifier.

    Usage
    -----
    >>> index = ClinVarProteinIndex()
    >>> index.load(Path("clinvar_protein_index.tsv"))
    >>> index.check_ps1("BRCA1", 1775, "M", "R", "chr17", 43091429, "T", "G")
    True  # different nucleotide produces same M1775R, known pathogenic
    """

    def __init__(self) -> None:
        # Key: (gene, amino_acid_position) -> ProteinPositionEntry
        self._index: dict[tuple[str, int], ProteinPositionEntry] = {}
        self._loaded: bool = False

    @classmethod
    def from_variants(cls, variants: list[PathogenicMissense]) -> "ClinVarProteinIndex":
        """Build an in-memory index from a list of PathogenicMissense entries.

        Useful for testing and programmatic construction without a TSV file.
        """
        instance = cls()
        for v in variants:
            key = (v.gene, v.position)
            if key not in instance._index:
                instance._index[key] = ProteinPositionEntry(
                    gene=v.gene, position=v.position
                )
            instance._index[key].variants.append(v)
        instance._loaded = True
        return instance

    @property
    def is_loaded(self) -> bool:
        """True if the index has been loaded from a reference file."""
        return self._loaded

    @property
    def variant_count(self) -> int:
        """Total pathogenic missense variants in the index."""
        return sum(len(entry.variants) for entry in self._index.values())

    def load(
        self,
        reference_path: Path,
        strict: bool = False,
        max_skipped_lines: int = 100,
    ) -> None:
        """Load the ClinVar protein index from TSV.

        Parameters
        ----------
        reference_path : Path
            Path to the protein index TSV with columns:
            gene, position, ref_aa, alt_aa, chrom, pos, ref, alt, significance
        strict : bool
            When True, any malformed line raises ReferenceFileError.
        max_skipped_lines : int
            Maximum number of malformed lines tolerated before raising
            ReferenceFileError. Ignored when strict is True. Default 100.

        Raises
        ------
        ReferenceFileError
            If the file cannot be read or parsed, or if malformed line
            thresholds are exceeded.
        """
        reference_path = resolve_path(reference_path)

        if not reference_path.exists():
            raise ReferenceFileError(
                f"{reference_path}: ClinVar protein index file not found"
            )

        try:
            self._parse_tsv(reference_path, strict=strict, max_skipped=max_skipped_lines)
        except ReferenceFileError:
            raise
        except Exception as exc:
            raise ReferenceFileError(
                f"{reference_path}: failed to parse protein index: {exc}"
            ) from exc

        self._loaded = True
        logger.info(
            "Loaded ClinVar protein index: %d positions, %d variants",
            len(self._index),
            self.variant_count,
        )

    def check_ps1(
        self,
        gene: str,
        aa_position: int,
        ref_aa: str,
        alt_aa: str,
        chrom: str,
        genomic_pos: int,
        ref_allele: str,
        alt_allele: str,
    ) -> bool:
        """Check if PS1 criterion is met.

        PS1: same amino acid change as a previously established pathogenic
        variant, but arising from a different nucleotide change.

        Parameters
        ----------
        gene : str
            Gene symbol.
        aa_position : int
            1-based amino acid position.
        ref_aa : str
            Reference amino acid (single letter).
        alt_aa : str
            Alternate amino acid (single letter).
        chrom : str
            Chromosome of the query variant.
        genomic_pos : int
            Genomic position of the query variant.
        ref_allele : str
            Reference allele of the query variant.
        alt_allele : str
            Alternate allele of the query variant.

        Returns
        -------
        bool
            True if a pathogenic variant with the same amino acid change
            exists at this position via a different nucleotide change.
        """
        if not self._loaded:
            return False

        entry = self._index.get((gene, aa_position))
        if entry is None:
            return False

        return entry.is_different_nucleotide(
            chrom=chrom,
            genomic_pos=genomic_pos,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
        )

    def check_pm5(
        self,
        gene: str,
        aa_position: int,
        ref_aa: str,
        alt_aa: str,
    ) -> bool:
        """Check if PM5 criterion is met.

        PM5: novel missense change at an amino acid position where a
        DIFFERENT missense change has been determined to be pathogenic.

        Parameters
        ----------
        gene : str
            Gene symbol.
        aa_position : int
            1-based amino acid position.
        ref_aa : str
            Reference amino acid (single letter).
        alt_aa : str
            Alternate amino acid being evaluated.

        Returns
        -------
        bool
            True if a pathogenic variant with a different amino acid change
            exists at this position.
        """
        if not self._loaded:
            return False

        entry = self._index.get((gene, aa_position))
        if entry is None:
            return False

        return entry.has_different_aa_change(ref_aa, alt_aa)

    def _parse_tsv(
        self, path: Path, strict: bool = False, max_skipped: int = 100
    ) -> None:
        """Parse the protein index TSV into the lookup dictionary."""
        skipped = 0

        with open(path, encoding="utf-8") as fh:
            for line_num, line in enumerate(fh, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                # Skip header (works regardless of position after comment lines)
                if stripped.lower().startswith("gene"):
                    continue

                parts = stripped.split("\t")
                if len(parts) < 8:
                    skipped += 1
                    if strict:
                        raise ReferenceFileError(
                            f"{path}:{line_num}: expected 8+ columns, got {len(parts)}"
                        )
                    logger.warning(
                        "Skipping line %d: expected 8+ columns, got %d",
                        line_num,
                        len(parts),
                    )
                    if skipped > max_skipped:
                        raise ReferenceFileError(
                            f"{path}: too many malformed lines ({skipped}), "
                            f"exceeds threshold of {max_skipped}"
                        )
                    continue

                gene = parts[0]
                try:
                    position = int(parts[1])
                except ValueError:
                    skipped += 1
                    if strict:
                        raise ReferenceFileError(
                            f"{path}:{line_num}: non-integer position '{parts[1]}'"
                        )
                    logger.warning(
                        "Skipping line %d: non-integer position '%s'",
                        line_num,
                        parts[1],
                    )
                    if skipped > max_skipped:
                        raise ReferenceFileError(
                            f"{path}: too many malformed lines ({skipped}), "
                            f"exceeds threshold of {max_skipped}"
                        )
                    continue

                ref_aa = parts[2]
                alt_aa = parts[3]
                chrom = parts[4]

                try:
                    genomic_pos = int(parts[5])
                except ValueError:
                    skipped += 1
                    if strict:
                        raise ReferenceFileError(
                            f"{path}:{line_num}: non-integer genomic position '{parts[5]}'"
                        )
                    logger.warning(
                        "Skipping line %d: non-integer genomic position '%s'",
                        line_num,
                        parts[5],
                    )
                    if skipped > max_skipped:
                        raise ReferenceFileError(
                            f"{path}: too many malformed lines ({skipped}), "
                            f"exceeds threshold of {max_skipped}"
                        )
                    continue

                ref_allele = parts[6]
                alt_allele = parts[7]

                variant = PathogenicMissense(
                    gene=gene,
                    position=position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    chrom=chrom,
                    genomic_pos=genomic_pos,
                    ref_allele=ref_allele,
                    alt_allele=alt_allele,
                )

                key = (gene, position)
                if key not in self._index:
                    self._index[key] = ProteinPositionEntry(
                        gene=gene, position=position
                    )
                self._index[key].variants.append(variant)
