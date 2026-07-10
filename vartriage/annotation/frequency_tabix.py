"""Tabix-indexed VCF backend for gnomAD frequency lookups.

Queries a bgzipped, tabix-indexed gnomAD VCF file on-the-fly using
pysam's TabixFile. The entire VCF is never loaded into memory, making
this backend suitable for whole-genome frequency references of
arbitrary size.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pysam

from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.warnings import MissingDataWarning

logger = logging.getLogger(__name__)


class TabixFrequencyDatabase:
    """On-the-fly gnomAD frequency lookup via tabix-indexed VCF.

    Satisfies the FrequencyDatabase protocol. Queries are resolved
    per-region using pysam.TabixFile — the entire VCF is never
    loaded into memory.

    Attributes
    ----------
    warnings : list[MissingDataWarning]
        Accumulated warnings for variants not found.
    """

    def __init__(self) -> None:
        self._tabix: Optional[pysam.TabixFile] = None
        self.warnings: list[MissingDataWarning] = []

    def load(self, reference_path: Path) -> None:
        """Open tabix-indexed VCF and verify index exists.

        Parameters
        ----------
        reference_path : Path
            Path to .vcf.bgz or .vcf.gz file.

        Raises
        ------
        ReferenceFileError
            If the file or its .tbi index is missing/unreadable.
        """
        if not reference_path.exists():
            raise ReferenceFileError(
                f"{reference_path}: file not found"
            )

        index_path = Path(str(reference_path) + ".tbi")
        if not index_path.exists():
            raise ReferenceFileError(
                f"{reference_path}: tabix index file not found "
                f"(expected {index_path})"
            )

        try:
            self._tabix = pysam.TabixFile(str(reference_path))
        except OSError as exc:
            raise ReferenceFileError(
                f"{reference_path}: cannot open tabix file: {exc}"
            ) from exc

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[float]]:
        """Query allele frequencies for a batch of variants.

        For each (chrom, pos, ref, alt):
        1. Query tabix for the region chrom:pos-1 to pos (0-based)
        2. Parse matching VCF records
        3. For multiallelic sites, find the matching ALT allele
        4. Extract the corresponding AF value

        Parameters
        ----------
        variants : list[tuple[str, int, str, str]]
            (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[Optional[float]]
            Allele frequencies, positionally matched. None for
            variants not found.
        """
        results: list[Optional[float]] = []

        for chrom, pos, ref, alt in variants:
            af = self._lookup_single(chrom, pos, ref, alt)
            if af is None:
                self.warnings.append(
                    MissingDataWarning(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        source="gnomAD",
                        reason="not_found",
                    )
                )
            results.append(af)

        return results

    def _lookup_single(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Optional[float]:
        """Query tabix for a single variant's allele frequency."""
        if self._tabix is None:
            return None

        try:
            records = self._tabix.fetch(
                chrom, pos - 1, pos
            )
        except ValueError:
            # Chromosome not in the tabix index
            return None

        for record_line in records:
            af = self._parse_af_from_record(record_line, ref, alt)
            if af is not None:
                return af

        return None

    def _parse_af_from_record(
        self, record_line: str, ref: str, alt: str
    ) -> Optional[float]:
        """Extract AF for a specific alt allele from a VCF record.

        Handles multiallelic records by splitting ALT and AF fields
        and matching by index.

        Parameters
        ----------
        record_line : str
            Raw VCF record line from tabix query.
        ref : str
            Reference allele to match.
        alt : str
            Alternate allele to match.

        Returns
        -------
        Optional[float]
            Allele frequency for the matching alt, or None.
        """
        fields = record_line.split("\t")
        if len(fields) < 8:
            return None

        record_ref = fields[3]
        record_alts = fields[4].split(",")
        info_field = fields[7]

        if record_ref != ref:
            return None

        if alt not in record_alts:
            return None

        alt_index = record_alts.index(alt)

        af_value = self._extract_af_from_info(info_field)
        if af_value is None:
            return None

        af_values = af_value.split(",")

        if alt_index >= len(af_values):
            logger.warning(
                "AF field has fewer values than ALT alleles "
                "in record: %s",
                record_line[:100],
            )
            return None

        try:
            return float(af_values[alt_index])
        except ValueError:
            logger.warning(
                "Malformed AF value '%s' in record: %s",
                af_values[alt_index],
                record_line[:100],
            )
            return None

    def _extract_af_from_info(self, info_field: str) -> Optional[str]:
        """Parse the AF key from the INFO column string.

        Parameters
        ----------
        info_field : str
            The semicolon-separated INFO column from a VCF record.

        Returns
        -------
        Optional[str]
            The raw AF value string, or None if AF is not present.
        """
        for entry in info_field.split(";"):
            if entry.startswith("AF="):
                return entry[3:]
        return None
