"""VCF file parser using pysam for memory-efficient streaming.

Wraps pysam's VariantFile to stream Variant records one at a time from
.vcf or .vcf.gz files. Compressed files require a corresponding .tbi
tabix index.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator, Optional

import pysam

from vartriage.io.exceptions import ParseError
from vartriage.models.variant import Variant


class VCFParser:
    """Stream Variant records from a VCF or compressed VCF file.

    Uses pysam's VariantFile for memory-efficient streaming. Records
    are yielded one at a time without loading the entire file into memory.

    Parameters
    ----------
    file_path : Path
        Path to a .vcf or .vcf.gz file. Compressed files require
        a corresponding .tbi tabix index.

    Raises
    ------
    FileNotFoundError
        If the file or required .tbi index does not exist.
    ParseError
        If the VCF header is malformed or cannot be parsed.

    Examples
    --------
    >>> with VCFParser(Path("sample.vcf")) as parser:
    ...     for variant in parser:
    ...         print(variant.chrom, variant.pos)
    """

    def __init__(self, file_path: Path) -> None:
        self._file_path = Path(file_path)
        self._vcf: Optional[pysam.VariantFile] = None
        self._closed: bool = False

        self._validate_file_exists()
        self._check_tabix_index()
        self._open_and_validate_header()

    def _validate_file_exists(self) -> None:
        """Check the VCF file exists and is readable."""
        if not self._file_path.exists():
            raise FileNotFoundError(
                f"VCF file not found: {self._file_path}"
            )
        if not self._file_path.is_file():
            raise FileNotFoundError(
                f"Path is not a file: {self._file_path}"
            )

    def _check_tabix_index(self) -> None:
        """For .vcf.gz files, verify a .tbi tabix index exists."""
        if self._file_path.suffix == ".gz" or str(self._file_path).endswith(".vcf.gz"):
            tbi_path = Path(str(self._file_path) + ".tbi")
            if not tbi_path.exists():
                raise FileNotFoundError(
                    f"Tabix index not found: {tbi_path}. "
                    f"Compressed VCF files require a .tbi index file."
                )

    def _open_and_validate_header(self) -> None:
        """Open the VCF file via pysam and validate header structure."""
        try:
            self._vcf = pysam.VariantFile(str(self._file_path), "r")
        except (ValueError, OSError) as exc:
            raise ParseError(
                line_number=1,
                detail=f"Failed to open VCF file: {exc}",
            ) from exc

        header = self._vcf.header
        if header is None:
            raise ParseError(
                line_number=1,
                detail="VCF file has no valid header",
            )

        header_str = str(header)

        if "##fileformat=" not in header_str:
            raise ParseError(
                line_number=1,
                detail="Missing ##fileformat declaration in VCF header",
            )

        if "#CHROM" not in header_str:
            raise ParseError(
                line_number=1,
                detail="Missing #CHROM column header line in VCF header",
            )

    def __iter__(self) -> Iterator[Variant]:
        """Yield Variant records one at a time from the VCF file.

        Yields
        ------
        Variant
            Parsed variant record with all mandatory VCF fields.

        Raises
        ------
        ParseError
            If a data line cannot be parsed due to missing mandatory
            columns or invalid field values.
        """
        if self._vcf is None or self._closed:
            return

        line_number = _count_header_lines(self._vcf)

        while True:
            try:
                record = next(self._vcf)
            except StopIteration:
                break
            except (OSError, ValueError) as exc:
                line_number += 1
                raise ParseError(
                    line_number=line_number,
                    detail=f"Failed to parse VCF data line: {exc}",
                ) from exc
            line_number += 1
            yield self._record_to_variant(record, line_number)

    def _record_to_variant(
        self, record: pysam.VariantRecord, line_number: int
    ) -> Variant:
        """Convert a pysam VariantRecord to our Variant dataclass.

        Parameters
        ----------
        record : pysam.VariantRecord
            Raw record from pysam iteration.
        line_number : int
            1-based line number for error reporting.

        Returns
        -------
        Variant
            Immutable variant record.

        Raises
        ------
        ParseError
            If a required field is missing or has an invalid value.
        """
        try:
            chrom = record.contig
        except (AttributeError, TypeError) as exc:
            raise ParseError(
                line_number=line_number,
                field="CHROM",
                detail=f"Invalid or missing CHROM value: {exc}",
            ) from exc

        if not chrom:
            raise ParseError(
                line_number=line_number,
                field="CHROM",
                detail="Empty CHROM value",
            )

        try:
            pos = record.pos
        except (AttributeError, TypeError) as exc:
            raise ParseError(
                line_number=line_number,
                field="POS",
                detail=f"Invalid POS value: {exc}",
            ) from exc

        if pos is None or pos < 1:
            raise ParseError(
                line_number=line_number,
                field="POS",
                detail=f"POS must be a positive integer, got {pos}",
            )

        variant_id: Optional[str] = record.id if record.id else None

        try:
            ref = record.ref
        except (AttributeError, TypeError) as exc:
            raise ParseError(
                line_number=line_number,
                field="REF",
                detail=f"Invalid REF value: {exc}",
            ) from exc

        if not ref:
            raise ParseError(
                line_number=line_number,
                field="REF",
                detail="Empty REF allele",
            )

        try:
            alts = record.alts
        except (AttributeError, TypeError) as exc:
            raise ParseError(
                line_number=line_number,
                field="ALT",
                detail=f"Invalid ALT value: {exc}",
            ) from exc

        alt = alts[0] if alts else "."

        qual: Optional[float] = None
        try:
            raw_qual = record.qual
            if raw_qual is not None:
                qual = float(raw_qual)
        except (TypeError, ValueError) as exc:
            raise ParseError(
                line_number=line_number,
                field="QUAL",
                detail=f"Non-numeric QUAL value: {exc}",
            ) from exc

        filter_status = self._extract_filter(record, line_number)
        info = self._extract_info(record, line_number)

        return Variant(
            chrom=chrom,
            pos=pos,
            id=variant_id,
            ref=ref,
            alt=alt,
            qual=qual,
            filter_status=filter_status,
            info=info,
        )

    def _extract_filter(
        self, record: pysam.VariantRecord, line_number: int
    ) -> str:
        """Extract the FILTER field value from a record.

        Parameters
        ----------
        record : pysam.VariantRecord
            The variant record.
        line_number : int
            Line number for error context.

        Returns
        -------
        str
            Filter status ("PASS", ".", or semicolon-joined filter names).
        """
        try:
            filters = list(record.filter)
        except (AttributeError, TypeError):
            return "."

        if not filters:
            return "."

        filter_keys = [f.name if hasattr(f, "name") else str(f) for f in filters]

        if "PASS" in filter_keys:
            return "PASS"

        return ";".join(filter_keys)

    def _extract_info(
        self, record: pysam.VariantRecord, line_number: int
    ) -> dict[str, Any]:
        """Extract INFO field key-value pairs from a record.

        Parameters
        ----------
        record : pysam.VariantRecord
            The variant record.
        line_number : int
            Line number for error context.

        Returns
        -------
        dict[str, Any]
            INFO key-value pairs with appropriate Python types.
        """
        info: dict[str, Any] = {}
        try:
            for key in record.info:
                value = record.info[key]
                if isinstance(value, tuple):
                    info[key] = list(value)
                else:
                    info[key] = value
        except (AttributeError, TypeError):
            pass

        return info

    def close(self) -> None:
        """Close the underlying file handle.

        Safe to call multiple times. After closing, iteration
        will yield no further records.
        """
        if self._vcf is not None and not self._closed:
            self._vcf.close()
            self._closed = True

    def __enter__(self) -> "VCFParser":
        """Enter context manager, returning self.

        Returns
        -------
        VCFParser
            This parser instance.
        """
        return self

    def __exit__(self, *exc: object) -> None:
        """Exit context manager, closing the file handle."""
        self.close()


def _count_header_lines(vcf: pysam.VariantFile) -> int:
    """Count header lines to establish a starting line number for data.

    Parameters
    ----------
    vcf : pysam.VariantFile
        An open pysam VariantFile.

    Returns
    -------
    int
        Number of header lines (## meta-info lines + the #CHROM line).
    """
    header_str = str(vcf.header)
    return header_str.count("\n")
