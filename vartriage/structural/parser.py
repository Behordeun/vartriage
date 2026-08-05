"""Structural variant VCF parser.

Streams SV records from VCF files produced by SV callers (Manta, Delly,
LUMPY, GRIDSS, etc.). Extracts SV-specific INFO fields and converts
symbolic ALT alleles (<DEL>, <DUP>, etc.) and BND notation into
StructuralVariant records.

Uses pysam for memory-efficient streaming, consistent with the
point-variant parser.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterator
from pathlib import Path

import pysam

from vartriage.io.exceptions import ParseError
from vartriage.structural.models import StructuralVariant, SVType

logger = logging.getLogger(__name__)

_SYMBOLIC_ALT_PATTERN = re.compile(r"<(\w+)>")

# BND ALT field formats per VCF 4.3 spec:
# t[p[ - forward strand, connects to p forward
# t]p] - forward strand, connects to p reverse
# ]p]t - reverse strand, connects to p forward
# [p[t - reverse strand, connects to p reverse
_BND_PATTERN = re.compile(
    r"^(?:"
    r"(?P<pre_seq>\w*)(?P<bracket1>[\[\]])(?P<mate_chrom>[^:\[\]]+):(?P<mate_pos>\d+)(?P<bracket2>[\[\]])"
    r"|"
    r"(?P<bracket3>[\[\]])(?P<mate_chrom2>[^:\[\]]+):(?P<mate_pos2>\d+)(?P<bracket4>[\[\]])(?P<post_seq>\w*)"
    r")$"
)

_SVTYPE_MAP: dict[str, SVType] = {
    "DEL": SVType.DEL,
    "DUP": SVType.DUP,
    "DUP:TANDEM": SVType.DUP,
    "DUP:DISPERSED": SVType.DUP,
    "INV": SVType.INV,
    "INS": SVType.INS,
    "INS:ME": SVType.INS,
    "INS:ME:ALU": SVType.INS,
    "INS:ME:LINE1": SVType.INS,
    "INS:ME:SVA": SVType.INS,
    "BND": SVType.BND,
    "CNV": SVType.CNV,
    "TRA": SVType.BND,
    # GATK-SV uses DEL:ME for mobile element deletions
    "DEL:ME": SVType.DEL,
    "DEL:ME:ALU": SVType.DEL,
    "DEL:ME:LINE1": SVType.DEL,
    "DEL:ME:SVA": SVType.DEL,
    # DELLY uses specific sub-types
    "DEL:TANDEM": SVType.DEL,
    "DUP:INT": SVType.DUP,
}

# Caller-specific INFO fields for SV length / end position.
# Each caller stores coordinates slightly differently.
# Priority order: caller-specific fields first, then standard fields.
_CALLER_END_FIELDS: tuple[str, ...] = (
    "END",  # VCF 4.3 standard
    "END2",  # GRIDSS
    "CHR2_POS",  # Some LUMPY versions
)

_CALLER_SVLEN_FIELDS: tuple[str, ...] = (
    "SVLEN",  # Standard
    "INSLEN",  # Manta insertion length
    "HOMLEN",  # Microhomology length (not SV length, but used as fallback)
)

_CALLER_COPY_NUMBER_FIELDS: tuple[str, ...] = (
    "CN",  # Standard
    "CNVAL",  # GATK-SV
    "TCN",  # Total copy number (some callers)
)

_CALLER_MATE_FIELDS: tuple[str, ...] = (
    "MATEID",  # Standard BND mate
    "PARID",  # GRIDSS partner ID
    "EVENT",  # Manta groups BND pairs by event
)


class SVParser:
    """Stream StructuralVariant records from a VCF file.

    Filters the VCF to SV records only (those with SVTYPE INFO field
    or symbolic ALT alleles). Non-SV records are silently skipped.

    Parameters
    ----------
    file_path : Path
        Path to a .vcf or .vcf.gz file containing SV calls.
    min_size : int
        Minimum SV length to emit. Smaller events are skipped.
        Default is 50.
    max_size : int
        Maximum SV length to emit. Larger events are skipped.
        Default is 0 (no upper limit).
    min_quality : float
        Minimum QUAL threshold. SVs below this are skipped.
        Default is 0.0 (no filtering).
    sv_types : Optional[set[SVType]]
        Restrict output to these SV types. None means all types.
    """

    def __init__(
        self,
        file_path: Path,
        min_size: int = 50,
        max_size: int = 0,
        min_quality: float = 0.0,
        sv_types: set[SVType] | None = None,
    ) -> None:
        self._file_path = Path(file_path)
        self._min_size = min_size
        self._max_size = max_size
        self._min_quality = min_quality
        self._sv_types = sv_types
        self._vcf: pysam.VariantFile | None = None
        self._closed: bool = False

        self._validate_file()
        self._open_vcf()

    def _validate_file(self) -> None:
        resolved = self._file_path.resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"VCF file not found: {self._file_path}")
        if not resolved.is_file():
            raise FileNotFoundError(f"Path is not a file: {self._file_path}")
        self._file_path = resolved

    def _open_vcf(self) -> None:
        try:
            self._vcf = pysam.VariantFile(str(self._file_path), "r")
        except (ValueError, OSError) as exc:
            raise ParseError(
                line_number=1,
                detail=f"Failed to open VCF file: {exc}",
            ) from exc

    def __iter__(self) -> Iterator[StructuralVariant]:
        """Yield StructuralVariant records, skipping non-SV entries.

        Yields
        ------
        StructuralVariant
            Parsed SV records passing size and quality filters.
        """
        if self._vcf is None or self._closed:
            return

        for record in self._vcf:
            sv = self._try_parse_sv(record)
            if sv is None:
                continue
            if not self._passes_filters(sv):
                continue
            yield sv

    def _try_parse_sv(self, record: pysam.VariantRecord) -> StructuralVariant | None:
        """Attempt to parse a VCF record as a structural variant.

        Returns None if the record is not an SV (no SVTYPE and no
        symbolic ALT allele).
        """
        sv_type = self._resolve_sv_type(record)
        if sv_type is None:
            return None

        chrom = record.contig
        start = record.pos
        variant_id = record.id if record.id != "." else None

        qual: float | None = None
        if record.qual is not None:
            qual = float(record.qual)

        filter_status = self._extract_filter(record)
        alt = self._extract_alt(record)
        info = dict(record.info)

        end = self._resolve_end(record, start, sv_type)
        svlen = self._resolve_svlen(record, start, end)
        copy_number = self._resolve_copy_number(record)
        start_ci = self._resolve_ci(record, "CIPOS")
        end_ci = self._resolve_ci(record, "CIEND")
        mate_id = self._resolve_mate_id(record)

        return StructuralVariant(
            chrom=chrom,
            start=start,
            end=end,
            sv_type=sv_type,
            id=variant_id,
            svlen=svlen,
            qual=qual,
            filter_status=filter_status,
            alt=alt,
            copy_number=copy_number,
            start_ci=start_ci,
            end_ci=end_ci,
            mate_id=mate_id,
            info=info,
        )

    def _resolve_sv_type(self, record: pysam.VariantRecord) -> SVType | None:
        """Determine SV type from SVTYPE INFO or symbolic ALT."""
        # Try INFO/SVTYPE first (most SV callers set this)
        svtype_raw = record.info.get("SVTYPE")
        if svtype_raw is not None:
            svtype_str = str(svtype_raw).upper()
            mapped = _SVTYPE_MAP.get(svtype_str)
            if mapped is not None:
                return mapped

        # Fall back to symbolic ALT allele
        alt = self._extract_alt(record)
        match = _SYMBOLIC_ALT_PATTERN.match(alt)
        if match:
            sym = match.group(1).upper()
            return _SVTYPE_MAP.get(sym)

        # Check for BND notation in ALT
        if _BND_PATTERN.match(alt):
            return SVType.BND

        return None

    def _resolve_end(
        self, record: pysam.VariantRecord, start: int, sv_type: SVType
    ) -> int:
        """Resolve end position from caller-specific or standard fields."""
        # Try each known END field across callers
        for field_name in _CALLER_END_FIELDS:
            end_val = record.info.get(field_name)
            if end_val is not None:
                return int(end_val)

        # For BND, end equals start (point breakend)
        if sv_type == SVType.BND:
            return start

        # Fall back to SVLEN variants
        for field_name in _CALLER_SVLEN_FIELDS:
            svlen_val = record.info.get(field_name)
            if svlen_val is not None:
                length = abs(int(svlen_val))
                if length > 0:
                    return start + length - 1

        # Last resort: use REF/ALT length difference for non-symbolic
        ref_len = len(record.ref) if record.ref else 1
        return start + ref_len - 1

    def _resolve_svlen(
        self,
        record: pysam.VariantRecord,
        start: int,
        end: int,
    ) -> int | None:
        """Extract SV length from caller-specific or standard fields."""
        for field_name in _CALLER_SVLEN_FIELDS:
            svlen_val = record.info.get(field_name)
            if svlen_val is not None:
                return int(svlen_val)
        # Compute from coordinates
        if end > start:
            return end - start + 1
        return None

    def _resolve_copy_number(self, record: pysam.VariantRecord) -> int | None:
        """Extract copy number from caller-specific fields."""
        for field_name in _CALLER_COPY_NUMBER_FIELDS:
            cn = record.info.get(field_name)
            if cn is not None:
                return int(cn)
        return None

    def _resolve_ci(
        self, record: pysam.VariantRecord, field_name: str
    ) -> tuple[int, int]:
        """Extract confidence interval (CIPOS or CIEND).

        These are stored as two integers: (lower, upper) offset
        from the position.
        """
        ci_val = record.info.get(field_name)
        if ci_val is not None:
            try:
                vals = list(ci_val)
                if len(vals) >= 2:
                    return (int(vals[0]), int(vals[1]))
            except (TypeError, ValueError):
                pass
        return (0, 0)

    def _resolve_mate_id(self, record: pysam.VariantRecord) -> str | None:
        """Extract mate/partner ID for BND records from caller-specific fields."""
        for field_name in _CALLER_MATE_FIELDS:
            mate = record.info.get(field_name)
            if mate is not None:
                if isinstance(mate, (list, tuple)):
                    return str(mate[0]) if mate else None
                return str(mate)
        return None

    def _extract_alt(self, record: pysam.VariantRecord) -> str:
        """Get the first ALT allele as a string."""
        alts = record.alts
        if alts and len(alts) > 0:
            return str(alts[0])
        return ""

    def _extract_filter(self, record: pysam.VariantRecord) -> str:
        """Get FILTER field value."""
        filters = list(record.filter)
        if not filters:
            return "."
        if len(filters) == 1 and filters[0] == "PASS":
            return "PASS"
        return ";".join(str(f) for f in filters)

    def _passes_filters(self, sv: StructuralVariant) -> bool:
        """Apply size, quality, and type filters."""
        # Type filter
        if self._sv_types is not None and sv.sv_type not in self._sv_types:
            return False

        # Size filter (skip for BND which are point breakends)
        if sv.sv_type != SVType.BND:
            if sv.length < self._min_size:
                return False
            if self._max_size > 0 and sv.length > self._max_size:
                return False

        # Quality filter
        return not (
            self._min_quality > 0.0
            and sv.qual is not None
            and sv.qual < self._min_quality
        )

    def close(self) -> None:
        """Release file handles."""
        if not hasattr(self, "_vcf"):
            return
        if self._vcf is not None and not self._closed:
            self._vcf.close()
            self._closed = True

    def __enter__(self) -> SVParser:
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()
