"""Pure-Python sorted interval tree using bisect for O(log n) lookups.

Implements the IntervalIndex protocol without external dependencies beyond
the standard library. Uses a sorted array of interval start positions with
binary search for efficient overlap queries.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

from vartriage.io.exceptions import ReferenceFileError


@dataclass(slots=True)
class GenomicInterval:
    """A single genomic interval with associated metadata.

    Parameters
    ----------
    chrom : str
        Chromosome name.
    start : int
        0-based start position (inclusive).
    end : int
        0-based end position (exclusive).
    feature_type : str
        GTF/GFF feature type (e.g., "exon", "CDS", "transcript").
    gene_name : str
        Gene name from the annotation.
    transcript_id : str
        Transcript identifier.
    strand : str
        Strand orientation ("+" or "-").
    """

    chrom: str
    start: int
    end: int
    feature_type: str
    gene_name: str
    transcript_id: str
    strand: str


@dataclass
class _ChromIndex:
    """Per-chromosome sorted interval index.

    Stores intervals sorted by start position for binary search lookups.
    """

    starts: list[int] = field(default_factory=list)
    ends: list[int] = field(default_factory=list)
    intervals: list[GenomicInterval] = field(default_factory=list)
    _sorted: bool = False

    def add(self, interval: GenomicInterval) -> None:
        """Add an interval to this chromosome index."""
        self.starts.append(interval.start)
        self.ends.append(interval.end)
        self.intervals.append(interval)
        self._sorted = False

    def finalize(self) -> None:
        """Sort intervals by start position for binary search."""
        if self._sorted:
            return
        indices = sorted(range(len(self.starts)), key=lambda i: (self.starts[i], self.ends[i]))
        self.starts = [self.starts[i] for i in indices]
        self.ends = [self.ends[i] for i in indices]
        self.intervals = [self.intervals[i] for i in indices]
        self._sorted = True

    def query(self, pos_start: int, pos_end: int) -> list[GenomicInterval]:
        """Find all intervals overlapping the given range [pos_start, pos_end).

        Uses binary search on sorted start positions to find candidate
        intervals, then filters by end position.

        Parameters
        ----------
        pos_start : int
            0-based start position (inclusive).
        pos_end : int
            0-based end position (exclusive).

        Returns
        -------
        list[GenomicInterval]
            All intervals that overlap the query range.
        """
        if not self._sorted:
            self.finalize()

        if not self.starts:
            return []

        # Find the rightmost interval whose start < pos_end
        right_idx = bisect.bisect_left(self.starts, pos_end)

        results: list[GenomicInterval] = []
        for i in range(right_idx):
            if self.ends[i] > pos_start:
                results.append(self.intervals[i])

        return results


class SortedArrayIntervalIndex:
    """Pure-Python sorted interval tree implementing the IntervalIndex protocol.

    Uses bisect-based binary search on sorted interval start positions for
    O(log n) candidate identification, followed by linear scan for overlap
    verification. Suitable for moderate-sized gene annotations without
    requiring pyranges.

    Parameters
    ----------
    None

    Examples
    --------
    >>> from pathlib import Path
    >>> index = SortedArrayIntervalIndex()
    >>> index.load(Path("gencode.v38.annotation.gtf"))
    >>> hits = index.overlap("chr1", 12345, "A", "T")
    """

    def __init__(self) -> None:
        self._chromosomes: dict[str, _ChromIndex] = {}
        self._loaded: bool = False
        self._exon_boundaries: dict[str, list[tuple[int, int, str]]] = {}

    def load(self, annotation_path: Path) -> None:
        """Load gene annotation from a GTF/GFF file.

        Parses the file and builds per-chromosome sorted interval indices
        for exon, CDS, and transcript features. Also builds an exon boundary
        index for splice site detection.

        Parameters
        ----------
        annotation_path : Path
            Path to the GTF or GFF gene annotation file.

        Raises
        ------
        FileNotFoundError
            If the annotation file does not exist.
        ReferenceFileError
            If the file cannot be parsed as valid GTF/GFF.
        """
        if not annotation_path.exists():
            raise FileNotFoundError(
                f"Gene annotation file not found: {annotation_path}"
            )

        try:
            self._parse_gtf(annotation_path)
        except (OSError, UnicodeDecodeError) as exc:
            raise ReferenceFileError(
                f"{annotation_path}: failed to read file - {exc}"
            ) from exc

        for chrom_idx in self._chromosomes.values():
            chrom_idx.finalize()

        self._loaded = True

    def _parse_gtf(self, path: Path) -> None:
        """Parse GTF/GFF file and populate internal indices."""
        open_func = _get_open_func(path)

        with open_func(path, "rt") as fh:
            for line_num, line in enumerate(fh, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                self._parse_gtf_line(line, line_num, path)

    def _parse_gtf_line(self, line: str, line_num: int, path: Path) -> None:
        """Parse a single GTF/GFF line and add to internal indices."""
        parts = line.split("\t")
        if len(parts) < 9:
            raise ReferenceFileError(
                f"{path}: line {line_num} has {len(parts)} columns, "
                f"expected 9 for GTF/GFF format"
            )

        chrom = parts[0]
        feature_type = parts[2]
        try:
            start = int(parts[3]) - 1  # GTF is 1-based, convert to 0-based
            end = int(parts[4])  # GTF end is inclusive, but we use exclusive
        except ValueError as exc:
            raise ReferenceFileError(
                f"{path}: line {line_num} has invalid coordinates - {exc}"
            ) from exc

        strand = parts[6]
        attributes = _parse_attributes(parts[8])

        gene_name = (
            attributes.get("gene_name")
            or attributes.get("gene_id")
            or "unknown"
        )
        transcript_id = (
            attributes.get("transcript_id")
            or attributes.get("transcript_name")
            or ""
        )

        if feature_type in ("exon", "CDS", "transcript", "gene"):
            interval = GenomicInterval(
                chrom=chrom,
                start=start,
                end=end,
                feature_type=feature_type,
                gene_name=gene_name,
                transcript_id=transcript_id,
                strand=strand,
            )
            if chrom not in self._chromosomes:
                self._chromosomes[chrom] = _ChromIndex()
            self._chromosomes[chrom].add(interval)

        if feature_type == "exon":
            if chrom not in self._exon_boundaries:
                self._exon_boundaries[chrom] = []
            self._exon_boundaries[chrom].append((start, end, transcript_id))

    def overlap(self, chrom: str, pos: int, ref: str, alt: str) -> list[dict[str, Any]]:
        """Return overlapping gene regions for a variant coordinate.

        Queries the sorted interval index and determines the functional
        consequence for each overlapping transcript. The variant position
        is 1-based (VCF convention); internally converted to 0-based for
        interval queries.

        Parameters
        ----------
        chrom : str
            Chromosome name (e.g., "chr1", "1").
        pos : int
            1-based genomic position.
        ref : str
            Reference allele.
        alt : str
            Alternate allele.

        Returns
        -------
        list[dict]
            List of overlapping regions with keys 'gene_name',
            'feature_type', 'transcript_id', 'consequence', and
            'is_splice_site'. Empty list when no overlaps found.
        """
        if not self._loaded:
            return []

        # Convert 1-based VCF position to 0-based interval coordinates
        var_start = pos - 1
        var_end = var_start + max(len(ref), len(alt))

        chrom_idx = self._chromosomes.get(chrom)
        if chrom_idx is None:
            return []

        hits = chrom_idx.query(var_start, var_end)
        if not hits:
            return []

        results: list[dict[str, Any]] = []
        for interval in hits:
            is_splice = self._is_splice_site(chrom, var_start, var_end)
            consequence = _determine_consequence(
                ref=ref,
                alt=alt,
                feature_type=interval.feature_type,
                is_splice_site=is_splice,
            )
            results.append({
                "gene_name": interval.gene_name,
                "feature_type": interval.feature_type,
                "transcript_id": interval.transcript_id,
                "consequence": consequence,
                "is_splice_site": is_splice,
            })

        return results

    def _is_splice_site(self, chrom: str, var_start: int, var_end: int) -> bool:
        """Check if variant falls within 2 bases of an exon-intron junction.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        var_start : int
            0-based variant start position.
        var_end : int
            0-based variant end position (exclusive).

        Returns
        -------
        bool
            True if the variant overlaps a splice site region.
        """
        exon_bounds = self._exon_boundaries.get(chrom)
        if not exon_bounds:
            return False

        for exon_start, exon_end, _ in exon_bounds:
            # Splice site: within 2 bases of exon-intron junction
            # Donor site: last 2 bases of exon + first 2 bases of intron
            # Acceptor site: last 2 bases of intron + first 2 bases of exon
            donor_start = exon_end - 2
            donor_end = exon_end + 2
            acceptor_start = exon_start - 2
            acceptor_end = exon_start + 2

            if (var_start < donor_end and var_end > donor_start) or (
                var_start < acceptor_end and var_end > acceptor_start
            ):
                return True

        return False


def _determine_consequence(
    ref: str,
    alt: str,
    feature_type: str,
    is_splice_site: bool,
) -> str:
    """Determine functional consequence based on variant type and genomic context.

    Uses simplified logic: classifies based on variant type (SNV vs indel)
    and position relative to gene features (coding region, splice site).

    Parameters
    ----------
    ref : str
        Reference allele.
    alt : str
        Alternate allele.
    feature_type : str
        GTF feature type of the overlapping region.
    is_splice_site : bool
        Whether the variant is at a splice site.

    Returns
    -------
    str
        The FunctionalConsequence value string.
    """
    from vartriage.models.variant import FunctionalConsequence

    # Splice site takes priority if within 2 bases of junction
    if is_splice_site:
        return FunctionalConsequence.SPLICE_SITE.value

    # Only CDS features represent coding regions
    is_coding = feature_type == "CDS"

    if not is_coding:
        # Exon but not CDS (e.g., UTR) or transcript-level hit
        if feature_type in ("exon", "transcript", "gene"):
            return FunctionalConsequence.SYNONYMOUS.value
        return FunctionalConsequence.INTERGENIC.value

    # In coding region: determine consequence by variant type
    ref_len = len(ref)
    alt_len = len(alt)
    is_snv = ref_len == 1 and alt_len == 1

    if is_snv:
        # SNV in coding region - simplified: classify as Missense
        # (full codon-level analysis would require sequence context)
        return FunctionalConsequence.MISSENSE.value

    # Insertion or deletion
    length_diff = alt_len - ref_len

    if length_diff == 0:
        # MNV (multi-nucleotide variant) in coding region
        return FunctionalConsequence.MISSENSE.value

    if length_diff % 3 != 0:
        # Not divisible by 3: frameshift
        return FunctionalConsequence.FRAMESHIFT.value

    # Divisible by 3: in-frame
    if length_diff > 0:
        return FunctionalConsequence.IN_FRAME_INSERTION.value
    return FunctionalConsequence.IN_FRAME_DELETION.value


def _parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse GTF/GFF attribute column into key-value pairs.

    Handles both GTF format (key "value";) and GFF3 format (key=value;).

    Parameters
    ----------
    attr_string : str
        The 9th column of a GTF/GFF line.

    Returns
    -------
    dict[str, str]
        Parsed attribute key-value pairs.
    """
    is_gtf = '="' not in attr_string and ('" ' in attr_string or '";' in attr_string)
    if is_gtf:
        return _parse_gtf_attributes(attr_string)
    return _parse_gff3_attributes(attr_string)


def _parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """Parse GTF-style attributes: key "value"; key "value";"""
    attributes: dict[str, str] = {}
    for item in attr_string.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) == 2:
            attributes[parts[0]] = parts[1].strip('"').strip("'")
    return attributes


def _parse_gff3_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3-style attributes: key=value;key=value;"""
    attributes: dict[str, str] = {}
    for item in attr_string.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, _, value = item.partition("=")
        attributes[key.strip()] = value.strip().strip('"')
    return attributes


def _get_open_func(path: Path) -> Callable[..., Any]:
    """Return the appropriate file opener based on file extension.

    Parameters
    ----------
    path : Path
        File path to check.

    Returns
    -------
    Callable[..., Any]
        Either gzip.open for .gz files or builtins.open for plain text.
    """
    if path.suffix == ".gz" or str(path).endswith(".gz"):
        import gzip
        return gzip.open
    return open
