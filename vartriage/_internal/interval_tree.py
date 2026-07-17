"""Pure-Python sorted interval tree using bisect for O(log n) lookups.

Implements the IntervalIndex protocol without external dependencies beyond
the standard library. Uses a sorted array of interval start positions with
binary search for efficient overlap queries.
"""

from __future__ import annotations

import bisect
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Optional

from vartriage._internal.cache import try_load_cache, try_write_cache
from vartriage.io.exceptions import ReferenceFileError

if TYPE_CHECKING:
    from vartriage.annotation.codon_resolver import CodonResolver
    from vartriage.annotation.transcript_index import TranscriptCDSIndex

logger = logging.getLogger(__name__)


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
        indices = sorted(
            range(len(self.starts)), key=lambda i: (self.starts[i], self.ends[i])
        )
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
        self._codon_resolver: Optional["CodonResolver"] = None
        self._transcript_index: Optional["TranscriptCDSIndex"] = None

    def load(self, annotation_path: Path) -> None:
        """Load gene annotation from a GTF/GFF file.

        Parses the file and builds per-chromosome sorted interval indices
        for exon, CDS, and transcript features. Also builds an exon boundary
        index for splice site detection.

        Uses pickle-based caching to skip re-parsing on subsequent loads
        when the source file has not changed.

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
        cached = try_load_cache(annotation_path)
        if cached is not None:
            self._chromosomes = cached["chromosomes"]
            self._exon_boundaries = cached["exon_boundaries"]
            self._loaded = True
            logger.debug("Loaded interval index from cache for %s", annotation_path)
            return

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

        # Finalize transcript index for codon resolution
        if self._transcript_index is not None:
            self._transcript_index.finalize()

        self._loaded = True

        try_write_cache(
            annotation_path,
            {
                "chromosomes": self._chromosomes,
                "exon_boundaries": self._exon_boundaries,
            },
        )

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
            attributes.get("gene_name") or attributes.get("gene_id") or "unknown"
        )
        transcript_id = (
            attributes.get("transcript_id") or attributes.get("transcript_name") or ""
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

        if feature_type == "CDS" and transcript_id:
            self._index_cds_exon(parts, transcript_id, gene_name, chrom, start, end, strand)

        if feature_type == "exon":
            self._exon_boundaries.setdefault(chrom, []).append((start, end, transcript_id))

    def _index_cds_exon(
        self,
        parts: list[str],
        transcript_id: str,
        gene_name: str,
        chrom: str,
        start: int,
        end: int,
        strand: str,
    ) -> None:
        """Add a CDS exon to the TranscriptCDSIndex."""
        if self._transcript_index is None:
            from vartriage.annotation.transcript_index import TranscriptCDSIndex
            self._transcript_index = TranscriptCDSIndex()
        try:
            frame = int(parts[7]) if parts[7] != "." else 0
        except (ValueError, IndexError):
            frame = 0
        self._transcript_index.add_cds_exon( 
            transcript_id=transcript_id,
            gene_name=gene_name,
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
            frame=frame,
        )

    def set_codon_resolver(self, resolver: CodonResolver) -> None:
        """Attach a CodonResolver for amino acid-level consequence calling.

        When set, SNVs in CDS regions use proper codon resolution
        instead of the positional heuristic. Requires a reference FASTA.
        """
        self._codon_resolver = resolver

    @property
    def transcript_index(self) -> Optional[TranscriptCDSIndex]:
        """Access the TranscriptCDSIndex built during GTF parsing."""
        return self._transcript_index

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
                codon_resolver=self._codon_resolver,
                chrom=chrom,
                pos=pos,
                transcript_id=interval.transcript_id,
            )
            results.append(
                {
                    "gene_name": interval.gene_name,
                    "feature_type": interval.feature_type,
                    "transcript_id": interval.transcript_id,
                    "consequence": consequence,
                    "is_splice_site": is_splice,
                }
            )

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


def _snv_consequence(
    codon_resolver: object,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    transcript_id: str,
) -> str:
    """Consequence for a coding SNV, using codon resolution when available."""
    from vartriage.models.variant import FunctionalConsequence

    if codon_resolver is not None and chrom and pos > 0:
        context = codon_resolver.resolve(chrom, pos, ref, alt, transcript_id or None)  # type: ignore[attr-defined]
        if context is not None:
            if context.is_nonsense:
                return FunctionalConsequence.NONSENSE.value
            if context.is_synonymous:
                return FunctionalConsequence.SYNONYMOUS.value
            return FunctionalConsequence.MISSENSE.value
    return FunctionalConsequence.MISSENSE.value


def _indel_consequence(ref: str, alt: str) -> str:
    """Consequence for a coding insertion, deletion, or MNV."""
    from vartriage.models.variant import FunctionalConsequence

    length_diff = len(alt) - len(ref)
    if length_diff == 0:
        return FunctionalConsequence.MISSENSE.value
    if length_diff % 3 != 0:
        return FunctionalConsequence.FRAMESHIFT.value
    if length_diff > 0:
        return FunctionalConsequence.IN_FRAME_INSERTION.value
    return FunctionalConsequence.IN_FRAME_DELETION.value


def _determine_consequence(
    ref: str,
    alt: str,
    feature_type: str,
    is_splice_site: bool,
    codon_resolver: Optional[CodonResolver] = None,
    chrom: str = "",
    pos: int = 0,
    transcript_id: str = "",
) -> str:
    """Determine functional consequence based on variant type and genomic context.

    When a CodonResolver is provided, SNVs in CDS regions get proper
    codon-level analysis (checking the actual amino acid change) instead
    of the simplified positional heuristic.

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
    codon_resolver : object, optional
        CodonResolver instance for amino acid-level consequence calling.
        When None, falls back to the positional heuristic.
    chrom : str
        Chromosome (needed for codon resolution).
    pos : int
        1-based position (needed for codon resolution).
    transcript_id : str
        Transcript ID for targeted resolution.

    Returns
    -------
    str
        The FunctionalConsequence value string.
    """
    from vartriage.models.variant import FunctionalConsequence

    if is_splice_site:
        return FunctionalConsequence.SPLICE_SITE.value

    if feature_type != "CDS":
        if feature_type in ("exon", "transcript", "gene"):
            return FunctionalConsequence.SYNONYMOUS.value
        return FunctionalConsequence.INTERGENIC.value

    if len(ref) == 1 and len(alt) == 1:
        return _snv_consequence(codon_resolver, chrom, pos, ref, alt, transcript_id)

    return _indel_consequence(ref, alt)


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
