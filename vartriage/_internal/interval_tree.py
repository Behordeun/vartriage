"""Pure-Python sorted interval tree using bisect for O(log n) lookups.

Implements the IntervalIndex protocol without external dependencies beyond
the standard library. Uses a sorted array of interval start positions with
binary search for efficient overlap queries.
"""

from __future__ import annotations

import bisect
import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from vartriage._internal.cache import try_load_cache, try_write_cache
from vartriage.io.exceptions import ReferenceFileError

if TYPE_CHECKING:
    from vartriage.annotation.codon_resolver import CodonContext, CodonResolver
    from vartriage.annotation.transcript_index import TranscriptCDSIndex

logger = logging.getLogger(__name__)

# Splice windows span 2 bp on each side of an exon boundary, so every window is
# exactly 4 bp wide. Used to bound the binary-search band in _is_splice_site.
_SPLICE_WINDOW_WIDTH = 4


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
    _max_end_tree: list[int] = field(default_factory=list)
    _tree_size: int = 0

    def add(self, interval: GenomicInterval) -> None:
        """Add an interval to this chromosome index."""
        self.starts.append(interval.start)
        self.ends.append(interval.end)
        self.intervals.append(interval)
        self._sorted = False

    def finalize(self) -> None:
        """Sort intervals by start position and build the max-end segment tree."""
        if self._sorted:
            return
        indices = sorted(
            range(len(self.starts)), key=lambda i: (self.starts[i], self.ends[i])
        )
        self.starts = [self.starts[i] for i in indices]
        self.ends = [self.ends[i] for i in indices]
        self.intervals = [self.intervals[i] for i in indices]
        self._build_max_end_tree()
        self._sorted = True

    def _build_max_end_tree(self) -> None:
        """Build an iterative segment tree over ends holding subtree maxima.

        Leaves live at [size, size+n); internal node i holds max(2i, 2i+1).
        Enables pruning the overlap scan to intervals whose end exceeds the
        query start, so query() is O(log n + k) instead of O(n).
        """
        n = len(self.ends)
        if n == 0:
            self._max_end_tree = []
            return
        size = 1
        while size < n:
            size *= 2
        tree = [0] * (2 * size)
        for i, e in enumerate(self.ends):
            tree[size + i] = e
        for i in range(size - 1, 0, -1):
            tree[i] = max(tree[2 * i], tree[2 * i + 1])
        self._max_end_tree = tree
        self._tree_size = size

    def query(self, pos_start: int, pos_end: int) -> list[GenomicInterval]:
        """Find all intervals overlapping the given range [pos_start, pos_end).

        Binary-searches the sorted starts for the candidate prefix (start <
        pos_end), then walks that prefix guided by a max-end segment tree so
        only intervals whose end exceeds pos_start are visited. O(log n + k).

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

        # Self-heal a cache pickled before the max-end tree existed: such an
        # object is _sorted=True but has no tree, so finalize() would be skipped.
        if not self._max_end_tree:
            self._build_max_end_tree()

        # Candidate prefix: intervals whose start < pos_end.
        right_idx = bisect.bisect_left(self.starts, pos_end)
        if right_idx <= 0:
            return []

        results: list[GenomicInterval] = []
        self._collect_overlaps(1, 0, self._tree_size, 0, right_idx, pos_start, results)
        return results

    def _collect_overlaps(
        self,
        node: int,
        node_lo: int,
        node_hi: int,
        q_lo: int,
        q_hi: int,
        pos_start: int,
        out: list[GenomicInterval],
    ) -> None:
        """Descend the segment tree over indices [q_lo, q_hi), collecting leaves
        whose interval end > pos_start. Subtrees whose max end <= pos_start are
        pruned."""
        if node_hi <= q_lo or q_hi <= node_lo:
            return
        if self._max_end_tree[node] <= pos_start:
            return
        if node_hi - node_lo == 1:
            # Leaf: node_lo is the interval index (may be padding beyond n).
            if node_lo < len(self.ends) and self.ends[node_lo] > pos_start:
                out.append(self.intervals[node_lo])
            return
        mid = (node_lo + node_hi) // 2
        self._collect_overlaps(2 * node, node_lo, mid, q_lo, q_hi, pos_start, out)
        self._collect_overlaps(2 * node + 1, mid, node_hi, q_lo, q_hi, pos_start, out)


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
        self._splice_windows_cache: dict[str, tuple[list[int], list[int]]] | None = None
        self._codon_resolver: CodonResolver | None = None
        self._transcript_index: TranscriptCDSIndex | None = None

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
            self._index_cds_exon(
                parts, transcript_id, gene_name, chrom, start, end, strand
            )

        if feature_type == "exon":
            self._exon_boundaries.setdefault(chrom, []).append(
                (start, end, transcript_id)
            )

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
    def transcript_index(self) -> TranscriptCDSIndex | None:
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
            consequence, codon_context = _determine_consequence(
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
                    "codon_context": codon_context,
                }
            )

        return results

    def _is_splice_site(self, chrom: str, var_start: int, var_end: int) -> bool:
        """Check if variant falls within 2 bases of an exon-intron junction.

        Splice windows (donor: exon_end +/- 2; acceptor: exon_start +/- 2) are
        precomputed once per chromosome into a start-sorted array and queried by
        binary search. Every window is exactly 4 bp wide, so a match can only
        involve windows whose start lies in [var_start - 4, var_end); bisecting
        to that band makes the check O(log n + k) with a tiny constant k instead
        of scanning every exon on the chromosome.

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
        windows = self._splice_window_starts.get(chrom)
        if not windows:
            return False
        starts, ends = windows

        # Windows are 4 bp wide, so any overlapping window has its start in
        # [var_start - 4, var_end). Bisect the upper bound, then walk back over
        # the bounded band checking the end condition.
        hi = bisect.bisect_left(starts, var_end)
        lo_bound = var_start - _SPLICE_WINDOW_WIDTH
        for i in range(hi - 1, -1, -1):
            w_start = starts[i]
            if w_start < lo_bound:
                break
            if ends[i] > var_start:
                return True
        return False

    @property
    def _splice_window_starts(
        self,
    ) -> dict[str, tuple[list[int], list[int]]]:
        """Per-chromosome start-sorted splice windows, built once on first use.

        Derived from ``_exon_boundaries``: each exon contributes an acceptor
        window (exon_start +/- 2) and a donor window (exon_end +/- 2).
        """
        cached = self._splice_windows_cache
        if cached is not None:
            return cached

        built: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, bounds in self._exon_boundaries.items():
            windows: list[tuple[int, int]] = []
            for exon_start, exon_end, _ in bounds:
                windows.append((exon_start - 2, exon_start + 2))
                windows.append((exon_end - 2, exon_end + 2))
            windows.sort()
            starts = [w[0] for w in windows]
            ends = [w[1] for w in windows]
            built[chrom] = (starts, ends)

        self._splice_windows_cache = built
        return built


def _snv_consequence(
    codon_resolver: object,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    transcript_id: str,
) -> tuple[str, CodonContext | None]:
    """Consequence for a coding SNV, using codon resolution when available.

    Returns a tuple of (consequence_string, CodonContext_or_None).
    """
    from vartriage.models.variant import FunctionalConsequence

    if codon_resolver is not None and chrom and pos > 0:
        context = codon_resolver.resolve(chrom, pos, ref, alt, transcript_id or None)  # type: ignore[attr-defined]
        if context is not None:
            if context.is_nonsense:
                return FunctionalConsequence.NONSENSE.value, context
            if context.is_synonymous:
                return FunctionalConsequence.SYNONYMOUS.value, context
            return FunctionalConsequence.MISSENSE.value, context
    return FunctionalConsequence.MISSENSE.value, None


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
    codon_resolver: CodonResolver | None = None,
    chrom: str = "",
    pos: int = 0,
    transcript_id: str = "",
) -> tuple[str, CodonContext | None]:
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
    tuple[str, Optional[CodonContext]]
        (FunctionalConsequence value string, CodonContext or None).
    """
    from vartriage.models.variant import FunctionalConsequence

    if is_splice_site:
        return FunctionalConsequence.SPLICE_SITE.value, None

    if feature_type != "CDS":
        if feature_type in ("exon", "transcript", "gene"):
            return FunctionalConsequence.SYNONYMOUS.value, None
        return FunctionalConsequence.INTERGENIC.value, None

    if len(ref) == 1 and len(alt) == 1:
        return _snv_consequence(codon_resolver, chrom, pos, ref, alt, transcript_id)

    return _indel_consequence(ref, alt), None


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
