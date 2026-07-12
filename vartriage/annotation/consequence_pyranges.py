"""Functional consequence assignment using pyranges for vectorized overlaps.

Optimized backend for consequence assignment using pyranges vectorized
genomic interval overlap queries. Activated automatically when pyranges
is installed; otherwise falls back to SortedArrayIntervalIndex.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.variant import (
    CONSEQUENCE_SEVERITY_ORDER,
    FunctionalConsequence,
    Variant,
)

try:
    import pyranges as pr
    import pandas as pd

    PYRANGES_AVAILABLE = True
except ImportError:
    PYRANGES_AVAILABLE = False


class PyRangesIntervalIndex:
    """Genomic interval index using pyranges for vectorized overlap queries.

    Implements the IntervalIndex protocol using pyranges for efficient
    batch overlap operations. Requires the pyranges optional dependency.

    Parameters
    ----------
    None

    Raises
    ------
    ImportError
        If pyranges is not installed.

    Examples
    --------
    >>> from pathlib import Path
    >>> index = PyRangesIntervalIndex()
    >>> index.load(Path("gencode.v38.annotation.gtf"))
    >>> hits = index.overlap("chr1", 12345, "A", "T")
    """

    def __init__(self) -> None:
        if not PYRANGES_AVAILABLE:
            raise ImportError(
                "pyranges is required for PyRangesIntervalIndex. "
                "Install with: pip install vartriage[accelerated]"
            )
        self._gr: Optional[pr.PyRanges] = None
        self._exon_gr: Optional[pr.PyRanges] = None
        self._loaded: bool = False

    def load(self, annotation_path: Path) -> None:
        """Load gene annotation from a GTF/GFF file into pyranges.

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
            gr = pr.read_gtf(str(annotation_path))
        except Exception as exc:
            raise ReferenceFileError(
                f"{annotation_path}: failed to parse GTF/GFF - {exc}"
            ) from exc

        # Filter to relevant features
        relevant_features = {"exon", "CDS", "transcript", "gene"}
        mask = gr.df["Feature"].isin(relevant_features)
        self._gr = gr[mask]

        # Separate exon intervals for splice site detection
        exon_mask = gr.df["Feature"] == "exon"
        self._exon_gr = gr[exon_mask]
        self._loaded = True

    def overlap(self, chrom: str, pos: int, ref: str, alt: str) -> list[dict[str, Any]]:
        """Return overlapping gene regions for a variant coordinate.

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
        list[dict[str, Any]]
            List of overlapping regions with keys 'gene_name',
            'feature_type', 'transcript_id', 'consequence', and
            'is_splice_site'. Empty list when no overlaps found.
        """
        if not self._loaded or self._gr is None:
            return []

        var_start = pos - 1
        var_end = var_start + max(len(ref), len(alt))

        query_df = pd.DataFrame({
            "Chromosome": [chrom],
            "Start": [var_start],
            "End": [var_end],
        })
        query_gr = pr.PyRanges(query_df)

        hits = self._gr.join(query_gr)
        if hits.df.empty:
            return []

        is_splice = self._check_splice_site(chrom, var_start, var_end)

        results: list[dict[str, Any]] = []
        for _, row in hits.df.iterrows():
            feature_type = row.get("Feature", "unknown")
            gene_name = row.get("gene_name", row.get("gene_id", "unknown"))
            transcript_id = row.get("transcript_id", "")

            consequence = _determine_consequence_pyranges(
                ref=ref,
                alt=alt,
                feature_type=feature_type,
                is_splice_site=is_splice,
            )

            results.append({
                "gene_name": gene_name,
                "feature_type": feature_type,
                "transcript_id": transcript_id,
                "consequence": consequence,
                "is_splice_site": is_splice,
            })

        return results

    def _check_splice_site(self, chrom: str, var_start: int, var_end: int) -> bool:
        """Check if variant falls within 2 bases of an exon-intron junction.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        var_start : int
            0-based variant start.
        var_end : int
            0-based variant end (exclusive).

        Returns
        -------
        bool
            True if variant overlaps a splice site region.
        """
        if self._exon_gr is None or self._exon_gr.df.empty:
            return False

        exon_df = self._exon_gr.df
        chrom_exons = exon_df[exon_df["Chromosome"] == chrom]

        if chrom_exons.empty:
            return False

        for _, exon in chrom_exons.iterrows():
            exon_start = exon["Start"]
            exon_end = exon["End"]

            # Donor site: 2 bases around exon end (exon-intron junction)
            donor_start = exon_end - 2
            donor_end = exon_end + 2
            # Acceptor site: 2 bases around exon start (intron-exon junction)
            acceptor_start = exon_start - 2
            acceptor_end = exon_start + 2

            if (var_start < donor_end and var_end > donor_start) or (
                var_start < acceptor_end and var_end > acceptor_start
            ):
                return True

        return False


class PyRangesConsequenceAnnotator:
    """Assign functional consequence using pyranges vectorized overlaps.

    Provides the same interface as ConsequenceAnnotator but uses pyranges
    for batch overlap operations, offering better performance on large
    datasets.

    Parameters
    ----------
    annotation_path : Path
        Path to the GTF or GFF gene annotation file.

    Raises
    ------
    ImportError
        If pyranges is not installed.
    FileNotFoundError
        If the annotation file does not exist.
    ReferenceFileError
        If the file cannot be parsed as valid GTF/GFF.
    """

    def __init__(self, annotation_path: Path) -> None:
        self._index = PyRangesIntervalIndex()
        self._index.load(annotation_path)

    def load(self, annotation_path: Path) -> None:
        """Load gene annotation from a GTF/GFF file.

        Parameters
        ----------
        annotation_path : Path
            Path to the GTF or GFF gene annotation file.
        """
        self._index.load(annotation_path)

    def overlap(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> list[dict[str, Any]]:
        """Return overlapping gene regions for a variant coordinate.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        pos : int
            1-based genomic position.
        ref : str
            Reference allele.
        alt : str
            Alternate allele.

        Returns
        -------
        list[dict[str, Any]]
            Overlapping regions from the interval index.
        """
        return self._index.overlap(chrom, pos, ref, alt)

    def assign(self, variant: Variant) -> FunctionalConsequence:
        """Assign the most severe functional consequence for a variant.

        Parameters
        ----------
        variant : Variant
            The variant to annotate.

        Returns
        -------
        FunctionalConsequence
            Most severe consequence, or INTERGENIC if no overlap.
        """
        overlaps = self._index.overlap(
            chrom=variant.chrom,
            pos=variant.pos,
            ref=variant.ref,
            alt=variant.alt,
        )

        if not overlaps:
            return FunctionalConsequence.INTERGENIC

        return _most_severe_consequence_pyranges(overlaps)

    def gene_names_batch(
        self, variants: list[Variant]
    ) -> list[Optional[str]]:
        """Extract gene names for a batch using a single vectorized join.

        Parameters
        ----------
        variants : list[Variant]
            Variants to look up gene names for.

        Returns
        -------
        list[Optional[str]]
            Gene names positionally matched to the batch.
        """
        if not variants:
            return []

        if not self._index._loaded or self._index._gr is None:
            return [None] * len(variants)

        records = []
        for i, v in enumerate(variants):
            var_start = v.pos - 1
            var_end = var_start + max(len(v.ref), len(v.alt))
            records.append({
                "Chromosome": v.chrom,
                "Start": var_start,
                "End": var_end,
                "_idx": i,
            })

        query_df = pd.DataFrame(records)
        query_gr = pr.PyRanges(query_df)

        hits = self._index._gr.join(query_gr)
        hits_df = hits.df

        gene_names: list[Optional[str]] = [None] * len(variants)

        if hits_df.empty:
            return gene_names

        # Pick the first hit per variant index
        gene_col = None
        for col in ("gene_name", "Gene", "gene_id"):
            if col in hits_df.columns:
                gene_col = col
                break

        if gene_col is None:
            return gene_names

        seen: set[int] = set()
        for _, row in hits_df.iterrows():
            var_idx = int(row["_idx"])
            if var_idx not in seen:
                seen.add(var_idx)
                val = row[gene_col]
                if pd.notna(val):
                    gene_names[var_idx] = str(val)

        return gene_names

    def assign_batch(self, variants: list[Variant]) -> list[FunctionalConsequence]:
        """Assign consequences to a batch of variants using vectorized join.

        Builds a single PyRanges from all variant positions, joins against
        the gene model in one operation, then maps consequences back.

        Parameters
        ----------
        variants : list[Variant]
            List of variants to annotate.

        Returns
        -------
        list[FunctionalConsequence]
            Consequences in the same order as input variants.
        """
        if not variants:
            return []

        if not self._index._loaded or self._index._gr is None:
            return [FunctionalConsequence.INTERGENIC] * len(variants)

        # Build a DataFrame of all variant positions at once
        records = []
        for i, v in enumerate(variants):
            var_start = v.pos - 1
            var_end = var_start + max(len(v.ref), len(v.alt))
            records.append({
                "Chromosome": v.chrom,
                "Start": var_start,
                "End": var_end,
                "_idx": i,
                "_ref": v.ref,
                "_alt": v.alt,
            })

        query_df = pd.DataFrame(records)
        query_gr = pr.PyRanges(query_df)

        # Single vectorized join against gene model
        hits = self._index._gr.join(query_gr)
        hits_df = hits.df

        # Initialize all as Intergenic
        results: list[FunctionalConsequence] = [
            FunctionalConsequence.INTERGENIC
        ] * len(variants)

        if hits_df.empty:
            return results

        # Build splice site lookup (vectorized)
        splice_positions = self._find_splice_positions(query_df)

        # Group hits by variant index, pick most severe consequence
        severity_rank = {
            c.value: idx for idx, c in enumerate(CONSEQUENCE_SEVERITY_ORDER)
        }

        for _, row in hits_df.iterrows():
            var_idx = int(row["_idx"])
            feature_type = row.get("Feature", "unknown")
            ref = row.get("_ref", "")
            alt = row.get("_alt", "")
            is_splice = var_idx in splice_positions

            consequence_str = _determine_consequence_pyranges(
                ref=ref, alt=alt,
                feature_type=feature_type,
                is_splice_site=is_splice,
            )

            new_rank = severity_rank.get(
                consequence_str, len(CONSEQUENCE_SEVERITY_ORDER)
            )
            current_rank = severity_rank.get(
                results[var_idx].value, len(CONSEQUENCE_SEVERITY_ORDER)
            )

            if new_rank < current_rank:
                results[var_idx] = FunctionalConsequence(consequence_str)

        return results

    def _find_splice_positions(
        self, query_df: pd.DataFrame,
    ) -> set[int]:
        """Identify variant indices that overlap splice sites."""
        splice_positions: set[int] = set()
        if self._index._exon_gr is None or self._index._exon_gr.df.empty:
            return splice_positions

        exon_df = self._index._exon_gr.df
        for _, v_row in query_df.iterrows():
            chrom = v_row["Chromosome"]
            chrom_exons = exon_df[exon_df["Chromosome"] == chrom]
            if chrom_exons.empty:
                continue

            vs, ve = v_row["Start"], v_row["End"]
            starts = chrom_exons["Start"].values
            ends = chrom_exons["End"].values

            donor_hit = ((vs < ends + 2) & (ve > ends - 2)).any()
            acceptor_hit = ((vs < starts + 2) & (ve > starts - 2)).any()
            if donor_hit or acceptor_hit:
                splice_positions.add(v_row["_idx"])

        return splice_positions


def _determine_consequence_pyranges(
    ref: str,
    alt: str,
    feature_type: str,
    is_splice_site: bool,
) -> str:
    """Determine functional consequence based on variant type and context.

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
    if is_splice_site:
        return FunctionalConsequence.SPLICE_SITE.value

    is_coding = feature_type == "CDS"

    if not is_coding:
        if feature_type in ("exon", "transcript", "gene"):
            return FunctionalConsequence.SYNONYMOUS.value
        return FunctionalConsequence.INTERGENIC.value

    ref_len = len(ref)
    alt_len = len(alt)
    is_snv = ref_len == 1 and alt_len == 1

    if is_snv:
        return FunctionalConsequence.MISSENSE.value

    length_diff = alt_len - ref_len

    if length_diff == 0:
        return FunctionalConsequence.MISSENSE.value

    if length_diff % 3 != 0:
        return FunctionalConsequence.FRAMESHIFT.value

    if length_diff > 0:
        return FunctionalConsequence.IN_FRAME_INSERTION.value
    return FunctionalConsequence.IN_FRAME_DELETION.value


def _most_severe_consequence_pyranges(
    overlaps: list[dict[str, Any]],
) -> FunctionalConsequence:
    """Select the most severe consequence from overlap results.

    Parameters
    ----------
    overlaps : list[dict[str, Any]]
        Overlap results containing 'consequence' keys.

    Returns
    -------
    FunctionalConsequence
        The most severe consequence found.
    """
    severity_rank = {c.value: idx for idx, c in enumerate(CONSEQUENCE_SEVERITY_ORDER)}

    best_consequence = FunctionalConsequence.INTERGENIC
    best_rank = severity_rank[best_consequence.value]

    for overlap in overlaps:
        consequence_str = overlap.get("consequence", FunctionalConsequence.INTERGENIC.value)
        rank = severity_rank.get(consequence_str, len(CONSEQUENCE_SEVERITY_ORDER))
        if rank < best_rank:
            best_rank = rank
            best_consequence = FunctionalConsequence(consequence_str)

    return best_consequence
