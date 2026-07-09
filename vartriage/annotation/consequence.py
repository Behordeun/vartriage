"""Functional consequence assignment using pure-Python interval tree.

Assigns the most severe FunctionalConsequence to each variant by checking
coordinate overlaps against gene annotation data. Uses SortedArrayIntervalIndex
as the default backend (always available without optional dependencies).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from vartriage._internal.interval_tree import SortedArrayIntervalIndex
from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.variant import (
    CONSEQUENCE_SEVERITY_ORDER,
    FunctionalConsequence,
    Variant,
)


class ConsequenceAnnotator:
    """Assign functional consequence to variants using sorted interval tree.

    Loads gene annotation from a GTF/GFF file into a pure-Python sorted
    interval tree and determines the most severe consequence for each
    variant based on overlapping transcript features.

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

    Examples
    --------
    >>> from pathlib import Path
    >>> annotator = ConsequenceAnnotator(Path("gencode.gtf"))
    >>> consequence = annotator.assign(variant)
    """

    def __init__(self, annotation_path: Path) -> None:
        self._index = SortedArrayIntervalIndex()
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

        Queries the interval index for overlapping gene regions and returns
        the most severe consequence according to the severity ranking. When
        multiple transcripts overlap with different consequences, the most
        severe is selected. When no transcripts overlap, returns INTERGENIC.

        Parameters
        ----------
        variant : Variant
            The variant to annotate.

        Returns
        -------
        FunctionalConsequence
            The most severe consequence from overlapping transcripts,
            or INTERGENIC if no overlap exists.
        """
        overlaps = self._index.overlap(
            chrom=variant.chrom,
            pos=variant.pos,
            ref=variant.ref,
            alt=variant.alt,
        )

        if not overlaps:
            return FunctionalConsequence.INTERGENIC

        return _most_severe_consequence(overlaps)

    def assign_batch(self, variants: list[Variant]) -> list[FunctionalConsequence]:
        """Assign consequences to a batch of variants.

        Parameters
        ----------
        variants : list[Variant]
            List of variants to annotate.

        Returns
        -------
        list[FunctionalConsequence]
            Consequences in the same order as input variants.
        """
        return [self.assign(v) for v in variants]


def _most_severe_consequence(
    overlaps: list[dict[str, Any]],
) -> FunctionalConsequence:
    """Select the most severe consequence from a list of overlapping regions.

    Parameters
    ----------
    overlaps : list[dict[str, Any]]
        List of overlap results from the interval index, each containing
        a 'consequence' key with the consequence value string.

    Returns
    -------
    FunctionalConsequence
        The most severe consequence found.
    """
    # Build severity lookup: lower index = more severe
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
