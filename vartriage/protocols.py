"""Protocol interfaces for pluggable pipeline backends.

The pipeline depends only on these protocols, not on concrete library types.
Each protocol has two implementations:

- Pure-Python fallback (sorted-array intervals, dict lookups, no PDF):
  works with just pysam + numpy.
- Optimized backend (pyranges, polars, reportlab): used when the
  optional extras are installed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Protocol

from vartriage.models.variant import (ClinVarAssertion, FunctionalConsequence,
                                      Variant)


class IntervalIndex(Protocol):
    """Genomic interval overlap queries.

    Implementations
    ---------------
    SortedArrayIntervalIndex : pure-Python, always available
    PyRangesIntervalIndex : requires pyranges extra
    """

    def load(self, annotation_path: Path) -> None:
        """Load gene annotation from a GTF/GFF file.

        Parameters
        ----------
        annotation_path : Path
            Path to the GTF or GFF gene annotation file.

        Raises
        ------
        FileNotFoundError
            If the annotation file does not exist.
        ValueError
            If the file cannot be parsed as valid GTF/GFF.
        """
        ...

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
            List of overlapping regions, each represented as a dictionary
            containing at minimum 'gene_name', 'feature_type', and
            'consequence' keys. Empty list when no overlaps are found.
        """
        ...

    def assign_batch(self, variants: list["Variant"]) -> list["FunctionalConsequence"]:
        """Assign consequences to a batch of variants.

        Returns a list of the same length as the input, positionally
        matched.

        Parameters
        ----------
        variants : list[Variant]
            Variants to annotate.

        Returns
        -------
        list[FunctionalConsequence]
            Consequences, same order as input.
        """
        ...


class FrequencyDatabase(Protocol):
    """Population allele frequency lookups.

    Implementations
    ---------------
    DictFrequencyDatabase : pure-Python dict, always available
    PolarsFrequencyDatabase : requires polars extra
    """

    def load(self, reference_path: Path) -> None:
        """Load reference frequency data from file.

        Parameters
        ----------
        reference_path : Path
            Path to the gnomAD or equivalent frequency reference file.

        Raises
        ------
        FileNotFoundError
            If the reference file does not exist.
        ValueError
            If the file cannot be parsed.
        """
        ...

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[float]]:
        """Batch lookup of allele frequencies by genomic coordinate.

        Parameters
        ----------
        variants : list[tuple[str, int, str, str]]
            List of (chrom, pos, ref, alt) tuples to look up.

        Returns
        -------
        list[Optional[float]]
            Allele frequencies in the same order as input. None for
            variants not found in the reference database.
        """
        ...


class ClinVarDatabase(Protocol):
    """ClinVar clinical significance lookups.

    Implementations
    ---------------
    DictClinVarDatabase : pure-Python dict, always available
    PolarsClinVarDatabase : requires polars extra
    """

    def load(self, reference_path: Path) -> None:
        """Load ClinVar reference data from file.

        Parameters
        ----------
        reference_path : Path
            Path to the ClinVar reference file.

        Raises
        ------
        FileNotFoundError
            If the reference file does not exist.
        ValueError
            If the file cannot be parsed.
        """
        ...

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[ClinVarAssertion]]:
        """Batch lookup of ClinVar clinical significance assertions.

        Parameters
        ----------
        variants : list[tuple[str, int, str, str]]
            List of (chrom, pos, ref, alt) tuples to look up.

        Returns
        -------
        list[Optional[ClinVarAssertion]]
            ClinVar assertions in the same order as input. None for
            variants not found in the ClinVar database.
        """
        ...


class PDFRenderer(Protocol):
    """PDF report generation.

    Implementations
    ---------------
    ReportlabPDFRenderer : requires reportlab extra
    PDFFallback : raises ImportError with install instructions
    """

    def render(self, variants: list[Any], output_path: Path) -> Path:
        """Render classified variants to a PDF clinical report.

        Parameters
        ----------
        variants : list[Any]
            List of ClassifiedVariant instances to include in the report.
        output_path : Path
            Filesystem path where the PDF should be written.

        Returns
        -------
        Path
            The path to the generated PDF file.

        Raises
        ------
        ImportError
            If the required PDF rendering backend is not installed.
        IOError
            If the file cannot be written to the specified path.
        """
        ...
