"""Abstract Protocol interfaces for pluggable pipeline backends.

Contracts that both pure-Python fallback and optimized backends must satisfy.
The pipeline depends only on these protocols — never on concrete library types.
Each protocol has two implementations:

- A pure-Python fallback (sorted-array interval tree, dictionary-based lookups,
  no PDF) that works with only pysam + numpy installed.
- An optimized backend (pyranges for intervals, polars for batch joins,
  reportlab for PDF) activated when optional extras are installed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Protocol

from vartriage.models.variant import ClinVarAssertion


class IntervalIndex(Protocol):
    """Interface for genomic interval overlap queries.

    Implementations must support loading gene annotation data from GTF/GFF
    files and performing coordinate-based overlap lookups for variant
    consequence assignment.

    Implementations
    ---------------
    SortedArrayIntervalIndex : pure-Python, always available
    PyRangesIntervalIndex : requires pyranges extra

    Methods
    -------
    load(annotation_path)
        Load gene annotation from GTF/GFF file.
    overlap(chrom, pos, ref, alt)
        Return overlapping gene regions for a variant coordinate.
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


class FrequencyDatabase(Protocol):
    """Interface for population frequency lookups.

    Implementations must support loading reference frequency data and
    performing batch lookups by genomic coordinate.

    Implementations
    ---------------
    DictFrequencyDatabase : pure-Python dict, always available
    PolarsFrequencyDatabase : requires polars extra

    Methods
    -------
    load(reference_path)
        Load reference frequency data from file.
    lookup_batch(variants)
        Batch lookup of allele frequencies by coordinate.
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
    """Interface for ClinVar clinical significance lookups.

    Implementations must support loading ClinVar reference data and
    performing batch lookups by genomic coordinate.

    Implementations
    ---------------
    DictClinVarDatabase : pure-Python dict, always available
    PolarsClinVarDatabase : requires polars extra

    Methods
    -------
    load(reference_path)
        Load ClinVar reference data from file.
    lookup_batch(variants)
        Batch lookup of ClinVar assertions by coordinate.
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
    """Interface for PDF report generation.

    Implementations must support rendering classified variant data into
    a formatted clinical PDF report.

    Implementations
    ---------------
    ReportlabPDFRenderer : requires reportlab extra
    PDFFallback : raises ImportError with install instructions

    Methods
    -------
    render(variants, output_path)
        Render classified variants to a PDF clinical report.
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
