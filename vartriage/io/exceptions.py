"""Exception hierarchy for the variant prioritization library."""

from typing import Optional


class VariantPrioritizationError(Exception):
    """Base exception for all library errors.

    All domain-specific exceptions in this library inherit from
    this class, making it easy for callers to catch any library
    error with a single except clause while still allowing
    fine-grained handling when needed.
    """


class ParseError(VariantPrioritizationError):
    """Raised when VCF content cannot be parsed.

    Covers both header-level violations (missing ``##fileformat``
    declaration, malformed ``##INFO``/``##FORMAT`` meta-information
    lines, absent ``#CHROM`` column header) and data-line violations
    (missing mandatory columns, invalid field values).

    Parameters
    ----------
    line_number : int
        1-based line number in the source file where the error
        occurred.
    detail : str
        Human-readable description of the parse failure.
    field : str, optional
        Name of the specific VCF field that failed validation
        (e.g., ``"QUAL"``, ``"POS"``). ``None`` when the error
        is not tied to a single field (e.g., a missing column
        header line).

    Attributes
    ----------
    line_number : int
        1-based line number where the error occurred.
    field : Optional[str]
        Name of the field that failed validation, or ``None``.
    detail : str
        Human-readable description of the parse failure.

    Examples
    --------
    >>> raise ParseError(
    ...     line_number=5,
    ...     detail="Missing mandatory QUAL column",
    ... )
    Traceback (most recent call last):
        ...
    ParseError: Line 5: Missing mandatory QUAL column

    >>> raise ParseError(
    ...     line_number=12,
    ...     field="POS",
    ...     detail="Non-integer value 'abc'",
    ... )
    Traceback (most recent call last):
        ...
    ParseError: Line 12, field 'POS': Non-integer value 'abc'
    """

    def __init__(
        self,
        line_number: int,
        detail: str,
        field: Optional[str] = None,
    ) -> None:
        self.line_number: int = line_number
        self.field: Optional[str] = field
        self.detail: str = detail

        if field is not None:
            message = (
                f"Line {line_number}, field '{field}': {detail}"
            )
        else:
            message = f"Line {line_number}: {detail}"

        super().__init__(message)


class ConfigurationError(VariantPrioritizationError):
    """Raised when pipeline configuration is invalid.

    Thrown during configuration construction or pipeline
    initialization when a parameter value falls outside its
    permitted range or a required value is absent. The exception
    message names the invalid parameter and describes the valid
    range or format.

    Examples
    --------
    >>> raise ConfigurationError(
    ...     "min_qual must be between 0 and 1000000, got -5"
    ... )
    Traceback (most recent call last):
        ...
    ConfigurationError: min_qual must be between 0 and 1000000...
    """


class ReferenceFileError(VariantPrioritizationError):
    """Raised when a reference file cannot be loaded or parsed.

    Thrown when the annotation or prioritization engine is unable
    to open, read, or parse a reference data file (gnomAD, ClinVar,
    GTF/GFF gene annotation, CADD scores, REVEL scores). The
    exception message includes the file path and a description of
    the failure so the caller can distinguish a missing file from a
    format error.

    Examples
    --------
    >>> raise ReferenceFileError(
    ...     "/data/gnomad.vcf.gz: file not found"
    ... )
    Traceback (most recent call last):
        ...
    ReferenceFileError: /data/gnomad.vcf.gz: file not found
    """
