"""Warning models for missing data and validation events."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True, slots=True)
class MissingDataWarning:
    """Emitted when a reference database query returns no result for a variant.

    Captures the full variant identity and the reference source that failed
    to produce data, allowing downstream consumers to audit which variants
    are missing annotations and from which databases.

    Parameters
    ----------
    chrom : str
        Chromosome name of the variant (e.g., ``"chr1"``, ``"X"``).
    pos : int
        1-based genomic position of the variant.
    ref : str
        Reference allele string.
    alt : str
        Alternate allele string.
    source : str
        Name of the reference source that returned no data (e.g.,
        ``"gnomAD"``, ``"ClinVar"``).
    reason : str, optional
        Additional context about why data is missing. Distinguishes a normal
        "not found" result from infrastructure failures such as
        ``"connection_timeout"`` or ``"parse_error"``.  Defaults to ``None``
        for a standard absence (variant simply not present in the database).

    Examples
    --------
    >>> warning = MissingDataWarning(
    ...     chrom="chr1",
    ...     pos=12345,
    ...     ref="A",
    ...     alt="T",
    ...     source="gnomAD",
    ... )
    >>> warning.source
    'gnomAD'
    >>> warning.reason is None
    True
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    source: str
    reason: Optional[str] = None
