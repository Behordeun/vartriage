"""Polars-based gnomAD frequency lookup (optional accelerated backend).

Implements FrequencyDatabase using polars LazyFrame batch left joins for
efficient lookups on large gnomAD reference datasets. Only available when
polars is installed (``pip install vartriage[accelerated]``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

try:
    import polars as pl

    POLARS_AVAILABLE = True
except ImportError:
    POLARS_AVAILABLE = False

from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.warnings import MissingDataWarning


class PolarsFrequencyDatabase:
    """Polars-based gnomAD frequency lookup using LazyFrame joins.

    Loads a gnomAD reference TSV file into a polars LazyFrame and
    performs batch left joins for frequency lookups. Significantly
    faster than the dict-based fallback for large datasets.

    Requires polars to be installed. Check ``POLARS_AVAILABLE``
    before instantiating.

    Parameters
    ----------
    None

    Attributes
    ----------
    warnings : list[MissingDataWarning]
        Accumulated warnings for variants not found in the database.

    Raises
    ------
    ImportError
        If polars is not installed.

    Examples
    --------
    >>> from vartriage.annotation.frequency_polars import (
    ...     PolarsFrequencyDatabase, POLARS_AVAILABLE,
    ... )
    >>> if POLARS_AVAILABLE:
    ...     db = PolarsFrequencyDatabase()
    ...     db.load(Path("gnomad_reference.tsv"))
    ...     results = db.lookup_batch([("chr1", 100, "A", "T")])
    """

    def __init__(self) -> None:
        if not POLARS_AVAILABLE:
            raise ImportError(
                "polars is required for PolarsFrequencyDatabase. "
                "Install with: pip install "
                "vartriage[accelerated]"
            )
        self._lazy_frame: Optional[pl.LazyFrame] = None
        self.warnings: list[MissingDataWarning] = []

    def load(self, reference_path: Path) -> None:
        """Load gnomAD reference data from a TSV file into a LazyFrame.

        The expected file format is tab-separated with columns:
        chrom, pos, ref, alt, af

        Parameters
        ----------
        reference_path : Path
            Path to the gnomAD reference TSV file.

        Raises
        ------
        ReferenceFileError
            If the file does not exist, cannot be read, or has
            an invalid format (missing required columns or
            unparseable values).
        """
        if not reference_path.exists():
            raise ReferenceFileError(f"{reference_path}: file not found")

        if not reference_path.is_file():
            raise ReferenceFileError(f"{reference_path}: not a regular file")

        try:
            df = pl.read_csv(
                reference_path,
                separator="\t",
                has_header=True,
                null_values=[".", ""],
                schema_overrides={
                    "chrom": pl.Utf8,
                    "pos": pl.Int64,
                    "ref": pl.Utf8,
                    "alt": pl.Utf8,
                    "af": pl.Float64,
                },
            )
        except Exception as exc:
            raise ReferenceFileError(
                f"{reference_path}: failed to parse with polars: " f"{exc}"
            ) from exc

        column_names = {col.lower() for col in df.columns}
        expected_columns = {"chrom", "pos", "ref", "alt", "af"}

        if not expected_columns.issubset(column_names):
            missing = expected_columns - column_names
            raise ReferenceFileError(
                f"{reference_path}: missing required columns: " f"{sorted(missing)}"
            )

        # Normalize column names to lowercase
        df = df.rename({col: col.lower() for col in df.columns})

        # Select only the columns we need and store as LazyFrame
        self._lazy_frame = df.select(["chrom", "pos", "ref", "alt", "af"]).lazy()

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[float]]:
        """Batch lookup of allele frequencies via polars left join.

        For each variant tuple not found in the loaded reference,
        a MissingDataWarning is appended to `self.warnings`.

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
        if self._lazy_frame is None:
            raise ReferenceFileError("No reference data loaded. Call load() first.")

        if not variants:
            return []

        # Build a query DataFrame from the input variants
        query_df = pl.DataFrame(
            {
                "chrom": [v[0] for v in variants],
                "pos": [v[1] for v in variants],
                "ref": [v[2] for v in variants],
                "alt": [v[3] for v in variants],
                "row_idx": list(range(len(variants))),
            }
        )

        # Left join against the reference LazyFrame
        result = (
            query_df.lazy()
            .join(
                self._lazy_frame,
                on=["chrom", "pos", "ref", "alt"],
                how="left",
            )
            .sort("row_idx")
            .collect()
        )

        af_series = result["af"]
        results: list[Optional[float]] = []

        for i, freq_val in enumerate(af_series):
            if freq_val is None:
                chrom, pos, ref, alt = variants[i]
                self.warnings.append(
                    MissingDataWarning(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        source="gnomAD",
                        reason="not_found",
                    )
                )
                results.append(None)
            else:
                results.append(float(freq_val))

        return results
