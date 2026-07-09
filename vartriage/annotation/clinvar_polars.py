"""Polars-based ClinVar clinical significance lookup.

Implements ClinVarDatabase using polars LazyFrames for batch left-join
lookups. Activated only when polars is installed (part of the
``accelerated`` optional extra).

The reference file format is TSV with columns:
    chrom, pos, ref, alt, clinical_significance

Clinical significance string values map to ClinVarAssertion enum members:
    "Pathogenic"           -> ClinVarAssertion.PATHOGENIC
    "Likely pathogenic"    -> ClinVarAssertion.LIKELY_PATHOGENIC
    "Uncertain significance" -> ClinVarAssertion.VUS
    "Likely benign"        -> ClinVarAssertion.LIKELY_BENIGN
    "Benign"               -> ClinVarAssertion.BENIGN
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

try:
    import polars as pl

    _POLARS_AVAILABLE = True
except ImportError:
    _POLARS_AVAILABLE = False

from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.variant import ClinVarAssertion


_SIGNIFICANCE_MAP: dict[str, ClinVarAssertion] = {
    "Pathogenic": ClinVarAssertion.PATHOGENIC,
    "Likely pathogenic": ClinVarAssertion.LIKELY_PATHOGENIC,
    "Uncertain significance": ClinVarAssertion.VUS,
    "Likely benign": ClinVarAssertion.LIKELY_BENIGN,
    "Benign": ClinVarAssertion.BENIGN,
}

_REVERSE_MAP: dict[ClinVarAssertion, str] = {
    v: k for k, v in _SIGNIFICANCE_MAP.items()
}


class PolarsClinVarDatabase:
    """Polars-based ClinVar lookup implementing ClinVarDatabase protocol.

    Uses polars LazyFrames for efficient batch left-join lookups against
    the ClinVar reference. This implementation is significantly faster
    than the dict-based fallback for large batch sizes, but requires
    polars to be installed.

    Parameters
    ----------
    None

    Raises
    ------
    ImportError
        If polars is not installed when the class is instantiated.

    Examples
    --------
    >>> db = PolarsClinVarDatabase()
    >>> db.load(Path("clinvar_reference.tsv"))
    >>> results = db.lookup_batch([("chr1", 12345, "A", "T")])
    """

    def __init__(self) -> None:
        if not _POLARS_AVAILABLE:
            raise ImportError(
                "polars is required for PolarsClinVarDatabase. "
                "Install with: pip install vartriage[accelerated]"
            )
        self._reference_df: Optional[pl.LazyFrame] = None
        self._loaded: bool = False

    def load(self, reference_path: Path) -> None:
        """Load ClinVar reference data from a TSV file into a LazyFrame.

        Parameters
        ----------
        reference_path : Path
            Path to the ClinVar reference file in TSV format with
            columns: chrom, pos, ref, alt, clinical_significance.

        Raises
        ------
        ReferenceFileError
            If the file does not exist, is not readable, or contains
            malformed data that cannot be parsed.
        """
        if not reference_path.exists():
            raise ReferenceFileError(
                f"{reference_path}: file not found"
            )

        if not reference_path.is_file():
            raise ReferenceFileError(
                f"{reference_path}: not a regular file"
            )

        try:
            df = pl.read_csv(
                reference_path,
                separator="\t",
                has_header=True,
                columns=[
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "clinical_significance",
                ],
                schema_overrides={
                    "chrom": pl.Utf8,
                    "pos": pl.Int64,
                    "ref": pl.Utf8,
                    "alt": pl.Utf8,
                    "clinical_significance": pl.Utf8,
                },
            )
        except Exception as exc:
            raise ReferenceFileError(
                f"{reference_path}: failed to parse ClinVar reference "
                f"with polars — {exc}"
            ) from exc

        # Filter to only rows with recognized significance values
        valid_significances = list(_SIGNIFICANCE_MAP.keys())
        df = df.filter(
            pl.col("clinical_significance").is_in(valid_significances)
        )

        self._reference_df = df.lazy()
        self._loaded = True

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[ClinVarAssertion]]:
        """Batch lookup of ClinVar assertions using polars left join.

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
        if not variants:
            return []

        if self._reference_df is None:
            return [None] * len(variants)

        # Build a query DataFrame from the input variants
        query_df = pl.DataFrame(
            {
                "chrom": [v[0] for v in variants],
                "pos": [v[1] for v in variants],
                "ref": [v[2] for v in variants],
                "alt": [v[3] for v in variants],
                "_idx": list(range(len(variants))),
            }
        ).lazy()

        # Left join against the reference
        joined = query_df.join(
            self._reference_df,
            on=["chrom", "pos", "ref", "alt"],
            how="left",
        ).sort("_idx")

        result_df = joined.collect()

        # Map significance strings back to enum values
        sig_column = result_df["clinical_significance"]
        results: list[Optional[ClinVarAssertion]] = []

        for value in sig_column:
            if value is None:
                results.append(None)
            else:
                results.append(_SIGNIFICANCE_MAP.get(value))

        return results
