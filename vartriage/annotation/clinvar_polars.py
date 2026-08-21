"""Polars-accelerated ClinVar clinical significance lookup.

Uses polars for fast TSV parsing, then converts to a Python dict for
O(1) per-variant lookups. The dict is cached to disk via pickle so
subsequent runs skip TSV parsing entirely (~0.03s load from cache vs
20+ seconds of parsing 4.4M rows).

Only available when polars is installed (part of the ``accelerated``
optional extra).

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

import logging
from pathlib import Path

try:
    import polars as pl

    _POLARS_AVAILABLE = True
except ImportError:
    _POLARS_AVAILABLE = False

from vartriage._internal.cache import try_load_cache, try_write_cache
from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.variant import ClinVarAssertion

logger = logging.getLogger(__name__)

_SIGNIFICANCE_MAP: dict[str, ClinVarAssertion] = {
    "Pathogenic": ClinVarAssertion.PATHOGENIC,
    "Likely pathogenic": ClinVarAssertion.LIKELY_PATHOGENIC,
    "Uncertain significance": ClinVarAssertion.VUS,
    "Likely benign": ClinVarAssertion.LIKELY_BENIGN,
    "Benign": ClinVarAssertion.BENIGN,
}

_REVERSE_MAP: dict[ClinVarAssertion, str] = {v: k for k, v in _SIGNIFICANCE_MAP.items()}


class PolarsClinVarDatabase:
    """ClinVar lookup using polars for parsing + dict for O(1) lookups.

    Parses ClinVar TSV via polars, converts to a Python dict keyed on
    (chrom, pos, ref, alt) -> ClinVarAssertion for O(1) per-variant
    lookups. The dict is pickle-cached so repeated runs load in ~0.03s.

    This replaces the previous approach of per-batch DataFrame.join()
    which was O(N*M) per batch (N=batch_size, M=reference_size).

    Parameters
    ----------
    None

    Raises
    ------
    ImportError
        If polars is not installed when the class is instantiated.
    """

    def __init__(self) -> None:
        if not _POLARS_AVAILABLE:
            raise ImportError(
                "polars is required for PolarsClinVarDatabase. "
                "Install with: pip install vartriage[accelerated]"
            )
        self._data: dict[tuple[str, int, str, str], ClinVarAssertion] = {}
        self._loaded: bool = False

    def load(self, reference_path: Path) -> None:
        """Load ClinVar reference data from a TSV file into a lookup dict.

        Checks for a pickle cache first. On cache hit, deserializes
        the dict directly. On cache miss, parses the TSV via polars,
        converts to dict, and writes the cache for next time.

        Parameters
        ----------
        reference_path : Path
            Path to the ClinVar reference file in TSV format.

        Raises
        ------
        ReferenceFileError
            If the file does not exist, is not readable, or contains
            malformed data that cannot be parsed.
        """
        if not reference_path.exists():
            raise ReferenceFileError(f"{reference_path}: file not found")

        if not reference_path.is_file():
            raise ReferenceFileError(f"{reference_path}: not a regular file")

        # Try loading from pickle cache
        cached = try_load_cache(reference_path)
        if cached is not None and isinstance(cached, dict):
            logger.info("ClinVar dict loaded from cache (%d entries)", len(cached))
            self._data = cached
            self._loaded = True
            return

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
                f"with polars: {exc}"
            ) from exc

        # Filter to only rows with recognized significance values
        valid_significances = list(_SIGNIFICANCE_MAP.keys())
        df = df.filter(pl.col("clinical_significance").is_in(valid_significances))

        # Convert to dict for O(1) lookups
        chroms = df["chrom"].to_list()
        positions = df["pos"].to_list()
        refs = df["ref"].to_list()
        alts = df["alt"].to_list()
        sigs = df["clinical_significance"].to_list()

        self._data = {
            (chroms[i], positions[i], refs[i], alts[i]): _SIGNIFICANCE_MAP[sigs[i]]
            for i in range(len(chroms))
        }

        self._loaded = True
        logger.info("ClinVar dict built: %d entries", len(self._data))

        # Write cache for next run
        try_write_cache(reference_path, self._data)

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[ClinVarAssertion | None]:
        """Batch lookup of ClinVar assertions via O(1) dict access.

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

        return [self._data.get(key) for key in variants]
