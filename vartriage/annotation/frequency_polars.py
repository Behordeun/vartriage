"""Polars-accelerated gnomAD frequency lookup.

Uses polars for fast TSV parsing, then converts to a Python dict for
O(1) per-variant lookups. The dict is cached to disk via pickle so
subsequent runs skip TSV parsing entirely (~0.3s load from cache vs
30+ seconds of parsing 4.8M rows).

Only available when polars is installed (``pip install vartriage[accelerated]``).
"""

from __future__ import annotations

import logging
from pathlib import Path

try:
    import polars as pl

    POLARS_AVAILABLE = True
except ImportError:
    POLARS_AVAILABLE = False

from vartriage._internal.cache import try_load_cache, try_write_cache
from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.resolution import Resolution
from vartriage.models.warnings import MissingDataWarning

logger = logging.getLogger(__name__)


class PolarsFrequencyDatabase:
    """gnomAD frequency lookup using polars for parsing + dict for lookups.

    Parses gnomAD TSV via polars (fast columnar read), converts to a
    Python dict keyed on (chrom, pos, ref, alt) for O(1) per-variant
    lookups. The dict is pickle-cached so repeated runs load in ~0.3s.

    This replaces the previous approach of per-batch DataFrame.join()
    which was O(N*M) per batch (N=batch_size, M=reference_size).

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
    """

    def __init__(self) -> None:
        if not POLARS_AVAILABLE:
            raise ImportError(
                "polars is required for PolarsFrequencyDatabase. "
                "Install with: pip install "
                "vartriage[accelerated]"
            )
        self._data: dict[tuple[str, int, str, str], float] = {}
        self._loaded_chroms: set[str] = set()
        self._loaded: bool = False
        self.warnings: list[MissingDataWarning] = []

    def load(self, reference_path: Path) -> None:
        """Load gnomAD reference data from a TSV file into a lookup dict.

        Checks for a pickle cache first. On cache hit, deserializes
        the dict directly. On cache miss, parses the TSV via polars,
        converts to dict, and writes the cache for next time.

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
            an invalid format.
        """
        if not reference_path.exists():
            raise ReferenceFileError(f"{reference_path}: file not found")

        if not reference_path.is_file():
            raise ReferenceFileError(f"{reference_path}: not a regular file")

        # Try loading from pickle cache
        cached = try_load_cache(reference_path)
        if cached is not None and isinstance(cached, dict):
            logger.info(
                "gnomAD frequency dict loaded from cache (%d entries)", len(cached)
            )
            self._data = cached
            self._loaded_chroms = {k[0] for k in self._data}
            self._loaded = True
            return

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
                f"{reference_path}: failed to parse with polars: {exc}"
            ) from exc

        column_names = {col.lower() for col in df.columns}
        expected_columns = {"chrom", "pos", "ref", "alt", "af"}

        if not expected_columns.issubset(column_names):
            missing = expected_columns - column_names
            raise ReferenceFileError(
                f"{reference_path}: missing required columns: {sorted(missing)}"
            )

        # Normalize column names to lowercase
        df = df.rename({col: col.lower() for col in df.columns})

        # Convert to dict for O(1) lookups
        df = df.select(["chrom", "pos", "ref", "alt", "af"]).drop_nulls(subset=["af"])
        chroms = df["chrom"].to_list()
        positions = df["pos"].to_list()
        refs = df["ref"].to_list()
        alts = df["alt"].to_list()
        afs = df["af"].to_list()

        self._data = {
            (chroms[i], positions[i], refs[i], alts[i]): afs[i]
            for i in range(len(chroms))
        }
        self._loaded_chroms = set(chroms)

        logger.info("gnomAD frequency dict built: %d entries", len(self._data))

        self._loaded = True

        # Write cache for next run
        try_write_cache(reference_path, self._data)

    def _is_chrom_represented(self, chrom: str) -> bool:
        """Check whether a chromosome is represented in the loaded reference."""
        if chrom in self._loaded_chroms:
            return True
        if chrom.startswith("chr") and chrom[3:] in self._loaded_chroms:
            return True
        return not chrom.startswith("chr") and f"chr{chrom}" in self._loaded_chroms

    def lookup_batch_resolutions(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Resolution[float]]:
        """Batch lookup of allele frequencies with tri-state resolution outcomes.

        For each variant tuple:
        - If the database is not loaded, returns Resolution.lookup_failed.
        - If the chromosome is not represented in the loaded reference, returns
          Resolution.lookup_failed.
        - If the chromosome is represented and the variant is found, returns
          Resolution.found(af).
        - If the chromosome is represented and the variant is absent, returns
          Resolution.confirmed_absent().
        """
        if not variants:
            return []

        if not self._loaded:
            return [
                Resolution.lookup_failed("No reference data loaded") for _ in variants
            ]

        resolutions: list[Resolution[float]] = []
        for chrom, pos, ref, alt in variants:
            if not self._is_chrom_represented(chrom):
                resolutions.append(
                    Resolution.lookup_failed(
                        f"Chromosome '{chrom}' not present in reference database"
                    )
                )
                continue

            freq = self._data.get((chrom, pos, ref, alt))
            if freq is None:
                alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                freq = self._data.get((alt_chrom, pos, ref, alt))

            if freq is not None:
                resolutions.append(Resolution.found(freq))
            else:
                resolutions.append(Resolution.confirmed_absent())

        return resolutions

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[float | None]:
        """Batch lookup of allele frequencies via O(1) dict access.

        For each variant tuple not found or unrepresented in the loaded reference,
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
        if not variants:
            return []

        if not self._loaded:
            raise ReferenceFileError(
                "No reference data loaded. Call load() before lookup_batch()."
            )

        resolutions = self.lookup_batch_resolutions(variants)
        results: list[float | None] = []

        for (chrom, pos, ref, alt), res in zip(variants, resolutions, strict=True):
            if res.is_found:
                results.append(res.value)
            else:
                self.warnings.append(
                    MissingDataWarning(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        source="gnomAD",
                        reason="not_found"
                        if res.is_confirmed_absent
                        else "lookup_failed",
                    )
                )
                results.append(None)

        return results
