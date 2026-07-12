"""Pure-Python dictionary-based gnomAD frequency lookup.

Dict-based FrequencyDatabase implementation. Loads the entire gnomAD
reference TSV into memory keyed on (chrom, pos, ref, alt) tuples.
Always available without optional dependencies.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.warnings import MissingDataWarning


class DictFrequencyDatabase:
    """Dictionary-based gnomAD frequency lookup.

    Loads a gnomAD reference TSV file into memory as a Python dict
    keyed on (chrom, pos, ref, alt) -> allele frequency. Suitable for
    datasets of moderate size; for whole-genome scale references, prefer
    the polars-based implementation.

    Parameters
    ----------
    None

    Attributes
    ----------
    warnings : list[MissingDataWarning]
        Accumulated warnings for variants not found in the database.

    Examples
    --------
    >>> db = DictFrequencyDatabase()
    >>> db.load(Path("gnomad_reference.tsv"))
    >>> results = db.lookup_batch([("chr1", 100, "A", "T")])
    """

    def __init__(self) -> None:
        self._data: dict[tuple[str, int, str, str], float] = {}
        self.warnings: list[MissingDataWarning] = []

    def load(self, reference_path: Path) -> None:
        """Load gnomAD reference data from a TSV file.

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
            with open(reference_path, "r", newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader, None)

                if header is None:
                    raise ReferenceFileError(f"{reference_path}: file is empty")

                expected_columns = {"chrom", "pos", "ref", "alt", "af"}
                header_lower = [col.lower().strip() for col in header]

                if not expected_columns.issubset(set(header_lower)):
                    missing = expected_columns - set(header_lower)
                    raise ReferenceFileError(
                        f"{reference_path}: missing required columns: "
                        f"{sorted(missing)}"
                    )

                col_idx = {name: header_lower.index(name) for name in expected_columns}

                for line_num, row in enumerate(reader, start=2):
                    if not row or all(field.strip() == "" for field in row):
                        continue

                    try:
                        chrom = row[col_idx["chrom"]].strip()
                        pos = int(row[col_idx["pos"]].strip())
                        ref = row[col_idx["ref"]].strip()
                        alt = row[col_idx["alt"]].strip()
                        af_str = row[col_idx["af"]].strip()

                        # '.' or empty means no frequency data, skip
                        if af_str in (".", ""):
                            continue

                        af = float(af_str)
                    except (IndexError, ValueError) as exc:
                        raise ReferenceFileError(
                            f"{reference_path}: parse error at "
                            f"line {line_num}: {exc}"
                        ) from exc

                    self._data[(chrom, pos, ref, alt)] = af

        except ReferenceFileError:
            raise
        except OSError as exc:
            raise ReferenceFileError(
                f"{reference_path}: cannot read file: {exc}"
            ) from exc

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[float]]:
        """Batch lookup of allele frequencies by genomic coordinate.

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
        results: list[Optional[float]] = []

        for chrom, pos, ref, alt in variants:
            freq = self._data.get((chrom, pos, ref, alt))
            if freq is None:
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
            results.append(freq)

        return results
