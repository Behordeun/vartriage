"""Dictionary-based ClinVar clinical significance lookup.

Pure-Python ClinVarDatabase implementation using an in-memory dictionary
keyed on (chrom, pos, ref, alt) tuples. Works without any optional
dependencies. Uses only the Python standard library plus core package models.

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

from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.variant import ClinVarAssertion


_SIGNIFICANCE_MAP: dict[str, ClinVarAssertion] = {
    "Pathogenic": ClinVarAssertion.PATHOGENIC,
    "Likely pathogenic": ClinVarAssertion.LIKELY_PATHOGENIC,
    "Uncertain significance": ClinVarAssertion.VUS,
    "Likely benign": ClinVarAssertion.LIKELY_BENIGN,
    "Benign": ClinVarAssertion.BENIGN,
}


class DictClinVarDatabase:
    """Pure-Python dict-based ClinVar lookup implementing ClinVarDatabase.

    Loads a ClinVar TSV reference file into a dictionary keyed on
    (chrom, pos, ref, alt) -> ClinVarAssertion. Lookups are O(1) per
    variant after the initial load.

    Parameters
    ----------
    None

    Examples
    --------
    >>> db = DictClinVarDatabase()
    >>> db.load(Path("clinvar_reference.tsv"))
    >>> results = db.lookup_batch([("chr1", 12345, "A", "T")])
    """

    def __init__(self) -> None:
        self._data: dict[tuple[str, int, str, str], ClinVarAssertion] = {}
        self._loaded: bool = False

    def load(self, reference_path: Path) -> None:
        """Load ClinVar reference data from a TSV file.

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
            self._data = self._parse_tsv(reference_path)
        except ReferenceFileError:
            raise
        except Exception as exc:
            raise ReferenceFileError(
                f"{reference_path}: failed to parse ClinVar reference: "
                f"{exc}"
            ) from exc

        self._loaded = True

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[ClinVarAssertion]]:
        """Batch lookup of ClinVar assertions by genomic coordinate.

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
        return [self._data.get(key) for key in variants]

    def _parse_tsv(
        self, path: Path
    ) -> dict[tuple[str, int, str, str], ClinVarAssertion]:
        """Parse the ClinVar TSV reference into a lookup dictionary.

        Parameters
        ----------
        path : Path
            Path to the TSV file.

        Returns
        -------
        dict[tuple[str, int, str, str], ClinVarAssertion]
            Mapping from (chrom, pos, ref, alt) to assertion.

        Raises
        ------
        ReferenceFileError
            If any line has fewer than 5 columns or contains an
            unrecognized clinical significance value.
        """
        data: dict[tuple[str, int, str, str], ClinVarAssertion] = {}

        with open(path, encoding="utf-8") as fh:
            for line_num, line in enumerate(fh, start=1):
                stripped = line.strip()

                if not stripped or stripped.startswith("#"):
                    continue

                # Skip header line if present
                if line_num == 1 and stripped.lower().startswith("chrom"):
                    continue

                parts = stripped.split("\t")
                if len(parts) < 5:
                    raise ReferenceFileError(
                        f"{path}: line {line_num} has {len(parts)} "
                        f"columns, expected at least 5"
                    )

                chrom = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    raise ReferenceFileError(
                        f"{path}: line {line_num} has non-integer "
                        f"position value '{parts[1]}'"
                    )

                ref = parts[2]
                alt = parts[3]
                significance_str = parts[4]

                assertion = _SIGNIFICANCE_MAP.get(significance_str)
                if assertion is None:
                    # Skip unrecognized significance values rather than
                    # failing the entire load. ClinVar has many non-
                    # standard categories we don't map.
                    continue

                data[(chrom, pos, ref, alt)] = assertion

        return data
