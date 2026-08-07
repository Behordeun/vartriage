"""MITOMAP pathogenic mutation database lookup.

Loads curated confirmed/reported pathogenic mtDNA mutations from the
shipped TSV and provides (position, ref, alt) keyed lookups. Data
sourced from MITOMAP (https://www.mitomap.org/) disease table.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class MitomapEntry:
    """A single MITOMAP pathogenic mutation record.

    Parameters
    ----------
    position
        1-based mtDNA position (rCRS).
    ref
        Reference allele.
    alt
        Alternate (pathogenic) allele.
    disease
        Associated disease phenotype(s).
    status
        Confirmation status: "Cfrm" (confirmed), "Reported", or "P.M."
    locus
        Mitochondrial gene name (e.g., "MT-TL1", "MT-ND1").
    """

    position: int
    ref: str
    alt: str
    disease: str
    status: str
    locus: str

    @property
    def is_confirmed(self) -> bool:
        """True if this mutation has confirmed pathogenicity in MITOMAP."""
        return self.status == "Cfrm"


class MitomapDatabase:
    """Lookup index for MITOMAP confirmed pathogenic mutations.

    Loads from the bundled TSV and provides O(1) dict-based lookups
    keyed on (position, ref, alt).

    Parameters
    ----------
    data_path
        Path to mitomap_pathogenic.tsv. If None, uses the
        package-bundled default.
    """

    def __init__(self, data_path: Path | None = None) -> None:
        self._entries: dict[tuple[int, str, str], MitomapEntry] = {}
        path = data_path or self._default_path()
        self._load(path)

    def lookup(self, pos: int, ref: str, alt: str) -> MitomapEntry | None:
        """Query for a known pathogenic mutation at a given position.

        Parameters
        ----------
        pos
            1-based mtDNA position.
        ref
            Reference allele (uppercase).
        alt
            Alternate allele (uppercase).

        Returns
        -------
        MitomapEntry or None
            The matching entry if found, None otherwise.
        """
        return self._entries.get((pos, ref.upper(), alt.upper()))

    @property
    def size(self) -> int:
        """Number of entries loaded."""
        return len(self._entries)

    def _load(self, path: Path) -> None:
        """Parse the MITOMAP TSV into the lookup dict."""
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    pos = int(row["position"])
                except (ValueError, KeyError):
                    continue

                entry = MitomapEntry(
                    position=pos,
                    ref=row["ref"].upper(),
                    alt=row["alt"].upper(),
                    disease=row["disease"],
                    status=row["status"],
                    locus=row["locus"],
                )
                self._entries[(pos, entry.ref, entry.alt)] = entry

        logger.info(
            "Loaded %d MITOMAP pathogenic entries from %s",
            len(self._entries),
            path,
        )

    @staticmethod
    def _default_path() -> Path:
        """Resolve the package-bundled mitomap_pathogenic.tsv."""
        data_dir = resources.files("vartriage") / "data" / "mito"
        return Path(str(data_dir / "mitomap_pathogenic.tsv"))
