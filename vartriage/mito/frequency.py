"""HelixMTdb mitochondrial population frequency lookup.

Provides allele frequency data for mtDNA variants from HelixMTdb,
a large-scale mitochondrial database. Used to distinguish common
haplogroup-defining polymorphisms from rare potentially pathogenic
variants.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class MtFrequencyEntry:
    """Population frequency record for a single mtDNA variant.

    Parameters
    ----------
    position
        1-based mtDNA position (rCRS).
    ref
        Reference allele.
    alt
        Alternate allele.
    af
        Allele frequency (0.0 to 1.0).
    allele_count
        Number of observed alleles in the dataset.
    """

    position: int
    ref: str
    alt: str
    af: float
    allele_count: int

    @property
    def is_common_haplogroup_marker(self) -> bool:
        """True if AF > 5%, indicating a haplogroup-defining polymorphism."""
        return self.af > 0.05

    @property
    def is_rare(self) -> bool:
        """True if AF < 0.01% (PM2 equivalent for mtDNA)."""
        return self.af < 0.0001


class HelixMTdbDatabase:
    """HelixMTdb population frequency lookup for mtDNA variants.

    Loads from the bundled TSV and provides O(1) dict-based lookups
    keyed on (position, ref, alt).

    Parameters
    ----------
    data_path
        Path to helixmtdb_frequency.tsv. If None, uses the
        package-bundled default.
    """

    def __init__(self, data_path: Path | None = None) -> None:
        self._entries: dict[tuple[int, str, str], MtFrequencyEntry] = {}
        path = data_path or self._default_path()
        self._load(path)

    def lookup(self, pos: int, ref: str, alt: str) -> MtFrequencyEntry | None:
        """Query population frequency for a mtDNA variant.

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
        MtFrequencyEntry or None
            Frequency data if found, None for novel variants.
        """
        return self._entries.get((pos, ref.upper(), alt.upper()))

    def get_af(self, pos: int, ref: str, alt: str) -> float | None:
        """Shortcut to get just the allele frequency value.

        Returns
        -------
        float or None
            Allele frequency, or None if the variant is not in the database.
        """
        entry = self.lookup(pos, ref, alt)
        return entry.af if entry is not None else None

    @property
    def size(self) -> int:
        """Number of entries loaded."""
        return len(self._entries)

    def _load(self, path: Path) -> None:
        """Parse the HelixMTdb TSV into the lookup dict."""
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    pos = int(row["position"])
                    af = float(row["af"])
                    allele_count = int(row["allele_count"])
                except (ValueError, KeyError):
                    continue

                entry = MtFrequencyEntry(
                    position=pos,
                    ref=row["ref"].upper(),
                    alt=row["alt"].upper(),
                    af=af,
                    allele_count=allele_count,
                )
                self._entries[(pos, entry.ref, entry.alt)] = entry

        logger.info(
            "Loaded %d HelixMTdb frequency entries from %s",
            len(self._entries),
            path,
        )

    @staticmethod
    def _default_path() -> Path:
        """Resolve the package-bundled helixmtdb_frequency.tsv."""
        data_dir = resources.files("vartriage") / "data" / "mito"
        return Path(str(data_dir / "helixmtdb_frequency.tsv"))
