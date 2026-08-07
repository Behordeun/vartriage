"""Mitochondrial gene map for position-based annotation.

Loads the 37 mitochondrial genes (+ D-loop control region) from the
shipped TSV and provides position-based lookups. Uses a sorted interval
list with binary search for fast queries against the 16,569bp genome.
"""

from __future__ import annotations

import bisect
import csv
import logging
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Literal

logger = logging.getLogger(__name__)

MtGeneType = Literal["protein_coding", "tRNA", "rRNA", "control_region", "intergenic"]


@dataclass(frozen=True, slots=True)
class MtGeneEntry:
    """A single mitochondrial gene interval.

    Parameters
    ----------
    gene_name
        Gene symbol (e.g., "MT-ND1", "MT-TL1", "MT-Dloop").
    start
        1-based start position on the rCRS reference.
    end
        1-based end position (inclusive).
    gene_type
        Functional category: protein_coding, tRNA, rRNA, or control_region.
    strand
        Transcription strand ("+" or "-").
    """

    gene_name: str
    start: int
    end: int
    gene_type: MtGeneType
    strand: str


@dataclass(frozen=True, slots=True)
class MtGeneContext:
    """Result of querying a position against the MT gene map.

    Parameters
    ----------
    gene_name
        Gene symbol at this position, or None if intergenic.
    gene_type
        Functional category at this position.
    is_in_coding_or_trna
        True if the position falls within a protein-coding gene or tRNA.
        These regions have higher pathogenicity prior for novel variants.
    """

    gene_name: str | None
    gene_type: MtGeneType
    is_in_coding_or_trna: bool


class MtGeneMap:
    """Mitochondrial gene interval lookup.

    Loads gene intervals from the shipped TSV and provides O(log n)
    position queries via binary search on sorted start coordinates.

    Parameters
    ----------
    data_path
        Path to the mt_gene_map.tsv file. If None, loads the
        package-bundled default.
    """

    def __init__(self, data_path: Path | None = None) -> None:
        self._entries: list[MtGeneEntry] = []
        self._starts: list[int] = []

        path = data_path or self._default_path()
        self._load(path)

    def query(self, position: int) -> MtGeneContext:
        """Look up which gene (if any) overlaps a given MT position.

        Parameters
        ----------
        position
            1-based position on the mitochondrial genome (1-16569).

        Returns
        -------
        MtGeneContext
            Gene context at the queried position. Returns intergenic
            if no gene interval contains the position.
        """
        # Binary search: find rightmost entry whose start <= position
        idx = bisect.bisect_right(self._starts, position) - 1

        # Check entries at and around the found index (overlapping genes exist)
        for check_idx in range(max(0, idx - 1), min(len(self._entries), idx + 3)):
            entry = self._entries[check_idx]
            if entry.start <= position <= entry.end:
                return MtGeneContext(
                    gene_name=entry.gene_name,
                    gene_type=entry.gene_type,
                    is_in_coding_or_trna=(
                        entry.gene_type in ("protein_coding", "tRNA")
                    ),
                )

        return MtGeneContext(
            gene_name=None,
            gene_type="intergenic",
            is_in_coding_or_trna=False,
        )

    @property
    def entries(self) -> list[MtGeneEntry]:
        """All loaded gene entries (sorted by start position)."""
        return self._entries

    def _load(self, path: Path) -> None:
        """Parse the TSV gene map into sorted interval entries."""
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene_type = row["type"]
                if gene_type not in (
                    "protein_coding",
                    "tRNA",
                    "rRNA",
                    "control_region",
                ):
                    logger.warning(
                        "Unknown gene type '%s' for %s, skipping",
                        gene_type,
                        row["gene_name"],
                    )
                    continue

                self._entries.append(
                    MtGeneEntry(
                        gene_name=row["gene_name"],
                        start=int(row["start"]),
                        end=int(row["end"]),
                        gene_type=gene_type,
                        strand=row["strand"],
                    )
                )

        # Sort by start for binary search
        self._entries.sort(key=lambda e: e.start)
        self._starts = [e.start for e in self._entries]

        logger.info("Loaded %d MT gene intervals from %s", len(self._entries), path)

    @staticmethod
    def _default_path() -> Path:
        """Resolve the package-bundled mt_gene_map.tsv path."""
        data_dir = resources.files("vartriage") / "data" / "mito"
        return Path(str(data_dir / "mt_gene_map.tsv"))
