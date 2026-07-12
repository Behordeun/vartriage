"""Score file loading for CADD, REVEL, and SpliceAI pathogenicity lookups.

Parses tab-separated score files into dicts keyed by (chrom, pos, ref, alt)
for O(1) per-variant lookups during batch scoring.

Expected TSV format:
    - Columns: chrom, pos, ref, alt, score (tab-separated)
    - Lines starting with ``#`` are skipped
    - Malformed lines are skipped with a logged warning
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

CoordinateKey = tuple[str, int, str, str]


class ScoreLoader:
    """Loads CADD/REVEL/SpliceAI TSV files into coordinate-keyed dicts.

    All three formats share the same column layout so parsing is unified
    internally. Batch lookups return scores (or None) by coordinate.
    """

    def load_cadd(self, path: Path) -> dict[CoordinateKey, float]:
        """Parse a CADD TSV into a coordinate→score dict.

        Parameters
        ----------
        path : Path
            Tab-separated CADD score file.

        Returns
        -------
        dict[CoordinateKey, float]
            (chrom, pos, ref, alt) → CADD Phred score.

        Raises
        ------
        ValueError
            If the file doesn't exist or isn't readable.
        """
        return self._load_tsv(path)

    def load_revel(self, path: Path) -> dict[CoordinateKey, float]:
        """Parse a REVEL TSV into a coordinate→score dict.

        Parameters
        ----------
        path : Path
            Tab-separated REVEL score file.

        Returns
        -------
        dict[CoordinateKey, float]
            (chrom, pos, ref, alt) → REVEL score.

        Raises
        ------
        ValueError
            If the file doesn't exist or isn't readable.
        """
        return self._load_tsv(path)

    def load_spliceai(self, path: Path) -> dict[CoordinateKey, float]:
        """Parse a SpliceAI TSV into a coordinate→score dict.

        Parameters
        ----------
        path : Path
            Tab-separated SpliceAI score file.

        Returns
        -------
        dict[CoordinateKey, float]
            (chrom, pos, ref, alt) → SpliceAI score.

        Raises
        ------
        ValueError
            If the file doesn't exist or isn't readable.
        """
        return self._load_tsv(path)

    def lookup_batch(
        self,
        variants: list[CoordinateKey],
        score_dict: dict[CoordinateKey, float],
    ) -> list[Optional[float]]:
        """Look up scores for a list of coordinates.

        Parameters
        ----------
        variants : list[CoordinateKey]
            (chrom, pos, ref, alt) tuples to query.
        score_dict : dict[CoordinateKey, float]
            Pre-loaded dict from ``load_cadd`` or ``load_revel``.

        Returns
        -------
        list[Optional[float]]
            Scores in the same order as input; None where not found.
        """
        return [score_dict.get(key) for key in variants]

    def _load_tsv(self, path: Path) -> dict[CoordinateKey, float]:
        """Parse a TSV score file, skipping comments and bad lines.

        Uses pickle-based caching to avoid re-parsing on subsequent
        runs. Falls through to fresh parsing when cache is absent
        or invalid.

        Parameters
        ----------
        path : Path
            TSV file to parse.

        Returns
        -------
        dict[CoordinateKey, float]
            Parsed coordinate→score mapping.

        Raises
        ------
        ValueError
            If the file doesn't exist or isn't readable.
        """
        from typing import cast

        from vartriage._internal.cache import try_load_cache, try_write_cache

        cached = try_load_cache(path)
        if cached is not None:
            return cast(dict[CoordinateKey, float], cached)

        self._validate_path(path)

        scores: dict[CoordinateKey, float] = {}

        with open(path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh, start=1):
                stripped = line.rstrip("\n\r")

                if stripped.startswith("#") or stripped == "":
                    continue

                parts = stripped.split("\t")
                if len(parts) < 5:
                    logger.warning(
                        "Skipping malformed line %d in %s: "
                        "expected 5 tab-separated columns, got %d",
                        lineno,
                        path,
                        len(parts),
                    )
                    continue

                chrom = parts[0]
                pos_str = parts[1]
                ref = parts[2]
                alt = parts[3]
                score_str = parts[4]

                try:
                    pos = int(pos_str)
                except ValueError:
                    logger.warning(
                        "Skipping malformed line %d in %s: "
                        "non-integer position value '%s'",
                        lineno,
                        path,
                        pos_str,
                    )
                    continue

                try:
                    score = float(score_str)
                except ValueError:
                    logger.warning(
                        "Skipping malformed line %d in %s: "
                        "non-numeric score value '%s'",
                        lineno,
                        path,
                        score_str,
                    )
                    continue

                key: CoordinateKey = (chrom, pos, ref, alt)
                scores[key] = score

        try_write_cache(path, scores)
        return scores

    def _validate_path(self, path: Path) -> None:
        """Raise ValueError if the file is missing or unreadable."""
        if not path.exists():
            raise ValueError(f"Score file not found: {path}")
        if not path.is_file():
            raise ValueError(f"Score file path is not a regular file: {path}")
        if not os.access(path, os.R_OK):
            raise ValueError(f"Score file not readable: {path}")
