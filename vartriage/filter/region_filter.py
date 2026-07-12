"""BED-based genomic region filter for variant streaming."""

from __future__ import annotations

import bisect
from pathlib import Path
from typing import Iterator

from vartriage.io.exceptions import ParseError
from vartriage.models.config import RegionFilterConfig
from vartriage.models.variant import Variant


class RegionFilter:
    """Filter variants to those overlapping BED file intervals.

    Builds a sorted interval index keyed by chromosome for
    binary-search overlap queries.

    Parameters
    ----------
    config : RegionFilterConfig
        Must contain the BED file path.

    Raises
    ------
    FileNotFoundError
        BED file does not exist.
    ParseError
        BED file contains malformed lines.
    """

    def __init__(self, config: RegionFilterConfig) -> None:
        self._intervals: dict[str, list[tuple[int, int]]] = {}
        self._starts: dict[str, list[int]] = {}
        self._load_bed(config.bed_path)

    def apply(
        self, variants: Iterator[Variant]
    ) -> Iterator[Variant]:
        """Yield variants overlapping at least one BED interval.

        Parameters
        ----------
        variants : Iterator[Variant]
            Input variant stream.

        Yields
        ------
        Variant
            Variants whose 0-based position falls within a BED interval.
        """
        for variant in variants:
            if self._overlaps(variant.chrom, variant.pos):
                yield variant

    def _overlaps(self, chrom: str, pos: int) -> bool:
        """Check if a 1-based VCF position overlaps any interval.

        Converts pos to 0-based (pos - 1) and checks against
        [start, end) intervals using bisect for O(log n) lookup.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        pos : int
            1-based VCF position.

        Returns
        -------
        bool
        """
        intervals = self._intervals.get(chrom)
        if not intervals:
            return False

        starts = self._starts[chrom]
        query_pos = pos - 1

        # Find the rightmost interval whose start <= query_pos.
        # bisect_right gives the insertion point; index - 1 is the
        # last interval with start <= query_pos.
        idx = bisect.bisect_right(starts, query_pos) - 1

        if idx < 0:
            return False

        # Scan backward: check each candidate whose start <= query_pos
        # to see if its end extends past query_pos.
        for i in range(idx, -1, -1):
            _, end = intervals[i]
            if end > query_pos:
                return True
            # For typical non-overlapping BED files this loop
            # terminates after 1 iteration.

        return False

    def _load_bed(self, bed_path: Path) -> None:
        """Parse BED file and build sorted interval index.

        Skips comments (#) and browser/track header lines.
        Validates each data line has >= 3 columns with valid coords.

        Parameters
        ----------
        bed_path : Path
            Path to the BED file.

        Raises
        ------
        FileNotFoundError
            File does not exist.
        ParseError
            Any data line is malformed.
        """
        if not bed_path.exists():
            raise FileNotFoundError(
                f"BED file not found: {bed_path}"
            )

        with open(bed_path, "r") as fh:
            for line_num, raw_line in enumerate(fh, start=1):
                line = raw_line.strip()
                if self._is_skippable_line(line):
                    continue

                chrom, start, end = self._parse_bed_line(
                    line, line_num
                )
                if chrom not in self._intervals:
                    self._intervals[chrom] = []
                self._intervals[chrom].append((start, end))

        # Sort intervals by (start, end) for binary search
        for chrom in self._intervals:
            self._intervals[chrom].sort()
            self._starts[chrom] = [
                iv[0] for iv in self._intervals[chrom]
            ]

    @staticmethod
    def _is_skippable_line(line: str) -> bool:
        """Return True if the line should be skipped during BED parsing."""
        return (
            not line
            or line.startswith("#")
            or line.startswith(("browser", "track"))
        )

    @staticmethod
    def _parse_bed_line(
        line: str, line_num: int,
    ) -> tuple[str, int, int]:
        """Parse a single BED data line into (chrom, start, end)."""
        fields = line.split("\t")
        if len(fields) < 3:
            raise ParseError(
                line_number=line_num,
                detail=(
                    "Expected at least 3 tab-separated "
                    f"columns, got {len(fields)}: "
                    f"{line!r}"
                ),
            )

        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            raise ParseError(
                line_number=line_num,
                detail=(
                    "Start and end must be non-negative "
                    f"integers: {fields[1]!r}, "
                    f"{fields[2]!r}"
                ),
            )

        if start < 0 or end < 0:
            raise ParseError(
                line_number=line_num,
                detail=(
                    "Coordinates must be non-negative, "
                    f"got start={start}, end={end}"
                ),
            )

        if start >= end:
            raise ParseError(
                line_number=line_num,
                detail=(
                    "Start must be less than end, "
                    f"got start={start}, end={end}"
                ),
            )

        return chrom, start, end
