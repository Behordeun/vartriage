"""Header validation shared by the knowledge-library TSV loaders.

A loader that reads by column name loads nothing when a column is renamed
upstream, and the only trace is an info log. Validating the header up front
turns that silent zero-load into a clear error naming the file and the
missing columns.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path


def require_columns(
    reader: csv.DictReader[str],
    expected: Iterable[str],
    source: Path,
) -> None:
    """Raise if the reader's header is missing any expected column.

    Parameters
    ----------
    reader : csv.DictReader
        A reader whose ``fieldnames`` reflect the file header.
    expected : Iterable[str]
        Column names the loader requires.
    source : Path
        The file being read, named in the error message.

    Raises
    ------
    ValueError
        If the header is absent or any expected column is missing.
    """
    fieldnames = set(reader.fieldnames or [])
    missing = [c for c in expected if c not in fieldnames]
    if missing:
        raise ValueError(
            f"{source}: knowledge TSV is missing required column(s): "
            f"{', '.join(missing)}. Found: "
            f"{', '.join(sorted(fieldnames)) or '(no header)'}"
        )
