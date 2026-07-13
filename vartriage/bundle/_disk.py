"""Disk space utilities for pre-flight validation."""

from __future__ import annotations

import shutil
from pathlib import Path


def available_space_bytes(path: Path) -> int:
    """Return available disk space in bytes at the given path.

    If the path doesn't exist, walks up to the nearest existing
    parent to check available space.

    Parameters
    ----------
    path : Path
        Directory path to check (or a path whose parent will be used).

    Returns
    -------
    int
        Available bytes on the filesystem containing path.
    """
    check_path = path
    while not check_path.exists():
        check_path = check_path.parent
        if check_path == check_path.parent:
            break

    usage = shutil.disk_usage(check_path)
    return usage.free


def format_bytes(size_bytes: int) -> str:
    """Format a byte count as a human-readable string.

    Parameters
    ----------
    size_bytes : int
        Size in bytes.

    Returns
    -------
    str
        Formatted string (e.g., "4.7 GB", "80 MB", "1.2 KB").
    """
    if size_bytes < 0:
        return "0 B"

    units = [("TB", 1 << 40), ("GB", 1 << 30), ("MB", 1 << 20), ("KB", 1 << 10)]

    for suffix, threshold in units:
        if size_bytes >= threshold:
            value = size_bytes / threshold
            return f"{value:.1f} {suffix}"

    return f"{size_bytes} B"


def check_disk_space(path: Path, required_bytes: int) -> None:
    """Verify sufficient disk space is available.

    Parameters
    ----------
    path : Path
        Directory where files will be written.
    required_bytes : int
        Minimum bytes needed.

    Raises
    ------
    OSError
        If available space is less than required.
    """
    available = available_space_bytes(path)
    if available < required_bytes:
        raise OSError(
            f"Insufficient disk space at {path}:\n"
            f"  required:  {format_bytes(required_bytes)}\n"
            f"  available: {format_bytes(available)}"
        )
