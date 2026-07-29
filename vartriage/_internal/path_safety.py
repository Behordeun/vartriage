"""Path validation utilities for preventing path traversal (CWE-22).

All file path operations in the library should resolve paths via these
helpers before opening files or passing paths to external tools.
"""

from __future__ import annotations

from pathlib import Path


def resolve_path(path: Path) -> Path:
    """Resolve a path to its canonical absolute form.

    Follows symlinks, eliminates '..' components, and returns an
    absolute path. Used before any open() or subprocess call to
    prevent path traversal attacks.

    Parameters
    ----------
    path : Path
        User-provided or config-provided file path.

    Returns
    -------
    Path
        Canonicalized absolute path.
    """
    return path.resolve()


def validate_readable_path(path: Path, label: str = "File") -> Path:
    """Resolve path and verify it exists and is a regular file.

    Parameters
    ----------
    path : Path
        Path to validate.
    label : str
        Description for error messages (e.g., "GTF annotation").

    Returns
    -------
    Path
        Resolved, validated path.

    Raises
    ------
    FileNotFoundError
        If the resolved path does not exist or is not a file.
    ValueError
        If the path contains traversal sequences.
    """
    if ".." in path.parts:
        raise ValueError(f"{label} path contains traversal: {path}")
    resolved = path.resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return resolved


def validate_writable_path(path: Path, label: str = "Output") -> Path:
    """Resolve path for writing, creating parent directories as needed.

    Parameters
    ----------
    path : Path
        Destination path for writing.
    label : str
        Description for error messages.

    Returns
    -------
    Path
        Resolved path with parent directory created.

    Raises
    ------
    ValueError
        If the path contains traversal sequences.
    """
    if ".." in path.parts:
        raise ValueError(f"{label} path contains traversal: {path}")
    resolved = path.resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved
