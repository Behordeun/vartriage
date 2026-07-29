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

    Raises
    ------
    ValueError
        If the original path contains '..' traversal sequences.
    """
    _reject_traversal(path)
    return path.resolve()


def safe_read_path(path: Path, label: str = "File") -> Path:
    """Validate and resolve a path for reading.

    Rejects traversal sequences, resolves to absolute, and verifies
    existence. Use before every open(..., 'r') on user-provided paths.

    Parameters
    ----------
    path : Path
        Path to validate.
    label : str
        Description for error messages.

    Returns
    -------
    Path
        Resolved, validated path safe to open for reading.

    Raises
    ------
    FileNotFoundError
        If the resolved path does not exist.
    ValueError
        If the path contains traversal sequences.
    """
    _reject_traversal(path)
    resolved = path.resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return resolved


def safe_write_path(path: Path, label: str = "Output") -> Path:
    """Validate and resolve a path for writing.

    Rejects traversal sequences, resolves to absolute, creates parent
    directories. Use before every open(..., 'w') on user-provided paths.

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
    _reject_traversal(path)
    resolved = path.resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


# Keep old names for backward compatibility
validate_readable_path = safe_read_path
validate_writable_path = safe_write_path


def _reject_traversal(path: Path) -> None:
    """Raise ValueError if path contains directory traversal."""
    path_str = str(path)
    if ".." in path.parts or "/.." in path_str or "\\.." in path_str:
        raise ValueError(
            f"Path contains directory traversal sequence: {path}"
        )
