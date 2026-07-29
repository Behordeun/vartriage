"""Path validation utilities for preventing path traversal (CWE-22).

Two levels of safety depending on trust level of the input:

- resolve_path(): for CLI/user-supplied paths. Resolves to canonical absolute
  form (symlinks followed, '..' normalized). Does NOT reject '..' — users
  legitimately pass relative paths like '../sample.vcf'. Use this for any
  path that comes from argparse or config files.

- safe_read_path() / safe_write_path(): for untrusted programmatic inputs
  (e.g. filenames derived from URLs, API responses, bundle downloads).
  Rejects '..' traversal sequences before resolving. Use these whenever
  the path component is computed from external data rather than typed by
  the user.

Usage guidelines:
  CLI args (--vcf, --output, etc.)     -> resolve_path()
  Bundle storage paths from user       -> _resolve_storage_path() in bundle/cli.py
  Filenames from URL splits            -> _sanitize_filename() in bundle/cli.py
  Transformer source/dest              -> safe_read_path() / safe_write_path()
  Cache file paths (internal)          -> safe_read_path() / safe_write_path()
"""

from __future__ import annotations

from pathlib import Path


def resolve_path(path: Path) -> Path:
    """Resolve a path to its canonical absolute form.

    Follows symlinks, normalizes '..' components, and returns an
    absolute path. Used for CLI-supplied paths where relative
    traversal (../file.vcf) is legitimate.

    Parameters
    ----------
    path : Path
        User-provided file path (may be relative).

    Returns
    -------
    Path
        Canonicalized absolute path.
    """
    return path.resolve()


def safe_read_path(path: Path, label: str = "File") -> Path:
    """Validate and resolve a path for reading (untrusted inputs).

    Rejects traversal sequences, resolves to absolute, and verifies
    existence. Use for paths derived from untrusted sources (URLs,
    API responses, computed filenames).

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
    """Validate and resolve a path for writing (untrusted inputs).

    Rejects traversal sequences, resolves to absolute, creates parent
    directories. Use for paths derived from untrusted sources.

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
