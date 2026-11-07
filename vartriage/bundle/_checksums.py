"""Streaming SHA-256 checksum computation and verification."""

from __future__ import annotations

import hashlib
from pathlib import Path

from vartriage._internal.path_safety import safe_read_path


def compute_sha256(path: Path, chunk_size: int = 65536) -> str:
    """Compute SHA-256 hex digest for a file.

    Reads the file in chunks to avoid loading large files into memory.

    Parameters
    ----------
    path : Path
        Path to the file to checksum.
    chunk_size : int
        Read buffer size in bytes (default: 64 KB).

    Returns
    -------
    str
        Hex-encoded SHA-256 digest prefixed with "sha256:".

    Raises
    ------
    OSError
        If the file cannot be read.
    """
    sha256 = hashlib.sha256()
    path = safe_read_path(path, "Checksum target")
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            sha256.update(chunk)
    return f"sha256:{sha256.hexdigest()}"


def verify_checksum(path: Path, expected: str, *, strict: bool = False) -> bool:
    """Verify a file's SHA-256 checksum against an expected value.

    Parameters
    ----------
    path : Path
        Path to the file to verify.
    expected : str
        Expected checksum in "sha256:<hex>" format. If empty
        string, verification is skipped (returns True) unless
        strict is True.
    strict : bool, optional
        When True, an empty or missing expected checksum raises
        UnverifiedChecksumError rather than returning True. Default is False.

    Returns
    -------
    bool
        True if checksum matches or (when strict is False) expected is empty.

    Raises
    ------
    UnverifiedChecksumError
        If strict is True and expected is empty.
    """
    if not expected:
        if strict:
            raise UnverifiedChecksumError(
                f"Checksum verification failed for {path}: no expected checksum provided"
            )
        return True

    actual = compute_sha256(path)
    return actual.lower() == expected.lower()


class ChecksumMismatchError(Exception):
    """Raised when a file's checksum doesn't match the expected value."""

    def __init__(self, path: Path, expected: str, actual: str) -> None:
        self.path = path
        self.expected = expected
        self.actual = actual
        super().__init__(
            f"Checksum mismatch for {path}:\n"
            f"  expected: {expected}\n"
            f"  actual:   {actual}"
        )


class UnverifiedChecksumError(Exception):
    """Raised when strict checksum verification fails due to missing checksum metadata."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
