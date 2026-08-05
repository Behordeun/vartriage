"""Pickle-based file caching with mtime invalidation and atomic writes.

Provides a shared cache infrastructure for serializing parsed reference
data (GTF interval trees, score dictionaries) to disk. Uses mtime-based
invalidation and version stamping to detect stale or incompatible caches.

All public functions handle errors gracefully. Cache failures never
propagate exceptions to callers.
"""

from __future__ import annotations

import contextlib
import logging
import os
import pickle
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from vartriage._internal.path_safety import safe_read_path, safe_write_path

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class CacheEnvelope:
    """Metadata wrapper around cached data.

    Attributes
    ----------
    vartriage_version : str
        Package version at serialization time.
    python_version : str
        Python major.minor at serialization time.
    source_mtime : float
        Source file mtime at serialization time.
    data : Any
        The actual cached object.
    """

    vartriage_version: str
    python_version: str
    source_mtime: float
    data: Any


def _current_python_version() -> str:
    """Return 'major.minor' string for current interpreter."""
    return f"{sys.version_info.major}.{sys.version_info.minor}"


def _current_vartriage_version() -> str:
    """Return current vartriage package version."""
    from vartriage import __version__

    return __version__


def cache_path_for(source_path: Path) -> Path:
    """Compute cache file path for a given source file.

    Parameters
    ----------
    source_path : Path
        Path to the original data file.

    Returns
    -------
    Path
        A sibling file in the same directory with '.vartriage.cache' appended
        to the source filename.
    """
    resolved = source_path.resolve()
    return resolved.parent / (resolved.name + ".vartriage.cache")


def try_load_cache(source_path: Path) -> Any | None:
    """Attempt to load cached data for source_path.

    Returns the cached data if:
    - Cache file exists and is readable
    - Pickle deserialization succeeds
    - vartriage_version matches current version
    - python_version matches current major.minor
    - source_mtime matches source file's current mtime

    On any failure, logs a warning, deletes the invalid cache
    (if possible), and returns None.

    Parameters
    ----------
    source_path : Path
        Path to the original data file.

    Returns
    -------
    Optional[Any]
        The cached data, or None on miss/failure.
    """
    cp = cache_path_for(source_path)

    if not cp.exists():
        return None

    try:
        cp = safe_read_path(cp, "Cache file")
        with open(cp, "rb") as f:
            envelope: CacheEnvelope = pickle.load(f)  # noqa: S301
    except OSError as exc:
        logger.warning("Cannot read cache file %s: %s", cp, exc)
        _delete_cache(cp)
        return None
    except Exception as exc:
        logger.warning("Failed to deserialize cache %s: %s", cp, exc)
        _delete_cache(cp)
        return None

    if not isinstance(envelope, CacheEnvelope):
        logger.warning(
            "Cache %s contains unexpected type %s",
            cp,
            type(envelope).__name__,
        )
        _delete_cache(cp)
        return None

    current_vt = _current_vartriage_version()
    if envelope.vartriage_version != current_vt:
        logger.info(
            "Cache %s has vartriage version %s, current is %s",
            cp,
            envelope.vartriage_version,
            current_vt,
        )
        _delete_cache(cp)
        return None

    current_py = _current_python_version()
    if envelope.python_version != current_py:
        logger.info(
            "Cache %s has Python version %s, current is %s",
            cp,
            envelope.python_version,
            current_py,
        )
        _delete_cache(cp)
        return None

    try:
        current_mtime = source_path.stat().st_mtime
    except OSError as exc:
        logger.warning("Cannot stat source file %s: %s", source_path, exc)
        return None

    if envelope.source_mtime != current_mtime:
        logger.debug(
            "Cache %s is stale (mtime %s vs %s)",
            cp,
            envelope.source_mtime,
            current_mtime,
        )
        _delete_cache(cp)
        return None

    logger.debug("Cache hit for %s", source_path)
    return envelope.data


def try_write_cache(source_path: Path, data: Any) -> None:
    """Serialize data to cache file with atomic write.

    Writes to a temporary file in the same directory, then
    renames atomically via os.replace(). On any failure, logs
    a warning and returns without raising.

    Parameters
    ----------
    source_path : Path
        Path to the original data file (determines cache path
        and mtime to stamp).
    data : Any
        The object to serialize via pickle.
    """
    cp = safe_write_path(cache_path_for(source_path), "Cache write")

    try:
        source_mtime = source_path.stat().st_mtime
    except OSError as exc:
        logger.warning(
            "Cannot stat source file %s for cache write: %s",
            source_path,
            exc,
        )
        return

    envelope = CacheEnvelope(
        vartriage_version=_current_vartriage_version(),
        python_version=_current_python_version(),
        source_mtime=source_mtime,
        data=data,
    )

    tmp_fd = None
    tmp_path = None
    try:
        tmp_fd = tempfile.NamedTemporaryFile(  # noqa: SIM115
            dir=cp.parent,
            prefix=".vartriage_cache_",
            suffix=".tmp",
            delete=False,
        )
        tmp_path = tmp_fd.name
        pickle.dump(envelope, tmp_fd, protocol=pickle.HIGHEST_PROTOCOL)
        tmp_fd.close()
        tmp_fd = None
        os.replace(tmp_path, cp)
        logger.debug("Cache written for %s", source_path)
    except (OSError, pickle.PicklingError) as exc:
        logger.warning("Failed to write cache for %s: %s", source_path, exc)
        if tmp_fd is not None:
            with contextlib.suppress(OSError):
                tmp_fd.close()
        if tmp_path is not None:
            with contextlib.suppress(OSError):
                os.unlink(tmp_path)


def _delete_cache(cache_path: Path) -> None:
    """Best-effort deletion of a cache file."""
    try:
        cache_path.unlink()
        logger.debug("Deleted invalid cache %s", cache_path)
    except OSError:
        pass
