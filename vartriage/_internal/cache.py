"""Pickle-based file caching with mtime invalidation and atomic writes.

Provides a shared cache infrastructure for serializing parsed reference
data (GTF interval trees, score dictionaries) to disk. Uses mtime-based
invalidation and version stamping to detect stale or incompatible caches.

All public functions handle errors gracefully. Cache failures never
propagate exceptions to callers.
"""

from __future__ import annotations

import logging
import os
import pickle
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

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
        source_path with '.vartriage.cache' appended.
    """
    return Path(str(source_path) + ".vartriage.cache")


def try_load_cache(source_path: Path) -> Optional[Any]:
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
        with open(cp, "rb") as f:
            envelope: CacheEnvelope = pickle.load(f)  # noqa: S301
    except (OSError, PermissionError) as exc:
        logger.warning(
            "Cannot read cache file %s: %s", cp, exc
        )
        _delete_cache(cp)
        return None
    except (pickle.UnpicklingError, Exception) as exc:
        logger.warning(
            "Failed to deserialize cache %s: %s", cp, exc
        )
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
        logger.warning(
            "Cannot stat source file %s: %s", source_path, exc
        )
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
    cp = cache_path_for(source_path)

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
        tmp_fd = tempfile.NamedTemporaryFile(
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
        logger.warning(
            "Failed to write cache for %s: %s", source_path, exc
        )
        if tmp_fd is not None:
            try:
                tmp_fd.close()
            except OSError:
                pass
        if tmp_path is not None:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


def _delete_cache(cache_path: Path) -> None:
    """Best-effort deletion of a cache file."""
    try:
        cache_path.unlink()
        logger.debug("Deleted invalid cache %s", cache_path)
    except OSError:
        pass
