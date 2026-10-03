"""Shared pytest fixtures and configuration for the test suite."""

from __future__ import annotations

import os
import sqlite3
import tempfile
from collections.abc import Generator
from pathlib import Path

import pytest
from hypothesis import settings

from vartriage.models.config import PrioritizationConfig, QualityFilterConfig

settings.register_profile("ci", max_examples=500, deadline=None)
settings.register_profile("dev", max_examples=50, deadline=None)
settings.register_profile("debug", max_examples=10, deadline=None)
settings.load_profile("dev")


def _can_host_wal_sqlite(directory: Path) -> bool:
    """Return True when a WAL-mode SQLite database opens under ``directory``.

    Some sandboxed or overlay filesystems (seen when TMPDIR points at an
    agent scratch mount) reject the -wal/-shm sidecar files SQLite creates in
    WAL mode, surfacing as "unable to open database file". The remote score
    cache uses WAL, so a temp base that cannot host it makes those tests fail
    for reasons unrelated to the code under test.
    """
    probe = directory / ".wal_probe.db"
    try:
        conn = sqlite3.connect(str(probe))
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("CREATE TABLE t (a)")
        conn.close()
        return True
    except sqlite3.OperationalError:
        return False
    finally:
        for suffix in ("", "-wal", "-shm"):
            Path(f"{probe}{suffix}").unlink(missing_ok=True)


def pytest_configure(config: pytest.Config) -> None:
    """Redirect the temp base to a WAL-capable location when needed.

    pytest derives tmp_path from TMPDIR. When that mount cannot host a WAL
    SQLite database, point --basetemp at a filesystem that can, so cache-backed
    tests exercise the code rather than the sandbox's filesystem limits.
    """
    if config.option.basetemp is not None:
        return

    inherited = Path(tempfile.gettempdir())
    if _can_host_wal_sqlite(inherited):
        return

    for candidate in (Path("/var/tmp"), Path("/tmp"), Path.home() / ".cache"):
        if candidate.is_dir() and os.access(candidate, os.W_OK):
            fallback = candidate / "vartriage-pytest-tmp"
            fallback.mkdir(parents=True, exist_ok=True)
            if _can_host_wal_sqlite(fallback):
                config.option.basetemp = str(fallback)
                return


@pytest.fixture
def tmp_dir() -> Generator[Path, None, None]:
    """Provide a temporary directory cleaned up after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def default_quality_config() -> QualityFilterConfig:
    """Standard quality filter configuration with default threshold."""
    return QualityFilterConfig(min_qual=20.0)


@pytest.fixture
def default_prioritization_config() -> PrioritizationConfig:
    """Standard prioritization configuration with default thresholds."""
    return PrioritizationConfig(max_allele_frequency=0.01)
