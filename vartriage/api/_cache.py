"""SQLite-backed response cache with per-source TTL expiry.

Caches API responses keyed by (source, build, chrom, pos, ref, alt)
to avoid redundant network calls. Supports configurable TTL per entry,
source version tracking for reproducibility, and pinned mode (TTL=-1)
for clinical labs needing bit-identical results across runs.
"""

from __future__ import annotations

import json
import logging
import sqlite3
import threading
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_SCHEMA_SQL = """\
CREATE TABLE IF NOT EXISTS api_cache (
    key TEXT PRIMARY KEY,
    source TEXT NOT NULL,
    genome_build TEXT NOT NULL,
    response_json TEXT NOT NULL,
    source_version TEXT,
    created_at TEXT NOT NULL,
    expires_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_cache_expires ON api_cache(expires_at);
CREATE INDEX IF NOT EXISTS idx_cache_source ON api_cache(source);
"""

_PINNED_EXPIRY = "9999-12-31T23:59:59+00:00"


@dataclass(frozen=True)
class CacheStats:
    """Snapshot of cache state for observability."""

    entry_count: int
    disk_size_bytes: int
    oldest_entry: str | None
    newest_entry: str | None
    is_pinned: bool
    entries_by_source: dict[str, int]


class ResponseCache:
    """SQLite response cache with TTL and source version tracking.

    Thread-safe via a per-instance lock around all database operations.
    The database file is created lazily on first write.

    Parameters
    ----------
    db_path
        Path to the SQLite database file. Parent directories created if absent.
    default_ttl_days
        Default time-to-live in days. Use -1 to pin entries indefinitely.
    """

    def __init__(self, db_path: Path, default_ttl_days: int = 7) -> None:
        self._db_path = Path(db_path).expanduser()
        self._default_ttl_days = default_ttl_days
        self._lock = threading.Lock()
        self._conn: sqlite3.Connection | None = None

    def get(self, key: str) -> dict[str, Any] | None:
        """Retrieve a cached response, returning None if expired or absent.

        Expired entries are deleted on access (lazy eviction).
        """
        with self._lock:
            conn = self._get_connection()
            cursor = conn.execute(
                "SELECT response_json, expires_at FROM api_cache WHERE key = ?",
                (key,),
            )
            row = cursor.fetchone()

            if row is None:
                return None

            response_json, expires_at = row
            now = datetime.now(timezone.utc).isoformat()

            # Check expiry (pinned entries have far-future expiry)
            if expires_at < now:
                conn.execute("DELETE FROM api_cache WHERE key = ?", (key,))
                conn.commit()
                return None

            try:
                result: dict[str, Any] = json.loads(response_json)
                return result
            except (json.JSONDecodeError, TypeError):
                # Corrupted entry, remove it
                conn.execute("DELETE FROM api_cache WHERE key = ?", (key,))
                conn.commit()
                logger.warning("Removed corrupted cache entry: %s", key)
                return None

    def put(
        self,
        key: str,
        value: dict[str, Any],
        source: str,
        genome_build: str,
        source_version: str | None = None,
        ttl_days: int | None = None,
    ) -> None:
        """Store a response in the cache.

        Parameters
        ----------
        key
            Cache key (typically source:build:chrom:pos:ref:alt).
        value
            Response data to serialize as JSON.
        source
            API source name (e.g., "vep", "clinvar").
        genome_build
            Genome build (e.g., "grch38").
        source_version
            Database version from response headers (e.g., "Ensembl 112").
        ttl_days
            Override TTL for this entry. None uses the instance default.
            -1 means pinned (never expires).
        """
        effective_ttl = ttl_days if ttl_days is not None else self._default_ttl_days
        now = datetime.now(timezone.utc)
        created_at = now.isoformat()

        if effective_ttl == -1:
            expires_at = _PINNED_EXPIRY
        else:
            expires_at = (now + timedelta(days=effective_ttl)).isoformat()

        response_json = json.dumps(value, ensure_ascii=False, separators=(",", ":"))

        with self._lock:
            conn = self._get_connection()
            conn.execute(
                """
                INSERT OR REPLACE INTO api_cache
                    (key, source, genome_build, response_json, source_version, created_at, expires_at)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    key,
                    source,
                    genome_build,
                    response_json,
                    source_version,
                    created_at,
                    expires_at,
                ),
            )
            conn.commit()

    def clear(self, source: str | None = None) -> int:
        """Delete cache entries. Returns count of deleted rows.

        Parameters
        ----------
        source
            If provided, only clear entries for this source. None clears everything.
        """
        with self._lock:
            conn = self._get_connection()
            if source is None:
                cursor = conn.execute("DELETE FROM api_cache")
            else:
                cursor = conn.execute(
                    "DELETE FROM api_cache WHERE source = ?", (source,)
                )
            conn.commit()
            deleted = cursor.rowcount
            logger.info(
                "Cleared %d cache entries (source=%s)", deleted, source or "all"
            )
            return deleted

    def evict_expired(self) -> int:
        """Remove all expired entries. Returns count evicted."""
        now = datetime.now(timezone.utc).isoformat()
        with self._lock:
            conn = self._get_connection()
            cursor = conn.execute("DELETE FROM api_cache WHERE expires_at < ?", (now,))
            conn.commit()
            return cursor.rowcount

    def stats(self) -> CacheStats:
        """Compute cache statistics."""
        with self._lock:
            conn = self._get_connection()

            entry_count = conn.execute("SELECT COUNT(*) FROM api_cache").fetchone()[0]

            oldest_row = conn.execute(
                "SELECT created_at FROM api_cache ORDER BY created_at ASC LIMIT 1"
            ).fetchone()
            newest_row = conn.execute(
                "SELECT created_at FROM api_cache ORDER BY created_at DESC LIMIT 1"
            ).fetchone()

            source_rows = conn.execute(
                "SELECT source, COUNT(*) FROM api_cache GROUP BY source"
            ).fetchall()

        # Disk size
        disk_size = 0
        if self._db_path.exists():
            disk_size = self._db_path.stat().st_size

        is_pinned = self._default_ttl_days == -1

        return CacheStats(
            entry_count=entry_count,
            disk_size_bytes=disk_size,
            oldest_entry=oldest_row[0] if oldest_row else None,
            newest_entry=newest_row[0] if newest_row else None,
            is_pinned=is_pinned,
            entries_by_source=dict(source_rows),
        )

    def close(self) -> None:
        """Close the database connection."""
        with self._lock:
            if self._conn is not None:
                self._conn.close()
                self._conn = None

    def _get_connection(self) -> sqlite3.Connection:
        """Lazy-init the SQLite connection and schema. Caller holds lock."""
        if self._conn is not None:
            return self._conn

        self._db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self._db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA busy_timeout=5000")
        self._conn.executescript(_SCHEMA_SQL)
        return self._conn

    @staticmethod
    def build_key(
        source: str, genome_build: str, chrom: str, pos: int, ref: str, alt: str
    ) -> str:
        """Build a deterministic cache key from variant coordinates.

        Format: source:build:chrom:pos:ref:alt
        """
        return f"{source}:{genome_build}:{chrom}:{pos}:{ref}:{alt}"
