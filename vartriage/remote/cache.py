"""SQLite-backed score cache for remote tabix lookups.

Stores individual variant scores keyed by (source, chrom, pos, ref, alt)
to avoid redundant network fetches on repeated runs. Supports configurable
TTL and pinned mode (TTL=-1) for clinical reproducibility.

This cache is distinct from the API ResponseCache which stores JSON blobs.
RemoteScoreCache stores scalar float scores with a flat schema optimized
for high-volume batch lookups.
"""

from __future__ import annotations

import logging
import sqlite3
import threading
import time
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

_SCHEMA_SQL = """\
CREATE TABLE IF NOT EXISTS remote_scores (
    source TEXT NOT NULL,
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    score REAL NOT NULL,
    fetched_at INTEGER NOT NULL,
    PRIMARY KEY (source, chrom, pos, ref, alt)
);
CREATE INDEX IF NOT EXISTS idx_remote_scores_fetched
    ON remote_scores(fetched_at);
"""


@dataclass(frozen=True)
class CacheStats:
    """Snapshot of remote score cache state."""

    entry_count: int
    disk_size_bytes: int
    entries_by_source: dict[str, int]


class RemoteScoreCache:
    """SQLite cache for remote tabix score lookups.

    Thread-safe via a per-instance lock. The database file is created
    lazily on first write. Expired entries are evicted on read (lazy
    eviction).

    Parameters
    ----------
    db_path : Path
        Path to the SQLite database file. Parent directories are
        created if absent.
    ttl_days : int
        Default time-to-live in days. Use -1 to pin entries
        indefinitely (never expire).
    """

    def __init__(self, db_path: Path, ttl_days: int = 30) -> None:
        self._db_path = Path(db_path).expanduser()
        self._ttl_days = ttl_days
        self._ttl_seconds = ttl_days * 86_400 if ttl_days > 0 else -1
        self._lock = threading.Lock()
        self._conn: sqlite3.Connection | None = None

    def get(
        self, source: str, chrom: str, pos: int, ref: str, alt: str
    ) -> float | None:
        """Retrieve a cached score. Returns None on miss or expiry.

        Expired entries are deleted on access (lazy eviction).
        """
        now = int(time.time())
        with self._lock:
            conn = self._ensure_connection()
            result = self._lookup_row(conn, source, chrom, pos, ref, alt, now)
            conn.commit()
            return result

    def get_batch(
        self,
        source: str,
        variants: list[tuple[str, int, str, str]],
    ) -> list[float | None]:
        """Batch lookup of cached scores.

        Expired entries are deleted on access (lazy eviction),
        consistent with single get() behavior.

        Parameters
        ----------
        source : str
            Score source identifier (e.g., "cadd-remote", "gnomad-remote").
        variants : list[tuple[str, int, str, str]]
            (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[float | None]
            Cached scores positionally matched. None for cache misses.
        """
        results: list[float | None] = [None] * len(variants)
        now = int(time.time())

        with self._lock:
            conn = self._ensure_connection()
            for i, (chrom, pos, ref, alt) in enumerate(variants):
                results[i] = self._lookup_row(conn, source, chrom, pos, ref, alt, now)
            conn.commit()

        return results

    def put(
        self, source: str, chrom: str, pos: int, ref: str, alt: str, score: float
    ) -> None:
        """Store a single score in the cache."""
        now = int(time.time())
        with self._lock:
            conn = self._ensure_connection()
            conn.execute(
                "INSERT OR REPLACE INTO remote_scores "
                "(source, chrom, pos, ref, alt, score, fetched_at) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (source, chrom, pos, ref, alt, score, now),
            )
            conn.commit()

    def put_batch(
        self,
        source: str,
        entries: list[tuple[str, int, str, str, float]],
    ) -> None:
        """Store multiple scores in a single transaction.

        Parameters
        ----------
        source : str
            Score source identifier.
        entries : list[tuple[str, int, str, str, float]]
            (chrom, pos, ref, alt, score) tuples to cache.
        """
        if not entries:
            return

        now = int(time.time())
        with self._lock:
            conn = self._ensure_connection()
            conn.executemany(
                "INSERT OR REPLACE INTO remote_scores "
                "(source, chrom, pos, ref, alt, score, fetched_at) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                [
                    (source, chrom, pos, ref, alt, score, now)
                    for chrom, pos, ref, alt, score in entries
                ],
            )
            conn.commit()

    def evict_expired(self) -> int:
        """Remove all expired entries. Returns count evicted.

        No-op when TTL is pinned (-1).
        """
        if self._ttl_seconds == -1:
            return 0

        cutoff = int(time.time()) - self._ttl_seconds
        with self._lock:
            conn = self._ensure_connection()
            cursor = conn.execute(
                "DELETE FROM remote_scores WHERE fetched_at < ?", (cutoff,)
            )
            conn.commit()
            return cursor.rowcount

    def clear(self, source: str | None = None) -> int:
        """Delete cache entries. Returns count deleted.

        Parameters
        ----------
        source : str | None
            If provided, only clear entries for this source.
            None clears everything.
        """
        with self._lock:
            conn = self._ensure_connection()
            if source is None:
                cursor = conn.execute("DELETE FROM remote_scores")
            else:
                cursor = conn.execute(
                    "DELETE FROM remote_scores WHERE source = ?", (source,)
                )
            conn.commit()
            deleted = cursor.rowcount
            logger.info(
                "Cleared %d remote score cache entries (source=%s)",
                deleted,
                source or "all",
            )
            return deleted

    def stats(self) -> CacheStats:
        """Compute cache statistics."""
        with self._lock:
            conn = self._ensure_connection()
            entry_count = conn.execute("SELECT COUNT(*) FROM remote_scores").fetchone()[
                0
            ]
            source_rows = conn.execute(
                "SELECT source, COUNT(*) FROM remote_scores GROUP BY source"
            ).fetchall()

        disk_size = 0
        if self._db_path.exists():
            disk_size = self._db_path.stat().st_size

        return CacheStats(
            entry_count=entry_count,
            disk_size_bytes=disk_size,
            entries_by_source=dict(source_rows),
        )

    def close(self) -> None:
        """Close the database connection."""
        with self._lock:
            if self._conn is not None:
                self._conn.close()
                self._conn = None

    def _is_expired(self, fetched_at: int, now: int | None = None) -> bool:
        """Check whether a cached entry has expired."""
        if self._ttl_seconds == -1:
            return False
        if now is None:
            now = int(time.time())
        return (now - fetched_at) > self._ttl_seconds

    def _lookup_row(
        self,
        conn: sqlite3.Connection,
        source: str,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        now: int,
    ) -> float | None:
        """Look up a single row and lazily evict if expired.

        Caller must hold self._lock. A single commit() after all
        lookups in a batch is sufficient.
        """
        cursor = conn.execute(
            "SELECT score, fetched_at FROM remote_scores "
            "WHERE source = ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?",
            (source, chrom, pos, ref, alt),
        )
        row = cursor.fetchone()
        if row is None:
            return None

        score, fetched_at = row
        if self._is_expired(fetched_at, now):
            conn.execute(
                "DELETE FROM remote_scores "
                "WHERE source = ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?",
                (source, chrom, pos, ref, alt),
            )
            return None

        return float(score)

    def _ensure_connection(self) -> sqlite3.Connection:
        """Lazy-init the SQLite connection and schema. Caller holds lock."""
        if self._conn is not None:
            return self._conn

        self._db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self._db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA busy_timeout=5000")
        self._conn.executescript(_SCHEMA_SQL)
        return self._conn
