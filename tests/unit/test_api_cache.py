"""Unit tests for the SQLite response cache."""

from __future__ import annotations

import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from vartriage.api._cache import ResponseCache


@pytest.fixture
def cache(tmp_path: Path) -> ResponseCache:
    """Fresh cache instance with short TTL for testing."""
    return ResponseCache(db_path=tmp_path / "test_cache.db", default_ttl_days=1)


@pytest.fixture
def pinned_cache(tmp_path: Path) -> ResponseCache:
    """Cache with TTL disabled (pinned mode)."""
    return ResponseCache(db_path=tmp_path / "pinned_cache.db", default_ttl_days=-1)


class TestPutAndGet:
    """Basic store and retrieve operations."""

    def test_put_and_get_returns_stored_value(self, cache: ResponseCache) -> None:
        cache.put(
            key="vep:grch38:chr22:100:A:T",
            value={"consequence": "missense_variant", "gene": "BRCA1"},
            source="vep",
            genome_build="grch38",
        )
        result = cache.get("vep:grch38:chr22:100:A:T")
        assert result == {"consequence": "missense_variant", "gene": "BRCA1"}

    def test_get_missing_key_returns_none(self, cache: ResponseCache) -> None:
        assert cache.get("nonexistent:key") is None

    def test_put_overwrites_existing_entry(self, cache: ResponseCache) -> None:
        cache.put(key="k1", value={"v": 1}, source="vep", genome_build="grch38")
        cache.put(key="k1", value={"v": 2}, source="vep", genome_build="grch38")
        assert cache.get("k1") == {"v": 2}

    def test_stores_source_version(self, cache: ResponseCache) -> None:
        cache.put(
            key="k1",
            value={"data": "test"},
            source="vep",
            genome_build="grch38",
            source_version="Ensembl 112",
        )
        # source_version is stored but not exposed via get()
        # (it's for audit/stats purposes)
        assert cache.get("k1") == {"data": "test"}


class TestTTLExpiry:
    """Cache entry expiration behavior."""

    def test_expired_entry_returns_none(self, tmp_path: Path) -> None:
        cache = ResponseCache(db_path=tmp_path / "expiry.db", default_ttl_days=1)
        cache.put(key="k1", value={"x": 1}, source="vep", genome_build="grch38")

        # Manually backdate the expiry
        conn = cache._get_connection()
        past = (datetime.now(timezone.utc) - timedelta(days=2)).isoformat()
        conn.execute("UPDATE api_cache SET expires_at = ? WHERE key = ?", (past, "k1"))
        conn.commit()

        assert cache.get("k1") is None

    def test_expired_entry_is_deleted_on_read(self, tmp_path: Path) -> None:
        cache = ResponseCache(db_path=tmp_path / "evict.db", default_ttl_days=1)
        cache.put(key="k1", value={"x": 1}, source="vep", genome_build="grch38")

        conn = cache._get_connection()
        past = (datetime.now(timezone.utc) - timedelta(days=2)).isoformat()
        conn.execute("UPDATE api_cache SET expires_at = ? WHERE key = ?", (past, "k1"))
        conn.commit()

        cache.get("k1")  # triggers lazy eviction
        count = conn.execute(
            "SELECT COUNT(*) FROM api_cache WHERE key = ?", ("k1",)
        ).fetchone()[0]
        assert count == 0

    def test_custom_ttl_per_entry(self, cache: ResponseCache) -> None:
        cache.put(
            key="short",
            value={"v": 1},
            source="clinvar",
            genome_build="grch38",
            ttl_days=30,
        )
        # Entry with 30-day TTL should still be valid
        assert cache.get("short") == {"v": 1}


class TestPinnedMode:
    """Cache pinning (ttl_days=-1) for clinical reproducibility."""

    def test_pinned_entries_never_expire(self, pinned_cache: ResponseCache) -> None:
        pinned_cache.put(
            key="k1", value={"pinned": True}, source="vep", genome_build="grch38"
        )
        # Even if we check far in the future, entry persists
        assert pinned_cache.get("k1") == {"pinned": True}

    def test_pinned_cache_stats_shows_pinned_flag(
        self, pinned_cache: ResponseCache
    ) -> None:
        pinned_cache.put(key="k1", value={"x": 1}, source="vep", genome_build="grch38")
        stats = pinned_cache.stats()
        assert stats.is_pinned is True


class TestClear:
    """Cache clearing operations."""

    def test_clear_all_returns_deleted_count(self, cache: ResponseCache) -> None:
        for i in range(5):
            cache.put(key=f"k{i}", value={"i": i}, source="vep", genome_build="grch38")
        deleted = cache.clear()
        assert deleted == 5

    def test_clear_by_source(self, cache: ResponseCache) -> None:
        cache.put(key="v1", value={"x": 1}, source="vep", genome_build="grch38")
        cache.put(key="c1", value={"x": 2}, source="clinvar", genome_build="grch38")
        cache.put(key="v2", value={"x": 3}, source="vep", genome_build="grch38")

        deleted = cache.clear(source="vep")
        assert deleted == 2
        assert cache.get("c1") == {"x": 2}
        assert cache.get("v1") is None

    def test_clear_empty_cache_returns_zero(self, cache: ResponseCache) -> None:
        assert cache.clear() == 0


class TestEvictExpired:
    """Bulk eviction of expired entries."""

    def test_evict_removes_only_expired(self, tmp_path: Path) -> None:
        cache = ResponseCache(db_path=tmp_path / "evict_bulk.db", default_ttl_days=1)
        cache.put(key="fresh", value={"f": 1}, source="vep", genome_build="grch38")
        cache.put(key="stale", value={"s": 1}, source="vep", genome_build="grch38")

        conn = cache._get_connection()
        past = (datetime.now(timezone.utc) - timedelta(days=2)).isoformat()
        conn.execute(
            "UPDATE api_cache SET expires_at = ? WHERE key = ?", (past, "stale")
        )
        conn.commit()

        evicted = cache.evict_expired()
        assert evicted == 1
        assert cache.get("fresh") == {"f": 1}
        assert cache.get("stale") is None


class TestStats:
    """Cache statistics."""

    def test_empty_cache_stats(self, cache: ResponseCache) -> None:
        stats = cache.stats()
        assert stats.entry_count == 0
        assert stats.oldest_entry is None
        assert stats.newest_entry is None
        assert stats.entries_by_source == {}

    def test_stats_tracks_entries_by_source(self, cache: ResponseCache) -> None:
        cache.put(key="v1", value={}, source="vep", genome_build="grch38")
        cache.put(key="v2", value={}, source="vep", genome_build="grch38")
        cache.put(key="c1", value={}, source="clinvar", genome_build="grch38")

        stats = cache.stats()
        assert stats.entry_count == 3
        assert stats.entries_by_source == {"vep": 2, "clinvar": 1}

    def test_stats_reports_disk_size(self, cache: ResponseCache) -> None:
        cache.put(
            key="k1", value={"data": "x" * 100}, source="vep", genome_build="grch38"
        )
        stats = cache.stats()
        assert stats.disk_size_bytes > 0


class TestBuildKey:
    """Static key generation helper."""

    def test_build_key_format(self) -> None:
        key = ResponseCache.build_key("vep", "grch38", "chr22", 17818804, "G", "A")
        assert key == "vep:grch38:chr22:17818804:G:A"

    def test_build_key_deterministic(self) -> None:
        k1 = ResponseCache.build_key("clinvar", "grch37", "1", 12345, "AT", "A")
        k2 = ResponseCache.build_key("clinvar", "grch37", "1", 12345, "AT", "A")
        assert k1 == k2


class TestCorruptedEntries:
    """Handling of malformed cache data."""

    def test_corrupted_json_returns_none_and_cleans_up(self, tmp_path: Path) -> None:
        cache = ResponseCache(db_path=tmp_path / "corrupt.db", default_ttl_days=1)

        # Insert a valid entry then corrupt the JSON
        cache.put(key="bad", value={"ok": True}, source="vep", genome_build="grch38")
        conn = cache._get_connection()
        conn.execute(
            "UPDATE api_cache SET response_json = ? WHERE key = ?",
            ("{invalid json", "bad"),
        )
        conn.commit()

        assert cache.get("bad") is None
        # Entry should be removed
        count = conn.execute(
            "SELECT COUNT(*) FROM api_cache WHERE key = ?", ("bad",)
        ).fetchone()[0]
        assert count == 0
