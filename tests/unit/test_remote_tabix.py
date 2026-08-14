"""Unit tests for the remote tabix score backend modules."""

from __future__ import annotations

import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vartriage.remote.cache import RemoteScoreCache
from vartriage.remote.circuit_breaker import CircuitBreaker, CircuitState
from vartriage.remote.config import RemoteTabixConfig
from vartriage.remote.presets import (
    get_preset,
    is_preset_name,
    list_presets,
    resolve_preset,
)

# ============================================================
# RemoteTabixConfig tests
# ============================================================


class TestRemoteTabixConfig:
    """Validation and property behavior for RemoteTabixConfig."""

    def test_defaults_have_no_remote_active(self) -> None:
        config = RemoteTabixConfig()
        assert not config.is_cadd_active
        assert not config.is_gnomad_active
        assert not config.has_any_remote

    def test_cadd_url_activates_cadd(self) -> None:
        config = RemoteTabixConfig(cadd_remote_url="cadd-v1.7-grch38")
        assert config.is_cadd_active
        assert not config.is_gnomad_active
        assert config.has_any_remote

    def test_gnomad_url_activates_gnomad(self) -> None:
        config = RemoteTabixConfig(gnomad_remote_url="gnomad-exomes-v4-grch38")
        assert not config.is_cadd_active
        assert config.is_gnomad_active
        assert config.has_any_remote

    def test_both_urls_active(self) -> None:
        config = RemoteTabixConfig(
            cadd_remote_url="cadd-v1.7-grch38",
            gnomad_remote_url="gnomad-exomes-v4-grch38",
        )
        assert config.is_cadd_active
        assert config.is_gnomad_active

    def test_cache_ttl_minus_one_is_pinned(self) -> None:
        config = RemoteTabixConfig(cache_ttl_days=-1)
        assert config.cache_ttl_days == -1

    def test_cache_ttl_below_minus_one_raises(self) -> None:
        with pytest.raises(ValueError, match="cache_ttl_days must be >= -1"):
            RemoteTabixConfig(cache_ttl_days=-2)

    def test_connect_timeout_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="connect_timeout must be > 0"):
            RemoteTabixConfig(connect_timeout=0)

    def test_read_timeout_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="read_timeout must be > 0"):
            RemoteTabixConfig(read_timeout=-1.0)

    def test_max_retries_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="max_retries must be >= 0"):
            RemoteTabixConfig(max_retries=-1)

    def test_batch_window_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="batch_window_bp must be >= 1"):
            RemoteTabixConfig(batch_window_bp=0)

    def test_frozen_dataclass(self) -> None:
        config = RemoteTabixConfig()
        with pytest.raises(AttributeError):
            config.cache_ttl_days = 60  # type: ignore[misc]


# ============================================================
# Presets tests
# ============================================================


class TestPresets:
    """Named preset registry resolution."""

    def test_resolve_known_cadd_preset(self) -> None:
        url = resolve_preset("cadd-v1.7-grch38")
        assert "krishna.gs.washington.edu" in url
        assert "GRCh38" in url

    def test_resolve_known_gnomad_preset(self) -> None:
        url = resolve_preset("gnomad-exomes-v4-grch38")
        assert "{chrom}" in url
        assert "gnomad" in url

    def test_resolve_raw_url_passes_through(self) -> None:
        raw = "https://example.com/scores.tsv.gz"
        assert resolve_preset(raw) == raw

    def test_resolve_http_url_passes_through(self) -> None:
        raw = "http://internal.server/scores.tsv.gz"
        assert resolve_preset(raw) == raw

    def test_resolve_unknown_preset_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown preset 'nonexistent'"):
            resolve_preset("nonexistent")

    def test_error_message_lists_available_presets(self) -> None:
        with pytest.raises(ValueError, match="cadd-v1.7-grch38"):
            resolve_preset("bogus")

    def test_get_preset_returns_entry(self) -> None:
        entry = get_preset("cadd-v1.7-grch38")
        assert entry is not None
        assert entry.source == "cadd"
        assert entry.genome_build == "grch38"

    def test_get_preset_returns_none_for_unknown(self) -> None:
        assert get_preset("nonexistent") is None

    def test_list_presets_returns_all(self) -> None:
        presets = list_presets()
        assert len(presets) == 4

    def test_list_presets_filtered_by_source(self) -> None:
        cadd_presets = list_presets(source="cadd")
        gnomad_presets = list_presets(source="gnomad")
        assert all(p.source == "cadd" for p in cadd_presets)
        assert all(p.source == "gnomad" for p in gnomad_presets)
        assert len(cadd_presets) == 3
        assert len(gnomad_presets) == 1

    def test_list_presets_sorted_by_name(self) -> None:
        presets = list_presets()
        names = [p.name for p in presets]
        assert names == sorted(names)

    def test_is_preset_name_true(self) -> None:
        assert is_preset_name("cadd-v1.7-grch38")

    def test_is_preset_name_false_for_url(self) -> None:
        assert not is_preset_name("https://example.com/file.tsv.gz")

    def test_is_preset_name_false_for_unknown(self) -> None:
        assert not is_preset_name("unknown-preset")


# ============================================================
# RemoteScoreCache tests
# ============================================================


@pytest.fixture
def score_cache(tmp_path: Path) -> RemoteScoreCache:
    """Fresh score cache with 1-day TTL."""
    c = RemoteScoreCache(db_path=tmp_path / "test_remote_cache.db", ttl_days=1)
    yield c
    c.close()


@pytest.fixture
def pinned_score_cache(tmp_path: Path) -> RemoteScoreCache:
    """Score cache with pinned mode (never expires)."""
    c = RemoteScoreCache(db_path=tmp_path / "pinned_remote_cache.db", ttl_days=-1)
    yield c
    c.close()


class TestRemoteScoreCache:
    """SQLite score cache with TTL and batch operations."""

    def test_put_and_get_single(self, score_cache: RemoteScoreCache) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        result = score_cache.get("cadd-remote", "chr1", 100, "A", "T")
        assert result == 23.5

    def test_get_missing_returns_none(self, score_cache: RemoteScoreCache) -> None:
        assert score_cache.get("cadd-remote", "chr1", 999, "G", "C") is None

    def test_put_overwrites(self, score_cache: RemoteScoreCache) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 10.0)
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 25.0)
        assert score_cache.get("cadd-remote", "chr1", 100, "A", "T") == 25.0

    def test_different_sources_are_independent(
        self, score_cache: RemoteScoreCache
    ) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 10.0)
        score_cache.put("gnomad-remote", "chr1", 100, "A", "T", 0.001)
        assert score_cache.get("cadd-remote", "chr1", 100, "A", "T") == 10.0
        assert score_cache.get("gnomad-remote", "chr1", 100, "A", "T") == 0.001

    def test_get_batch_returns_positional_results(
        self, score_cache: RemoteScoreCache
    ) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        score_cache.put("cadd-remote", "chr1", 200, "G", "C", 15.0)

        variants = [
            ("chr1", 100, "A", "T"),
            ("chr1", 150, "C", "G"),  # not cached
            ("chr1", 200, "G", "C"),
        ]
        results = score_cache.get_batch("cadd-remote", variants)
        assert results == [23.5, None, 15.0]

    def test_put_batch(self, score_cache: RemoteScoreCache) -> None:
        entries = [
            ("chr1", 100, "A", "T", 23.5),
            ("chr1", 200, "G", "C", 15.0),
            ("chr2", 50, "T", "A", 8.2),
        ]
        score_cache.put_batch("cadd-remote", entries)
        assert score_cache.get("cadd-remote", "chr1", 100, "A", "T") == 23.5
        assert score_cache.get("cadd-remote", "chr2", 50, "T", "A") == 8.2

    def test_put_batch_empty_is_noop(self, score_cache: RemoteScoreCache) -> None:
        score_cache.put_batch("cadd-remote", [])
        stats = score_cache.stats()
        assert stats.entry_count == 0

    def test_expired_entry_returns_none(self, tmp_path: Path) -> None:
        # Use 1-day TTL and insert a record with timestamp far in the past
        cache = RemoteScoreCache(db_path=tmp_path / "short_ttl.db", ttl_days=1)
        # Manually insert with very old timestamp (epoch 0 = 1970)
        conn = cache._ensure_connection()
        conn.execute(
            "INSERT INTO remote_scores VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("cadd-remote", "chr1", 100, "A", "T", 23.5, 0),
        )
        conn.commit()
        assert cache.get("cadd-remote", "chr1", 100, "A", "T") is None
        cache.close()

    def test_pinned_mode_never_expires(
        self, pinned_score_cache: RemoteScoreCache
    ) -> None:
        # Insert with a very old timestamp
        conn = pinned_score_cache._ensure_connection()
        conn.execute(
            "INSERT INTO remote_scores VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("cadd-remote", "chr1", 100, "A", "T", 23.5, 1000000),
        )
        conn.commit()
        assert pinned_score_cache.get("cadd-remote", "chr1", 100, "A", "T") == 23.5

    def test_evict_expired(self, tmp_path: Path) -> None:
        cache = RemoteScoreCache(db_path=tmp_path / "evict.db", ttl_days=1)
        # Insert with current time
        cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        # Insert with old time
        conn = cache._ensure_connection()
        conn.execute(
            "INSERT OR REPLACE INTO remote_scores VALUES (?, ?, ?, ?, ?, ?, ?)",
            ("cadd-remote", "chr1", 200, "G", "C", 15.0, 0),
        )
        conn.commit()
        evicted = cache.evict_expired()
        assert evicted == 1
        assert cache.get("cadd-remote", "chr1", 100, "A", "T") == 23.5
        assert cache.get("cadd-remote", "chr1", 200, "G", "C") is None
        cache.close()

    def test_evict_expired_pinned_is_noop(
        self, pinned_score_cache: RemoteScoreCache
    ) -> None:
        pinned_score_cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        assert pinned_score_cache.evict_expired() == 0

    def test_clear_all(self, score_cache: RemoteScoreCache) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        score_cache.put("gnomad-remote", "chr1", 100, "A", "T", 0.001)
        deleted = score_cache.clear()
        assert deleted == 2
        assert score_cache.stats().entry_count == 0

    def test_clear_by_source(self, score_cache: RemoteScoreCache) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        score_cache.put("gnomad-remote", "chr1", 100, "A", "T", 0.001)
        deleted = score_cache.clear(source="cadd-remote")
        assert deleted == 1
        assert score_cache.get("gnomad-remote", "chr1", 100, "A", "T") == 0.001

    def test_stats(self, score_cache: RemoteScoreCache) -> None:
        score_cache.put("cadd-remote", "chr1", 100, "A", "T", 23.5)
        score_cache.put("gnomad-remote", "chr1", 100, "A", "T", 0.001)
        stats = score_cache.stats()
        assert stats.entry_count == 2
        assert stats.entries_by_source == {"cadd-remote": 1, "gnomad-remote": 1}
        assert stats.disk_size_bytes > 0


# ============================================================
# CircuitBreaker tests
# ============================================================


class TestCircuitBreaker:
    """Circuit breaker state machine transitions."""

    def test_starts_closed(self) -> None:
        cb = CircuitBreaker()
        assert cb.state == CircuitState.CLOSED
        assert not cb.is_open

    def test_single_failure_stays_closed(self) -> None:
        cb = CircuitBreaker(failure_threshold=5)
        cb.record_failure()
        assert cb.state == CircuitState.CLOSED

    def test_opens_at_threshold(self) -> None:
        cb = CircuitBreaker(failure_threshold=3, failure_window_seconds=60.0)
        cb.record_failure()
        cb.record_failure()
        cb.record_failure()
        assert cb.state == CircuitState.OPEN
        assert cb.is_open

    def test_success_resets_failure_count(self) -> None:
        cb = CircuitBreaker(failure_threshold=3)
        cb.record_failure()
        cb.record_failure()
        cb.record_success()
        cb.record_failure()
        cb.record_failure()
        # Still needs 3 consecutive — success reset the counter
        assert cb.state == CircuitState.CLOSED

    def test_transitions_to_half_open_after_recovery(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=2,
            recovery_seconds=0.01,
        )
        cb.record_failure()
        cb.record_failure()
        assert cb.state == CircuitState.OPEN

        time.sleep(0.02)
        assert cb.state == CircuitState.HALF_OPEN
        assert not cb.is_open  # half-open allows probe

    def test_success_in_half_open_closes_circuit(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=2,
            recovery_seconds=0.01,
        )
        cb.record_failure()
        cb.record_failure()
        time.sleep(0.02)
        assert cb.state == CircuitState.HALF_OPEN

        cb.record_success()
        assert cb.state == CircuitState.CLOSED

    def test_failure_in_half_open_reopens(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=2,
            recovery_seconds=0.01,
        )
        cb.record_failure()
        cb.record_failure()
        time.sleep(0.02)
        assert cb.state == CircuitState.HALF_OPEN

        cb.record_failure()
        assert cb.state == CircuitState.OPEN

    def test_old_failures_outside_window_are_forgotten(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=3,
            failure_window_seconds=0.01,
        )
        cb.record_failure()
        cb.record_failure()
        time.sleep(0.02)
        # Old failures expired
        cb.record_failure()
        assert cb.state == CircuitState.CLOSED


# ============================================================
# RemoteTabixCADD tests (mocked pysam)
# ============================================================


class TestRemoteTabixCADD:
    """CADD backend with mocked pysam.TabixFile."""

    def _make_config(self, tmp_path: Path) -> RemoteTabixConfig:
        return RemoteTabixConfig(
            cadd_remote_url="cadd-v1.7-grch38",
            cache_path=tmp_path / "test_cadd_cache.db",
            cache_ttl_days=30,
        )

    @patch("vartriage.remote.cadd.pysam.TabixFile")
    def test_lookup_batch_parses_cadd_records(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.cadd import RemoteTabixCADD

        config = self._make_config(tmp_path)

        # Mock tabix returns CADD-format records
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter(
            [
                "22\t100\tA\tT\t3.14\t23.5",
                "22\t100\tA\tG\t2.10\t15.0",
                "22\t100\tA\tC\t1.50\t10.0",
            ]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixCADD(config)
        variants = [("chr22", 100, "A", "T")]
        results = backend.lookup_batch(variants)

        assert results == {("chr22", 100, "A", "T"): 23.5}
        backend.close()

    @patch("vartriage.remote.cadd.pysam.TabixFile")
    def test_lookup_batch_strips_chr_prefix(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.cadd import RemoteTabixCADD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter([])
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixCADD(config)
        backend.lookup_batch([("chr22", 100, "A", "T")])

        # Should query with bare "22", not "chr22"
        mock_tabix.fetch.assert_called_once_with("22", 99, 100)
        backend.close()

    @patch("vartriage.remote.cadd.pysam.TabixFile")
    def test_cache_hit_skips_network(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.cadd import RemoteTabixCADD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter(["22\t100\tA\tT\t3.14\t23.5"])
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixCADD(config)

        # First call: network fetch
        results1 = backend.lookup_batch([("chr22", 100, "A", "T")])
        assert results1 == {("chr22", 100, "A", "T"): 23.5}
        assert backend.network_fetches == 1

        # Second call: cache hit
        mock_tabix.fetch.return_value = iter([])
        results2 = backend.lookup_batch([("chr22", 100, "A", "T")])
        assert results2 == {("chr22", 100, "A", "T"): 23.5}
        assert backend.cache_hits == 1
        backend.close()

    @patch("vartriage.remote.cadd.pysam.TabixFile")
    def test_network_failure_returns_empty_gracefully(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.cadd import RemoteTabixCADD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.side_effect = OSError("Connection refused")
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixCADD(config)
        results = backend.lookup_batch([("chr22", 100, "A", "T")])

        # No crash, just empty results
        assert results == {}
        backend.close()

    @patch("vartriage.remote.cadd.pysam.TabixFile")
    def test_circuit_breaker_stops_queries(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.cadd import RemoteTabixCADD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.side_effect = OSError("timeout")
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixCADD(config)

        # Trip the circuit breaker (5 failures)
        for _ in range(5):
            backend.lookup_batch([("chr22", i * 20000, "A", "T") for i in range(1)])

        # Breaker open — no more fetch calls
        call_count_before = mock_tabix.fetch.call_count
        backend.lookup_batch([("chr22", 999999, "G", "C")])
        assert mock_tabix.fetch.call_count == call_count_before
        backend.close()

    @patch("vartriage.remote.cadd.pysam.TabixFile")
    def test_batch_window_groups_nearby_variants(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.cadd import RemoteTabixCADD

        config = RemoteTabixConfig(
            cadd_remote_url="cadd-v1.7-grch38",
            cache_path=tmp_path / "batch_test.db",
            batch_window_bp=10_000,
        )
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter(
            [
                "22\t100\tA\tT\t3.14\t23.5",
                "22\t200\tG\tC\t2.50\t18.0",
            ]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixCADD(config)
        # Two variants within 10kb should produce one fetch call
        variants = [("chr22", 100, "A", "T"), ("chr22", 200, "G", "C")]
        backend.lookup_batch(variants)

        assert mock_tabix.fetch.call_count == 1
        backend.close()


# ============================================================
# RemoteTabixGnomAD tests (mocked pysam)
# ============================================================


class TestRemoteTabixGnomAD:
    """gnomAD backend with mocked pysam.TabixFile."""

    def _make_config(self, tmp_path: Path) -> RemoteTabixConfig:
        return RemoteTabixConfig(
            gnomad_remote_url="https://example.com/{chrom}.vcf.bgz",
            cache_path=tmp_path / "test_gnomad_cache.db",
            cache_ttl_days=30,
        )

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_lookup_batch_parses_gnomad_vcf(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        # gnomAD VCF record: CHROM POS ID REF ALT QUAL FILTER INFO
        mock_tabix.fetch.return_value = iter(
            [
                "chr22\t100\t.\tA\tT\t.\tPASS\tAF=0.0015;AN=100000",
            ]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        results = backend.lookup_batch([("chr22", 100, "A", "T")])

        assert results == [0.0015]
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_handles_multiallelic_vcf(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        # Multi-allelic: A>T,G with two AF values
        mock_tabix.fetch.return_value = iter(
            [
                "chr22\t100\t.\tA\tT,G\t.\tPASS\tAF=0.0015,0.05;AN=100000",
            ]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        results = backend.lookup_batch([("chr22", 100, "A", "G")])

        assert results == [0.05]
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_variant_not_found_returns_none(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter([])
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        results = backend.lookup_batch([("chr22", 100, "A", "T")])

        assert results == [None]
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_cache_hit_skips_network(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter(
            ["chr22\t100\t.\tA\tT\t.\tPASS\tAF=0.0015;AN=100000"]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)

        # First call: network
        results1 = backend.lookup_batch([("chr22", 100, "A", "T")])
        assert results1 == [0.0015]
        assert backend.network_fetches == 1

        # Second call: cache hit
        mock_tabix.fetch.return_value = iter([])
        results2 = backend.lookup_batch([("chr22", 100, "A", "T")])
        assert results2 == [0.0015]
        assert backend.cache_hits == 1
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_network_failure_returns_none_gracefully(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.side_effect = OSError("Connection refused")
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        results = backend.lookup_batch([("chr22", 100, "A", "T")])

        assert results == [None]
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_per_chromosome_url_formatting(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = self._make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter([])
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        backend.lookup_batch([("chr22", 100, "A", "T")])

        # Should open with {chrom} replaced
        mock_tabix_cls.assert_called_with("https://example.com/chr22.vcf.bgz")
        backend.close()


# ============================================================
# CLI integration tests
# ============================================================


class TestCLIRemoteFlags:
    """CLI flag parsing for remote tabix arguments."""

    def test_build_remote_config_returns_none_when_no_flags(self) -> None:
        from vartriage.cli import _build_remote_config

        args = MagicMock()
        args.cadd_remote = None
        args.gnomad_remote = None
        args.remote_cache_ttl = 30
        args.cadd_scores = None
        args.gnomad = None

        assert _build_remote_config(args) is None

    def test_build_remote_config_with_cadd_preset(self) -> None:
        from vartriage.cli import _build_remote_config

        args = MagicMock()
        args.cadd_remote = "cadd-v1.7-grch38"
        args.gnomad_remote = None
        args.remote_cache_ttl = 60
        args.cadd_scores = None
        args.gnomad = None

        config = _build_remote_config(args)
        assert config is not None
        assert config.cadd_remote_url == "cadd-v1.7-grch38"
        assert config.cache_ttl_days == 60

    def test_local_cadd_overrides_remote(self) -> None:
        from vartriage.cli import _build_remote_config

        args = MagicMock()
        args.cadd_remote = "cadd-v1.7-grch38"
        args.gnomad_remote = "gnomad-exomes-v4-grch38"
        args.remote_cache_ttl = 30
        args.cadd_scores = Path("/some/local.tsv")  # local wins
        args.gnomad = None

        config = _build_remote_config(args)
        assert config is not None
        assert config.cadd_remote_url is None  # overridden by local
        assert config.gnomad_remote_url == "gnomad-exomes-v4-grch38"

    def test_local_gnomad_overrides_remote(self) -> None:
        from vartriage.cli import _build_remote_config

        args = MagicMock()
        args.cadd_remote = "cadd-v1.7-grch38"
        args.gnomad_remote = "gnomad-exomes-v4-grch38"
        args.remote_cache_ttl = 30
        args.cadd_scores = None
        args.gnomad = Path("/some/gnomad.vcf.bgz")  # local wins

        config = _build_remote_config(args)
        assert config is not None
        assert config.cadd_remote_url == "cadd-v1.7-grch38"
        assert config.gnomad_remote_url is None  # overridden by local

    def test_both_local_override_returns_none(self) -> None:
        from vartriage.cli import _build_remote_config

        args = MagicMock()
        args.cadd_remote = "cadd-v1.7-grch38"
        args.gnomad_remote = "gnomad-exomes-v4-grch38"
        args.remote_cache_ttl = 30
        args.cadd_scores = Path("/some/local.tsv")
        args.gnomad = Path("/some/gnomad.vcf.bgz")

        config = _build_remote_config(args)
        assert config is None  # both overridden, no remote needed
