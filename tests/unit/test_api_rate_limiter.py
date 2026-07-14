"""Unit tests for the token bucket rate limiter."""

from __future__ import annotations

import threading
import time
from unittest.mock import patch

import pytest

from vartriage.api._rate_limiter import (DailyLimitExhausted, RateLimiter,
                                         RateLimitExhausted)


class TestTokenBucketBasics:
    """Core token bucket behavior."""

    def test_acquire_single_token_succeeds_immediately(self) -> None:
        limiter = RateLimiter(tokens_per_second=10.0, service_name="test")
        limiter.acquire()  # should not raise

    def test_burst_capacity_allows_multiple_immediate_acquires(self) -> None:
        limiter = RateLimiter(tokens_per_second=5.0, burst=5, service_name="test")
        for _ in range(5):
            limiter.acquire()

    def test_exceeding_burst_blocks_until_refill(self) -> None:
        limiter = RateLimiter(tokens_per_second=100.0, burst=2, service_name="test")
        limiter.acquire()
        limiter.acquire()
        # Third acquire needs a refill (at 100/sec, ~10ms for one token)
        start = time.monotonic()
        limiter.acquire()
        elapsed = time.monotonic() - start
        assert elapsed < 0.1  # should refill fast at 100/sec

    def test_remaining_reports_available_tokens(self) -> None:
        limiter = RateLimiter(tokens_per_second=10.0, burst=10, service_name="test")
        initial = limiter.remaining()
        assert initial == pytest.approx(10.0, abs=0.5)
        limiter.acquire(3)
        after = limiter.remaining()
        assert after == pytest.approx(7.0, abs=0.5)

    def test_tokens_refill_over_time(self) -> None:
        limiter = RateLimiter(tokens_per_second=100.0, burst=5, service_name="test")
        limiter.acquire(5)
        time.sleep(0.03)  # 30ms at 100/sec = ~3 tokens
        available = limiter.remaining()
        assert available >= 2.0

    def test_tokens_never_exceed_burst_capacity(self) -> None:
        limiter = RateLimiter(tokens_per_second=1000.0, burst=3, service_name="test")
        time.sleep(0.05)  # would refill 50 tokens, but burst caps at 3
        available = limiter.remaining()
        assert available <= 3.0


class TestRateLimitExhausted:
    """Timeout behavior when tokens unavailable."""

    def test_raises_after_timeout(self) -> None:
        limiter = RateLimiter(tokens_per_second=1.0, burst=1, service_name="vep")
        limiter.acquire()  # drain the bucket

        with pytest.raises(RateLimitExhausted) as exc_info:
            limiter.acquire(timeout=0.05)

        assert exc_info.value.service == "vep"
        assert exc_info.value.wait_seconds > 0

    def test_error_message_includes_service_name(self) -> None:
        limiter = RateLimiter(tokens_per_second=1.0, burst=1, service_name="clinvar")
        limiter.acquire()

        with pytest.raises(RateLimitExhausted, match="clinvar"):
            limiter.acquire(timeout=0.01)


class TestDailyLimit:
    """Daily request cap behavior."""

    def test_daily_limit_exhaustion_raises(self) -> None:
        limiter = RateLimiter(
            tokens_per_second=1000.0, burst=100, daily_limit=5, service_name="vep"
        )
        for _ in range(5):
            limiter.acquire()

        with pytest.raises(DailyLimitExhausted) as exc_info:
            limiter.acquire()

        assert exc_info.value.service == "vep"
        assert exc_info.value.limit == 5

    def test_daily_remaining_tracks_usage(self) -> None:
        limiter = RateLimiter(
            tokens_per_second=1000.0, burst=100, daily_limit=10, service_name="test"
        )
        assert limiter.daily_remaining() == 10
        limiter.acquire(3)
        assert limiter.daily_remaining() == 7

    def test_daily_remaining_returns_none_when_no_limit(self) -> None:
        limiter = RateLimiter(tokens_per_second=10.0, service_name="test")
        assert limiter.daily_remaining() is None

    def test_daily_count_resets_on_new_utc_day(self) -> None:
        limiter = RateLimiter(
            tokens_per_second=1000.0, burst=100, daily_limit=5, service_name="test"
        )
        for _ in range(5):
            limiter.acquire()

        # Simulate crossing UTC midnight by backdating the reset date
        from datetime import datetime as dt
        from datetime import timedelta, timezone

        yesterday_utc = dt.now(timezone.utc).date() - timedelta(days=1)
        with limiter._lock:
            limiter._daily_reset_date = yesterday_utc

        # Should reset and allow new requests
        limiter.acquire()
        assert limiter.daily_remaining() == 4


class TestStats:
    """Observability via stats()."""

    def test_stats_reflects_current_state(self) -> None:
        limiter = RateLimiter(
            tokens_per_second=10.0, burst=10, daily_limit=100, service_name="test"
        )
        limiter.acquire(3)
        stats = limiter.stats()

        assert stats.tokens_per_second == 10.0
        assert stats.burst_capacity == 10
        assert stats.daily_limit == 100
        assert stats.daily_count == 3
        assert stats.daily_remaining == 97
        assert stats.tokens_available == pytest.approx(7.0, abs=1.0)

    def test_stats_without_daily_limit(self) -> None:
        limiter = RateLimiter(tokens_per_second=5.0, service_name="test")
        stats = limiter.stats()
        assert stats.daily_limit is None
        assert stats.daily_remaining is None


class TestValidation:
    """Constructor validation."""

    def test_rejects_zero_rate(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            RateLimiter(tokens_per_second=0.0)

    def test_rejects_negative_rate(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            RateLimiter(tokens_per_second=-1.0)


class TestThreadSafety:
    """Concurrent access from multiple threads."""

    def test_concurrent_acquires_respect_burst(self) -> None:
        limiter = RateLimiter(
            tokens_per_second=1000.0, burst=10, daily_limit=100, service_name="test"
        )
        results: list[bool] = []
        errors: list[Exception] = []

        def worker() -> None:
            try:
                limiter.acquire(timeout=1.0)
                results.append(True)
            except (RateLimitExhausted, DailyLimitExhausted) as e:
                errors.append(e)

        threads = [threading.Thread(target=worker) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5.0)

        # All 10 should succeed (burst=10, rate=1000/sec)
        assert len(results) == 10
        assert len(errors) == 0
