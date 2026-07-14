"""Token bucket rate limiter for API request throttling.

Synchronous implementation with no asyncio dependency. Each service
gets its own RateLimiter instance configured with the appropriate
requests-per-second limit and optional daily cap.
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from datetime import date
from datetime import datetime as dt
from datetime import timezone


class RateLimitExhausted(Exception):
    """Raised when tokens cannot be acquired within the wait timeout."""

    def __init__(self, service: str, wait_seconds: float) -> None:
        self.service = service
        self.wait_seconds = wait_seconds
        super().__init__(
            f"Rate limit exhausted for {service}. "
            f"Next token available in {wait_seconds:.1f}s."
        )


class DailyLimitExhausted(Exception):
    """Raised when a service's daily request cap has been reached."""

    def __init__(self, service: str, limit: int, resets_at: str) -> None:
        self.service = service
        self.limit = limit
        self.resets_at = resets_at
        super().__init__(
            f"Daily limit ({limit}) exhausted for {service}. " f"Resets at {resets_at}."
        )


@dataclass
class RateLimiterStats:
    """Snapshot of rate limiter state for observability."""

    tokens_available: float
    tokens_per_second: float
    burst_capacity: int
    daily_count: int
    daily_limit: int | None
    daily_remaining: int | None


class RateLimiter:
    """Token bucket rate limiter with optional daily cap.

    Thread-safe. Callers block on `acquire()` until a token becomes
    available or the wait timeout is exceeded.

    Parameters
    ----------
    tokens_per_second
        Sustained refill rate.
    burst
        Maximum bucket capacity. Defaults to ceil(tokens_per_second).
    daily_limit
        Maximum requests per UTC day. None means unlimited.
    service_name
        Identifier for error messages and logging.
    """

    def __init__(
        self,
        tokens_per_second: float,
        burst: int | None = None,
        daily_limit: int | None = None,
        service_name: str = "unknown",
    ) -> None:
        if tokens_per_second <= 0:
            raise ValueError("tokens_per_second must be positive")

        self._rate = tokens_per_second
        self._burst = burst if burst is not None else max(1, int(tokens_per_second) + 1)
        self._daily_limit = daily_limit
        self._service_name = service_name

        self._tokens = float(self._burst)
        self._last_refill = time.monotonic()
        self._lock = threading.Lock()

        # Daily tracking resets at UTC midnight
        self._daily_count = 0
        self._daily_reset_date = dt.now(timezone.utc).date()

    def acquire(self, tokens: int = 1, timeout: float = 60.0) -> None:
        """Block until the requested number of tokens are available.

        Parameters
        ----------
        tokens
            Number of tokens to consume (typically 1 per request).
        timeout
            Maximum seconds to wait. Raises RateLimitExhausted if exceeded.

        Raises
        ------
        RateLimitExhausted
            If tokens cannot be acquired within the timeout.
        DailyLimitExhausted
            If the daily request cap has been reached.
        """
        deadline = time.monotonic() + timeout

        while True:
            with self._lock:
                self._maybe_reset_daily()
                self._refill()

                if (
                    self._daily_limit is not None
                    and self._daily_count >= self._daily_limit
                ):
                    tomorrow = dt.now(timezone.utc).date().isoformat()
                    raise DailyLimitExhausted(
                        self._service_name, self._daily_limit, f"{tomorrow} UTC"
                    )

                if self._tokens >= tokens:
                    self._tokens -= tokens
                    self._daily_count += tokens
                    return

                # How long until enough tokens accumulate?
                deficit = tokens - self._tokens
                wait_needed = deficit / self._rate

            remaining = deadline - time.monotonic()
            if remaining <= 0 or wait_needed > remaining:
                raise RateLimitExhausted(self._service_name, wait_needed)

            # Sleep for the shorter of: time to next token, or remaining timeout
            sleep_time = min(wait_needed, remaining)
            time.sleep(sleep_time)

    def remaining(self) -> float:
        """Current available token count (approximate, for observability)."""
        with self._lock:
            self._refill()
            return self._tokens

    def daily_remaining(self) -> int | None:
        """Requests remaining today, or None if no daily limit."""
        if self._daily_limit is None:
            return None
        with self._lock:
            self._maybe_reset_daily()
            return max(0, self._daily_limit - self._daily_count)

    def stats(self) -> RateLimiterStats:
        """Snapshot of current limiter state."""
        with self._lock:
            self._maybe_reset_daily()
            self._refill()
            daily_remaining = (
                max(0, self._daily_limit - self._daily_count)
                if self._daily_limit is not None
                else None
            )
            return RateLimiterStats(
                tokens_available=self._tokens,
                tokens_per_second=self._rate,
                burst_capacity=self._burst,
                daily_count=self._daily_count,
                daily_limit=self._daily_limit,
                daily_remaining=daily_remaining,
            )

    def _refill(self) -> None:
        """Add tokens based on elapsed time since last refill. Caller holds lock."""
        now = time.monotonic()
        elapsed = now - self._last_refill
        if elapsed <= 0:
            return
        self._tokens = min(self._burst, self._tokens + elapsed * self._rate)
        self._last_refill = now

    def _maybe_reset_daily(self) -> None:
        """Reset daily counter if we've crossed UTC midnight. Caller holds lock."""
        today = dt.now(timezone.utc).date()
        if today > self._daily_reset_date:
            self._daily_count = 0
            self._daily_reset_date = today
