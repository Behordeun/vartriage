"""Circuit breaker for external API call protection.

Prevents cascading failures by tracking consecutive errors per service.
After a threshold of failures, the breaker opens and rejects requests
for a recovery period before allowing a probe request through.

State machine:
    CLOSED  -> (failure_threshold reached) -> OPEN
    OPEN    -> (recovery_timeout elapsed)  -> HALF_OPEN
    HALF_OPEN -> (probe succeeds)          -> CLOSED
    HALF_OPEN -> (probe fails)             -> OPEN
"""

from __future__ import annotations

import threading
import time
from enum import Enum


class CircuitState(Enum):
    """Circuit breaker states."""

    CLOSED = "closed"
    OPEN = "open"
    HALF_OPEN = "half_open"


class CircuitBreakerOpen(Exception):
    """Raised when a request is rejected because the circuit is open."""

    def __init__(self, service: str, seconds_until_probe: float) -> None:
        self.service = service
        self.seconds_until_probe = seconds_until_probe
        super().__init__(
            f"Circuit breaker open for {service}. "
            f"Next probe in {seconds_until_probe:.1f}s."
        )


class CircuitBreaker:
    """Three-state circuit breaker for external service calls.

    Thread-safe. Callers check `allow_request()` before making a call,
    then report the outcome via `record_success()` or `record_failure()`.

    Parameters
    ----------
    failure_threshold
        Consecutive failures needed to trip the breaker open.
    recovery_timeout
        Seconds the breaker stays open before allowing a probe.
    service_name
        Identifier for error messages and logging.
    """

    def __init__(
        self,
        failure_threshold: int = 5,
        recovery_timeout: float = 60.0,
        service_name: str = "unknown",
    ) -> None:
        if failure_threshold < 1:
            raise ValueError("failure_threshold must be at least 1")
        if recovery_timeout <= 0:
            raise ValueError("recovery_timeout must be positive")

        self._failure_threshold = failure_threshold
        self._recovery_timeout = recovery_timeout
        self._service_name = service_name

        self._state = CircuitState.CLOSED
        self._consecutive_failures = 0
        self._opened_at: float = 0.0
        self._lock = threading.Lock()

    @property
    def state(self) -> CircuitState:
        """Current circuit state (may transition OPEN -> HALF_OPEN on read)."""
        with self._lock:
            self._maybe_transition_to_half_open()
            return self._state

    @property
    def is_open(self) -> bool:
        """True if the circuit is in OPEN state (rejecting requests)."""
        return self.state == CircuitState.OPEN

    @property
    def failure_count(self) -> int:
        """Current consecutive failure count."""
        with self._lock:
            return self._consecutive_failures

    def allow_request(self) -> bool:
        """Check whether a request should proceed.

        Returns True for CLOSED (normal operation) and HALF_OPEN (probe).

        Raises
        ------
        CircuitBreakerOpen
            When the circuit is open and the recovery period hasn't elapsed.
        """
        with self._lock:
            self._maybe_transition_to_half_open()

            if self._state == CircuitState.OPEN:
                elapsed = time.monotonic() - self._opened_at
                remaining = self._recovery_timeout - elapsed
                raise CircuitBreakerOpen(self._service_name, max(0.0, remaining))

            return True

    def record_success(self) -> None:
        """Record a successful call. Resets failure count, closes circuit."""
        with self._lock:
            self._consecutive_failures = 0
            self._maybe_transition_to_half_open()
            if self._state == CircuitState.HALF_OPEN:
                self._state = CircuitState.CLOSED

    def record_failure(self) -> None:
        """Record a failed call. May trip the breaker open."""
        with self._lock:
            self._consecutive_failures += 1

            if self._state == CircuitState.HALF_OPEN:
                # Probe failed, go back to OPEN
                self._state = CircuitState.OPEN
                self._opened_at = time.monotonic()
                return

            if self._consecutive_failures >= self._failure_threshold:
                self._state = CircuitState.OPEN
                self._opened_at = time.monotonic()

    def reset(self) -> None:
        """Force-reset the breaker to CLOSED state. For testing and manual recovery."""
        with self._lock:
            self._state = CircuitState.CLOSED
            self._consecutive_failures = 0
            self._opened_at = 0.0

    def _maybe_transition_to_half_open(self) -> None:
        """Move from OPEN to HALF_OPEN if recovery timeout elapsed. Caller holds lock."""
        if self._state != CircuitState.OPEN:
            return
        elapsed = time.monotonic() - self._opened_at
        if elapsed >= self._recovery_timeout:
            self._state = CircuitState.HALF_OPEN
