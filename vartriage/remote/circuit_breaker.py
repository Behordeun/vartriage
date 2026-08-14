"""Circuit breaker for remote tabix connections.

Tracks consecutive failures and opens the circuit when a threshold is
reached, preventing further network attempts for the remainder of the
run (or until the half-open recovery window passes).

State machine:
    CLOSED  -> failures >= threshold within window -> OPEN
    OPEN    -> recovery_seconds elapsed            -> HALF_OPEN
    HALF_OPEN -> next call succeeds                -> CLOSED
    HALF_OPEN -> next call fails                   -> OPEN
"""

from __future__ import annotations

import logging
import time
from enum import Enum

logger = logging.getLogger(__name__)


class CircuitState(Enum):
    """Circuit breaker states."""

    CLOSED = "closed"
    OPEN = "open"
    HALF_OPEN = "half_open"


class CircuitBreaker:
    """Failure-tracking circuit breaker for remote connections.

    Opens after `failure_threshold` consecutive failures within
    `failure_window_seconds`. Transitions to half-open after
    `recovery_seconds`, allowing a single probe request.

    Parameters
    ----------
    failure_threshold : int
        Number of consecutive failures to trigger open state.
    failure_window_seconds : float
        Time window for counting consecutive failures. Failures
        older than this are forgotten.
    recovery_seconds : float
        Time to wait in open state before allowing a probe.
    """

    def __init__(
        self,
        failure_threshold: int = 5,
        failure_window_seconds: float = 60.0,
        recovery_seconds: float = 30.0,
    ) -> None:
        self._failure_threshold = failure_threshold
        self._failure_window = failure_window_seconds
        self._recovery_seconds = recovery_seconds

        self._state = CircuitState.CLOSED
        self._failure_count: int = 0
        self._first_failure_time: float | None = None
        self._opened_at: float = 0.0

    @property
    def state(self) -> CircuitState:
        """Current circuit state, accounting for recovery timeout."""
        self._check_recovery()
        return self._state

    @property
    def is_open(self) -> bool:
        """True when the circuit is open (no requests should be made).

        Returns False for both CLOSED and HALF_OPEN states, since
        half-open allows a probe request.
        """
        self._check_recovery()
        return self._state == CircuitState.OPEN

    def record_success(self) -> None:
        """Record a successful operation.

        Resets the failure counter. If in HALF_OPEN state,
        transitions back to CLOSED.
        """
        if self._state == CircuitState.HALF_OPEN:
            logger.info("Circuit breaker closing after successful probe")
            self._state = CircuitState.CLOSED

        self._failure_count = 0
        self._first_failure_time = None

    def record_failure(self) -> None:
        """Record a failed operation.

        If in HALF_OPEN state, immediately re-opens the circuit.
        In CLOSED state, adds to the failure count and opens if
        the threshold is reached within the window.
        """
        now = time.monotonic()

        if self._state == CircuitState.HALF_OPEN:
            self._open_circuit(now)
            return

        if self._first_failure_time is None:
            self._first_failure_time = now
            self._failure_count = 1
        elif now - self._first_failure_time > self._failure_window:
            # Window expired — start a new window
            self._first_failure_time = now
            self._failure_count = 1
        else:
            self._failure_count += 1

        if self._failure_count >= self._failure_threshold:
            self._open_circuit(now)

    def _check_recovery(self) -> None:
        """Transition from OPEN to HALF_OPEN if recovery period elapsed."""
        if self._state == CircuitState.OPEN:
            elapsed = time.monotonic() - self._opened_at
            if elapsed >= self._recovery_seconds:
                self._state = CircuitState.HALF_OPEN
                logger.info(
                    "Circuit breaker transitioning to HALF_OPEN "
                    "after %.1fs recovery period",
                    elapsed,
                )

    def _open_circuit(self, now: float) -> None:
        """Transition to OPEN state."""
        self._state = CircuitState.OPEN
        self._opened_at = now
        self._failure_count = 0
        self._first_failure_time = None
        logger.warning(
            "Circuit breaker OPEN — remote queries disabled for %.0fs",
            self._recovery_seconds,
        )
