"""Unit tests for the circuit breaker state machine."""

from __future__ import annotations

import time

import pytest

from vartriage.api._circuit_breaker import (CircuitBreaker, CircuitBreakerOpen,
                                            CircuitState)


class TestClosedState:
    """Normal operation (CLOSED state)."""

    def test_starts_in_closed_state(self) -> None:
        cb = CircuitBreaker(service_name="test")
        assert cb.state == CircuitState.CLOSED

    def test_allows_requests_when_closed(self) -> None:
        cb = CircuitBreaker(service_name="test")
        assert cb.allow_request() is True

    def test_stays_closed_below_failure_threshold(self) -> None:
        cb = CircuitBreaker(failure_threshold=5, service_name="test")
        for _ in range(4):
            cb.record_failure()
        assert cb.state == CircuitState.CLOSED
        assert cb.failure_count == 4

    def test_success_resets_failure_count(self) -> None:
        cb = CircuitBreaker(failure_threshold=5, service_name="test")
        cb.record_failure()
        cb.record_failure()
        cb.record_success()
        assert cb.failure_count == 0
        assert cb.state == CircuitState.CLOSED


class TestOpenState:
    """Failure isolation (OPEN state)."""

    def test_trips_open_at_failure_threshold(self) -> None:
        cb = CircuitBreaker(failure_threshold=3, service_name="test")
        for _ in range(3):
            cb.record_failure()
        assert cb.state == CircuitState.OPEN

    def test_rejects_requests_when_open(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=2, recovery_timeout=60.0, service_name="vep"
        )
        cb.record_failure()
        cb.record_failure()

        with pytest.raises(CircuitBreakerOpen) as exc_info:
            cb.allow_request()

        assert exc_info.value.service == "vep"
        assert exc_info.value.seconds_until_probe > 0

    def test_is_open_property(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=60.0, service_name="test"
        )
        assert cb.is_open is False
        cb.record_failure()
        assert cb.is_open is True


class TestHalfOpenState:
    """Recovery probing (HALF_OPEN state)."""

    def test_transitions_to_half_open_after_recovery_timeout(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=0.05, service_name="test"
        )
        cb.record_failure()
        assert cb.state == CircuitState.OPEN

        time.sleep(0.06)
        assert cb.state == CircuitState.HALF_OPEN

    def test_half_open_allows_single_probe_request(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=0.05, service_name="test"
        )
        cb.record_failure()
        time.sleep(0.06)

        assert cb.allow_request() is True

    def test_probe_success_closes_circuit(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=0.05, service_name="test"
        )
        cb.record_failure()
        time.sleep(0.06)

        cb.record_success()
        assert cb.state == CircuitState.CLOSED
        assert cb.failure_count == 0

    def test_probe_failure_reopens_circuit(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=0.05, service_name="test"
        )
        cb.record_failure()
        time.sleep(0.06)

        # Probe fails
        cb.record_failure()
        assert cb.state == CircuitState.OPEN


class TestReset:
    """Manual reset for testing and recovery."""

    def test_reset_returns_to_closed(self) -> None:
        cb = CircuitBreaker(failure_threshold=1, service_name="test")
        cb.record_failure()
        assert cb.state == CircuitState.OPEN

        cb.reset()
        assert cb.state == CircuitState.CLOSED
        assert cb.failure_count == 0

    def test_reset_from_half_open(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=0.01, service_name="test"
        )
        cb.record_failure()
        time.sleep(0.02)
        assert cb.state == CircuitState.HALF_OPEN

        cb.reset()
        assert cb.state == CircuitState.CLOSED


class TestValidation:
    """Constructor validation."""

    def test_rejects_zero_failure_threshold(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            CircuitBreaker(failure_threshold=0)

    def test_rejects_negative_recovery_timeout(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            CircuitBreaker(recovery_timeout=-1.0)

    def test_rejects_zero_recovery_timeout(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            CircuitBreaker(recovery_timeout=0.0)


class TestErrorMessages:
    """Exception content for debugging."""

    def test_open_exception_includes_time_until_probe(self) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=10.0, service_name="spliceai"
        )
        cb.record_failure()

        with pytest.raises(CircuitBreakerOpen) as exc_info:
            cb.allow_request()

        # Should be close to 10 seconds (just opened)
        assert exc_info.value.seconds_until_probe > 9.0
        assert "spliceai" in str(exc_info.value)
