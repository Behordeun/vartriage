"""A systematic client error trips the circuit breaker.

A 4xx that is not retryable (schema drift, a revoked key, a moved endpoint)
used to be recorded as a circuit-breaker success, so a total upstream
misconfiguration looked healthy while every call failed. A non-retryable
4xx now records a failure, so a run of them opens the breaker.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.api._base import APIClientError, BaseAPIClient
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache  # noqa: E402


def _mock_transport(status: int):  # noqa: ANN202
    def handler(request):  # noqa: ANN001, ANN202
        return httpx.Response(status, text="error")

    return httpx.MockTransport(handler)


def _client(status: int, breaker: CircuitBreaker, tmp_path: Path) -> BaseAPIClient:
    client = BaseAPIClient(
        base_url="https://test.example.com",
        rate_limiter=RateLimiter(
            tokens_per_second=1000.0, burst=100, service_name="test"
        ),
        cache=ResponseCache(db_path=tmp_path / "cache.db", default_ttl_days=1),
        circuit_breaker=breaker,
        service_name="test",
        max_retries=1,
    )
    client._client = httpx.Client(
        base_url="https://test.example.com",
        transport=_mock_transport(status),
        timeout=httpx.Timeout(5.0),
    )
    return client


class TestFourxxTripsBreaker:
    def test_repeated_403_open_the_breaker(self, tmp_path: Path) -> None:
        breaker = CircuitBreaker(failure_threshold=3, service_name="test")
        client = _client(403, breaker, tmp_path)

        for _ in range(3):
            with pytest.raises(APIClientError):
                client.request("GET", "/forbidden")

        assert breaker.is_open

    def test_single_404_records_a_failure_not_a_success(self, tmp_path: Path) -> None:
        breaker = CircuitBreaker(failure_threshold=5, service_name="test")
        client = _client(404, breaker, tmp_path)

        with pytest.raises(APIClientError):
            client.request("GET", "/missing")

        assert breaker.failure_count == 1
