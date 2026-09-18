"""Unit tests for GnomADClient not-found miss caching.

gnomAD signals a genuinely-absent variant with a GraphQL error whose message
is "Variant not found" (HTTP 200, data.variant == null). That absence is
cacheable: re-querying it on every run is what made large validation runs pay
the per-request rate limit repeatedly. Transient GraphQL errors (rate limit,
timeout, schema drift) must stay uncached so a later run can retry them.
"""

from __future__ import annotations

from pathlib import Path

import pytest

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.gnomad_client import GnomADClient


@pytest.fixture
def rate_limiter() -> RateLimiter:
    return RateLimiter(tokens_per_second=1000.0, burst=100, service_name="gnomad")


@pytest.fixture
def circuit_breaker() -> CircuitBreaker:
    return CircuitBreaker(
        failure_threshold=5, recovery_timeout=60.0, service_name="gnomad"
    )


@pytest.fixture
def cache(tmp_path: Path) -> ResponseCache:
    return ResponseCache(db_path=tmp_path / "gnomad_test.db", default_ttl_days=30)


def _build_client(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
    responses: list[dict],
) -> tuple[GnomADClient, list[int]]:
    """Build a GnomADClient whose HTTP transport returns canned GraphQL bodies."""
    call_count = [0]

    def handler(_request: httpx.Request) -> httpx.Response:
        idx = min(call_count[0], len(responses) - 1)
        call_count[0] += 1
        return httpx.Response(200, json=responses[idx])

    client = GnomADClient(
        rate_limiter=rate_limiter,
        cache=cache,
        circuit_breaker=circuit_breaker,
        dataset="gnomad_r4",
        prefer_source="combined",
        max_retries=2,
    )
    client._http._client = httpx.Client(
        base_url="https://gnomad.broadinstitute.org/api",
        transport=httpx.MockTransport(handler),
        timeout=httpx.Timeout(5.0),
    )
    return client, call_count


_NOT_FOUND_BODY = {
    "errors": [{"message": "Variant not found"}],
    "data": {"variant": None},
}
_TRANSIENT_BODY = {
    "errors": [{"message": "Request timed out"}],
    "data": {"variant": None},
}


def test_variant_not_found_is_cached_as_miss(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
) -> None:
    client, calls = _build_client(
        rate_limiter, circuit_breaker, cache, [_NOT_FOUND_BODY]
    )

    first = client.lookup_frequency("1", 99999999, "A", "T")
    assert first is None
    assert calls[0] == 1

    second = client.lookup_frequency("1", 99999999, "A", "T")
    assert second is None
    assert calls[0] == 1, (
        "a genuine not-found must be served from cache, not re-queried"
    )


def test_transient_error_is_not_cached(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
) -> None:
    client, calls = _build_client(
        rate_limiter, circuit_breaker, cache, [_TRANSIENT_BODY]
    )

    first = client.lookup_frequency("1", 12345, "C", "G")
    assert first is None
    assert calls[0] == 1

    second = client.lookup_frequency("1", 12345, "C", "G")
    assert second is None
    assert calls[0] == 2, "a transient error must be retried, never cached as a miss"


_MALFORMED_BODY = {
    "errors": ["a bare string, not an object"],
    "data": {"variant": None},
}
_MIXED_BODY = {
    "errors": [{"message": "Variant not found"}, {"message": "Rate limit exceeded"}],
    "data": {"variant": None},
}


def test_malformed_error_entry_is_not_cached(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
) -> None:
    # A non-dict error entry must not crash and must not be cached as a miss.
    client, calls = _build_client(
        rate_limiter, circuit_breaker, cache, [_MALFORMED_BODY]
    )
    assert client.lookup_frequency("1", 222, "C", "G") is None
    assert calls[0] == 1
    assert client.lookup_frequency("1", 222, "C", "G") is None
    assert calls[0] == 2, "a malformed error response must be retried, not cached"


def test_not_found_mixed_with_transient_is_not_cached(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
) -> None:
    # "Variant not found" alongside a transient error is not a clean absence.
    client, calls = _build_client(rate_limiter, circuit_breaker, cache, [_MIXED_BODY])
    assert client.lookup_frequency("1", 333, "A", "G") is None
    assert calls[0] == 1
    assert client.lookup_frequency("1", 333, "A", "G") is None
    assert calls[0] == 2, "a mixed not-found/transient response must be retried"
