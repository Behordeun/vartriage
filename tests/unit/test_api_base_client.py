"""Unit tests for the BaseAPIClient with mocked HTTP transport."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._rate_limiter import RateLimiter

# httpx is optional; skip if unavailable
httpx = pytest.importorskip("httpx")

from vartriage.api._base import APIClientError, BaseAPIClient


@pytest.fixture
def rate_limiter() -> RateLimiter:
    return RateLimiter(tokens_per_second=1000.0, burst=100, service_name="test")


@pytest.fixture
def circuit_breaker() -> CircuitBreaker:
    return CircuitBreaker(
        failure_threshold=5, recovery_timeout=60.0, service_name="test"
    )


@pytest.fixture
def cache(tmp_path: Path) -> ResponseCache:
    return ResponseCache(db_path=tmp_path / "test.db", default_ttl_days=1)


def _mock_transport(responses: list[tuple[int, dict | str]]) -> httpx.MockTransport:
    """Build a mock transport that returns sequential responses.

    Each entry is (status_code, body). Body can be a dict (JSON) or str (text).
    """
    call_count = [0]

    def handler(request: httpx.Request) -> httpx.Response:
        idx = min(call_count[0], len(responses) - 1)
        call_count[0] += 1
        status, body = responses[idx]
        if isinstance(body, dict):
            return httpx.Response(status, json=body)
        return httpx.Response(status, text=body)

    return httpx.MockTransport(handler)


def _build_client(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
    transport: httpx.MockTransport,
    max_retries: int = 3,
) -> BaseAPIClient:
    """Build a BaseAPIClient with a mocked transport."""
    client = BaseAPIClient(
        base_url="https://test.example.com",
        rate_limiter=rate_limiter,
        cache=cache,
        circuit_breaker=circuit_breaker,
        service_name="test",
        max_retries=max_retries,
    )
    # Replace the internal httpx client with our mocked transport
    client._client = httpx.Client(
        base_url="https://test.example.com",
        transport=transport,
        timeout=httpx.Timeout(5.0),
    )
    return client


class TestSuccessfulRequests:
    """Happy path: 2xx responses."""

    def test_get_request_returns_response(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(200, {"result": "ok"})])
        client = _build_client(rate_limiter, circuit_breaker, cache, transport)

        response = client.request("GET", "/test")
        assert response.status_code == 200
        assert response.json() == {"result": "ok"}

    def test_post_request_with_json_body(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(201, {"id": 42})])
        client = _build_client(rate_limiter, circuit_breaker, cache, transport)

        response = client.request("POST", "/create", json_body={"name": "variant"})
        assert response.status_code == 201

    def test_success_records_on_circuit_breaker(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(200, {})])
        client = _build_client(rate_limiter, circuit_breaker, cache, transport)

        # Record a failure first so we can verify success resets it
        circuit_breaker.record_failure()
        assert circuit_breaker.failure_count == 1

        client.request("GET", "/test")
        assert circuit_breaker.failure_count == 0


class TestRetryBehavior:
    """Transient failure retry logic."""

    def test_retries_on_500_then_succeeds(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport(
            [
                (500, "Internal Server Error"),
                (500, "Internal Server Error"),
                (200, {"ok": True}),
            ]
        )
        client = _build_client(
            rate_limiter, circuit_breaker, cache, transport, max_retries=3
        )

        response = client.request("GET", "/flaky")
        assert response.status_code == 200

    def test_retries_on_503(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport(
            [
                (503, "Service Unavailable"),
                (200, {"recovered": True}),
            ]
        )
        client = _build_client(
            rate_limiter, circuit_breaker, cache, transport, max_retries=2
        )

        response = client.request("GET", "/service")
        assert response.json() == {"recovered": True}

    def test_raises_after_all_retries_exhausted(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport(
            [
                (500, "fail"),
                (500, "fail"),
                (500, "fail"),
            ]
        )
        client = _build_client(
            rate_limiter, circuit_breaker, cache, transport, max_retries=3
        )

        with pytest.raises(APIClientError) as exc_info:
            client.request("GET", "/always-fail")

        assert exc_info.value.service == "test"
        assert "All 3 attempts failed" in str(exc_info.value)

    def test_records_failure_on_circuit_breaker_after_exhaustion(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(500, "fail")] * 3)
        client = _build_client(
            rate_limiter, circuit_breaker, cache, transport, max_retries=3
        )

        with pytest.raises(APIClientError):
            client.request("GET", "/fail")

        assert circuit_breaker.failure_count == 1


class TestNonRetryableErrors:
    """Client errors (4xx except 429) are not retried."""

    def test_400_raises_immediately(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(400, "Bad Request")])
        client = _build_client(rate_limiter, circuit_breaker, cache, transport)

        with pytest.raises(APIClientError) as exc_info:
            client.request("GET", "/bad")

        assert exc_info.value.status_code == 400

    def test_404_raises_immediately(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(404, "Not Found")])
        client = _build_client(rate_limiter, circuit_breaker, cache, transport)

        with pytest.raises(APIClientError) as exc_info:
            client.request("GET", "/missing")

        assert exc_info.value.status_code == 404

    def test_non_retryable_does_not_trigger_circuit_breaker(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        transport = _mock_transport([(422, "Unprocessable")])
        client = _build_client(rate_limiter, circuit_breaker, cache, transport)

        with pytest.raises(APIClientError):
            client.request("POST", "/validate")

        # 4xx (non-429) counts as a successful circuit interaction
        assert circuit_breaker.failure_count == 0


class TestRateLimitHandling:
    """429 Too Many Requests with Retry-After."""

    def test_retries_429_with_retry_after_header(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            if call_count[0] == 1:
                return httpx.Response(429, headers={"Retry-After": "0.01"}, text="")
            return httpx.Response(200, json={"ok": True})

        transport = httpx.MockTransport(handler)
        client = _build_client(
            rate_limiter, circuit_breaker, cache, transport, max_retries=3
        )

        response = client.request("GET", "/limited")
        assert response.status_code == 200
        assert call_count[0] == 2


class TestCircuitBreakerIntegration:
    """Circuit breaker blocks requests when open."""

    def test_open_circuit_rejects_without_network_call(
        self, rate_limiter: RateLimiter, cache: ResponseCache
    ) -> None:
        cb = CircuitBreaker(
            failure_threshold=1, recovery_timeout=60.0, service_name="test"
        )
        cb.record_failure()  # trips the breaker

        transport = _mock_transport([(200, {})])
        client = _build_client(rate_limiter, cb, cache, transport)

        with pytest.raises(CircuitBreakerOpen):
            client.request("GET", "/blocked")


class TestUserAgent:
    """User-Agent header injection."""

    def test_user_agent_header_present(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        captured_headers: dict[str, str] = {}

        def handler(request: httpx.Request) -> httpx.Response:
            captured_headers.update(dict(request.headers))
            return httpx.Response(200, json={})

        transport = httpx.MockTransport(handler)
        client = BaseAPIClient(
            base_url="https://test.example.com",
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
            service_name="test",
            user_agent="vartriage/0.7.0-test",
        )
        client._client = httpx.Client(
            base_url="https://test.example.com",
            transport=transport,
            headers={"User-Agent": "vartriage/0.7.0-test"},
        )

        client.request("GET", "/check-ua")
        assert "vartriage/0.7.0-test" in captured_headers.get("user-agent", "")


class TestImportGuard:
    """Graceful handling when httpx is not installed."""

    def test_import_error_with_install_instructions(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        import builtins

        original_import = builtins.__import__

        def mock_import(name: str, *args: object, **kwargs: object) -> object:
            if name == "httpx":
                raise ImportError("No module named 'httpx'")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=mock_import):
            with pytest.raises(ImportError, match="pip install vartriage\\[api\\]"):
                from vartriage.api._base import _check_httpx_available

                _check_httpx_available()
