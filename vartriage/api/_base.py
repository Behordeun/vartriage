"""Base HTTP client with retry, rate limiting, circuit breaking, and caching.

All service-specific API clients inherit from BaseAPIClient. This class
handles cross-cutting concerns so individual clients focus only on
request construction and response parsing.
"""

from __future__ import annotations

import logging
import time
from typing import Any

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._rate_limiter import DailyLimitExhausted, RateLimiter

logger = logging.getLogger(__name__)

# Transient HTTP status codes worth retrying
_RETRYABLE_STATUS_CODES = frozenset({429, 500, 502, 503, 504})


def _check_httpx_available() -> None:
    """Raise ImportError with install instructions if httpx is missing."""
    try:
        import httpx  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "API mode requires the 'httpx' package. "
            "Install with: pip install vartriage[api]"
        ) from exc


class APIClientError(Exception):
    """Base exception for API client failures."""

    def __init__(
        self, service: str, message: str, status_code: int | None = None
    ) -> None:
        self.service = service
        self.status_code = status_code
        super().__init__(f"[{service}] {message}")


class BaseAPIClient:
    """HTTP client with retry, rate limiting, circuit breaker, and structured logging.

    Subclasses override request construction and response parsing.
    This base handles all resilience and observability plumbing.

    Parameters
    ----------
    base_url
        Root URL for the API (no trailing slash).
    rate_limiter
        Token bucket for request throttling.
    cache
        Response cache for deduplication.
    circuit_breaker
        Circuit breaker for failure isolation.
    service_name
        Identifier for logging and error messages.
    timeout
        Tuple of (connect_timeout, read_timeout) in seconds.
    max_retries
        Maximum retry attempts for transient failures.
    user_agent
        User-Agent header value.
    proxy_url
        Optional HTTP/HTTPS proxy URL.
    """

    def __init__(
        self,
        base_url: str,
        rate_limiter: RateLimiter,
        cache: ResponseCache,
        circuit_breaker: CircuitBreaker,
        service_name: str = "api",
        timeout: tuple[float, float] = (10.0, 30.0),
        max_retries: int = 3,
        user_agent: str = "vartriage/0.7.0 (https://github.com/Behordeun/vartriage)",
        proxy_url: str | None = None,
    ) -> None:
        _check_httpx_available()

        self._base_url = base_url.rstrip("/")
        self._rate_limiter = rate_limiter
        self._cache = cache
        self._circuit_breaker = circuit_breaker
        self._service_name = service_name
        self._timeout = timeout
        self._max_retries = max_retries
        self._user_agent = user_agent
        self._proxy_url = proxy_url
        self._client = self._build_client()

    def _build_client(self) -> Any:
        """Construct the httpx.Client with connection pooling and proxy."""
        import httpx

        transport_kwargs: dict[str, Any] = {"retries": 0}  # we handle retries ourselves
        proxy_param: Any = self._proxy_url if self._proxy_url else None

        return httpx.Client(
            base_url=self._base_url,
            timeout=httpx.Timeout(
                connect=self._timeout[0],
                read=self._timeout[1],
                write=self._timeout[1],
                pool=self._timeout[0],
            ),
            headers={"User-Agent": self._user_agent},
            proxy=proxy_param,
            transport=httpx.HTTPTransport(**transport_kwargs),
            follow_redirects=True,
        )

    def request(
        self,
        method: str,
        path: str,
        json_body: Any = None,
        params: dict[str, str] | None = None,
        headers: dict[str, str] | None = None,
    ) -> Any:
        """Execute an HTTP request with full resilience stack.

        Flow: rate_limiter.acquire() -> circuit_breaker.allow_request()
        -> httpx request -> retry on transient failure -> record outcome.

        Parameters
        ----------
        method
            HTTP method (GET, POST, etc.).
        path
            URL path appended to base_url.
        json_body
            JSON-serializable request body (for POST/PUT).
        params
            URL query parameters.
        headers
            Additional request headers (merged with defaults).

        Returns
        -------
        httpx.Response
            The successful HTTP response.

        Raises
        ------
        APIClientError
            On non-retryable failure after all attempts exhausted.
        CircuitBreakerOpen
            If the service circuit is open.
        DailyLimitExhausted
            If the daily request cap is reached.
        """
        import httpx

        # Circuit breaker check (raises CircuitBreakerOpen if tripped)
        self._circuit_breaker.allow_request()

        last_error: Exception | None = None
        last_status: int | None = None

        for attempt in range(1, self._max_retries + 1):
            # Rate limiting (blocks until token available)
            try:
                self._rate_limiter.acquire()
            except DailyLimitExhausted:
                raise

            start_time = time.monotonic()
            try:
                response = self._client.request(
                    method=method,
                    url=path,
                    json=json_body,
                    params=params,
                    headers=headers,
                )
                elapsed = time.monotonic() - start_time
                last_status = response.status_code

                logger.info(
                    "API %s %s %s status=%d latency=%.2fs attempt=%d/%d",
                    self._service_name,
                    method,
                    path,
                    response.status_code,
                    elapsed,
                    attempt,
                    self._max_retries,
                )

                if response.status_code < 400:
                    self._circuit_breaker.record_success()
                    return response

                if response.status_code not in _RETRYABLE_STATUS_CODES:
                    # Non-retryable client error (4xx except 429)
                    self._circuit_breaker.record_success()
                    raise APIClientError(
                        self._service_name,
                        f"HTTP {response.status_code}: {response.text[:200]}",
                        status_code=response.status_code,
                    )

                # Retryable status code
                last_error = APIClientError(
                    self._service_name,
                    f"HTTP {response.status_code} (attempt {attempt}/{self._max_retries})",
                    status_code=response.status_code,
                )

                # Respect Retry-After header on 429
                if response.status_code == 429:
                    retry_after = self._parse_retry_after(response)
                    if retry_after and retry_after < 120:
                        logger.info(
                            "Rate limited by %s, waiting %.1fs (Retry-After)",
                            self._service_name,
                            retry_after,
                        )
                        time.sleep(retry_after)
                        continue

                # Exponential backoff: 1s, 2s, 4s
                backoff = min(2 ** (attempt - 1), 8)
                time.sleep(backoff)

            except httpx.TimeoutException as exc:
                elapsed = time.monotonic() - start_time
                logger.warning(
                    "API %s %s %s timeout after %.2fs attempt=%d/%d",
                    self._service_name,
                    method,
                    path,
                    elapsed,
                    attempt,
                    self._max_retries,
                )
                last_error = exc
                backoff = min(2 ** (attempt - 1), 8)
                time.sleep(backoff)

            except httpx.ConnectError as exc:
                logger.warning(
                    "API %s %s %s connection error attempt=%d/%d: %s",
                    self._service_name,
                    method,
                    path,
                    attempt,
                    self._max_retries,
                    str(exc)[:100],
                )
                last_error = exc
                backoff = min(2 ** (attempt - 1), 8)
                time.sleep(backoff)

        # All retries exhausted
        self._circuit_breaker.record_failure()
        raise APIClientError(
            self._service_name,
            f"All {self._max_retries} attempts failed. Last status: {last_status}. "
            f"Error: {last_error}",
            status_code=last_status,
        )

    def close(self) -> None:
        """Close the underlying HTTP client and release connections."""
        self._client.close()

    def _parse_retry_after(self, response: Any) -> float | None:
        """Extract Retry-After header value in seconds."""
        retry_after = response.headers.get("Retry-After")
        if retry_after is None:
            return None
        try:
            return float(retry_after)
        except ValueError:
            return None
