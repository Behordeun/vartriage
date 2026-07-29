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


class _AttemptKind:
    """Outcome categories for a single HTTP attempt."""
    SUCCESS = "success"
    RETRYABLE = "retryable"
    NETWORK_ERROR = "network_error"


class _AttemptResult:
    """Structured result from _execute_attempt.

    Encodes what happened so the retry loop can make clear decisions
    without interpreting overloaded return types.
    """

    __slots__ = ("kind", "response", "error", "status_code", "already_delayed")

    def __init__(
        self,
        kind: str,
        response: Any = None,
        error: Exception | None = None,
        status_code: int | None = None,
        already_delayed: bool = False,
    ) -> None:
        self.kind = kind
        self.response = response
        self.error = error
        self.status_code = status_code
        self.already_delayed = already_delayed


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

    def _handle_retry_after(self, response: Any) -> bool:
        """Sleep for Retry-After duration on 429 if within bounds. Returns True if slept."""
        import httpx  # noqa: F401 — imported for type context only
        retry_after = self._parse_retry_after(response)
        if retry_after and retry_after < 120:
            logger.info(
                "Rate limited by %s, waiting %.1fs (Retry-After)",
                self._service_name,
                retry_after,
            )
            time.sleep(retry_after)
            return True
        return False

    def _execute_attempt(
        self,
        method: str,
        path: str,
        json_body: Any,
        params: dict[str, str] | None,
        headers: dict[str, str] | None,
        attempt: int,
    ) -> _AttemptResult:
        """Run one HTTP attempt and return a structured result."""
        import httpx

        start_time = time.monotonic()
        try:
            response = self._client.request(
                method=method, url=path, json=json_body,
                params=params, headers=headers,
            )
            elapsed = time.monotonic() - start_time
            logger.info(
                "API %s %s %s status=%d latency=%.2fs attempt=%d/%d",
                self._service_name, method, path,
                response.status_code, elapsed, attempt, self._max_retries,
            )
            if response.status_code < 400:
                self._circuit_breaker.record_success()
                return _AttemptResult(
                    kind=_AttemptKind.SUCCESS,
                    response=response,
                    status_code=response.status_code,
                )
            if response.status_code not in _RETRYABLE_STATUS_CODES:
                self._circuit_breaker.record_success()
                raise APIClientError(
                    self._service_name,
                    f"HTTP {response.status_code}: {response.text[:200]}",
                    status_code=response.status_code,
                )
            # Retryable HTTP error
            already_delayed = False
            if response.status_code == 429 and self._handle_retry_after(response):
                already_delayed = True
            return _AttemptResult(
                kind=_AttemptKind.RETRYABLE,
                status_code=response.status_code,
                already_delayed=already_delayed,
            )
        except (httpx.TimeoutException, httpx.ConnectError) as exc:
            elapsed = time.monotonic() - start_time
            logger.warning(
                "API %s %s %s network error after %.2fs attempt=%d/%d: %s",
                self._service_name, method, path,
                elapsed, attempt, self._max_retries, str(exc)[:100],
            )
            return _AttemptResult(
                kind=_AttemptKind.NETWORK_ERROR,
                error=exc,
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

        # Circuit breaker check (raises CircuitBreakerOpen if tripped)
        self._circuit_breaker.allow_request()

        last_error: Exception | None = None
        last_status: int | None = None

        for attempt in range(1, self._max_retries + 1):
            self._rate_limiter.acquire()

            result = self._execute_attempt(
                method, path, json_body, params, headers, attempt
            )
            if result.status_code is not None:
                last_status = result.status_code
            if result.kind == _AttemptKind.SUCCESS:
                return result.response
            if result.kind == _AttemptKind.NETWORK_ERROR:
                last_error = result.error
            else:
                last_error = APIClientError(
                    self._service_name,
                    f"HTTP {result.status_code} (attempt {attempt}/{self._max_retries})",
                    status_code=result.status_code,
                )
            # Skip extra backoff if Retry-After already applied
            if not result.already_delayed:
                backoff = min(2 ** (attempt - 1), 8)
                time.sleep(backoff)

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
