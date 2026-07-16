"""Unit tests for the ClinVar E-utilities client."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.clinvar_client import ClinVarClient, _map_significance
from vartriage.models.variant import ClinVarAssertion


@pytest.fixture
def rate_limiter() -> RateLimiter:
    return RateLimiter(tokens_per_second=1000.0, burst=100, service_name="clinvar")


@pytest.fixture
def circuit_breaker() -> CircuitBreaker:
    return CircuitBreaker(
        failure_threshold=5, recovery_timeout=60.0, service_name="clinvar"
    )


@pytest.fixture
def cache(tmp_path: Path) -> ResponseCache:
    return ResponseCache(db_path=tmp_path / "clinvar_test.db", default_ttl_days=30)


def _build_client(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
    handler: object,
    api_key: str | None = None,
) -> ClinVarClient:
    """Build a ClinVarClient with mocked HTTP transport."""
    client = ClinVarClient(
        rate_limiter=rate_limiter,
        cache=cache,
        circuit_breaker=circuit_breaker,
        ncbi_api_key=api_key,
        max_retries=2,
    )
    client._http._client = httpx.Client(
        base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
        transport=httpx.MockTransport(handler),
        timeout=httpx.Timeout(5.0),
    )
    return client


class TestSignificanceMapping:
    """Clinical significance string to enum mapping."""

    def test_pathogenic(self) -> None:
        assert _map_significance("Pathogenic") == ClinVarAssertion.PATHOGENIC

    def test_likely_pathogenic(self) -> None:
        assert (
            _map_significance("Likely pathogenic") == ClinVarAssertion.LIKELY_PATHOGENIC
        )

    def test_pathogenic_likely_pathogenic_compound(self) -> None:
        assert (
            _map_significance("Pathogenic/Likely pathogenic")
            == ClinVarAssertion.LIKELY_PATHOGENIC
        )

    def test_uncertain_significance(self) -> None:
        assert _map_significance("Uncertain significance") == ClinVarAssertion.VUS

    def test_likely_benign(self) -> None:
        assert _map_significance("Likely benign") == ClinVarAssertion.LIKELY_BENIGN

    def test_benign(self) -> None:
        assert _map_significance("Benign") == ClinVarAssertion.BENIGN

    def test_benign_likely_benign_compound(self) -> None:
        assert _map_significance("Benign/Likely benign") == ClinVarAssertion.BENIGN

    def test_conflicting_interpretations(self) -> None:
        result = _map_significance("Conflicting interpretations of pathogenicity")
        assert result == ClinVarAssertion.VUS

    def test_empty_string_returns_none(self) -> None:
        assert _map_significance("") is None

    def test_unknown_significance_returns_none(self) -> None:
        assert _map_significance("not provided") is None

    def test_case_insensitive(self) -> None:
        assert _map_significance("PATHOGENIC") == ClinVarAssertion.PATHOGENIC
        assert _map_significance("benign") == ClinVarAssertion.BENIGN


class TestSuccessfulLookup:
    """Happy path: esearch finds ID, esummary returns significance."""

    def test_pathogenic_variant_lookup(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            url = str(request.url)
            if "esearch" in url:
                return httpx.Response(
                    200, json={"esearchresult": {"idlist": ["12345"]}}
                )
            if "esummary" in url:
                return httpx.Response(
                    200,
                    json={
                        "result": {
                            "12345": {
                                "clinical_significance": {
                                    "description": "Pathogenic",
                                    "review_status": "criteria provided, multiple submitters, no conflicts",
                                    "last_evaluated": "2025-01-15",
                                }
                            }
                        }
                    },
                )
            return httpx.Response(404)

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr17", 43094452, "G", "A")])

        assert len(results) == 1
        assert results[0] == ClinVarAssertion.PATHOGENIC

    def test_benign_variant_lookup(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            url = str(request.url)
            if "esearch" in url:
                return httpx.Response(
                    200, json={"esearchresult": {"idlist": ["99999"]}}
                )
            if "esummary" in url:
                return httpx.Response(
                    200,
                    json={
                        "result": {
                            "99999": {
                                "clinical_significance": {
                                    "description": "Benign",
                                    "review_status": "criteria provided, single submitter",
                                }
                            }
                        }
                    },
                )
            return httpx.Response(404)

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr22", 17818804, "G", "A")])

        assert results[0] == ClinVarAssertion.BENIGN


class TestNoResults:
    """Variants not found in ClinVar."""

    def test_esearch_returns_empty_list(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr1", 99999999, "A", "T")])

        assert results[0] is None

    def test_cache_miss_is_persisted(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        """A miss should be cached to avoid repeated lookups."""
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        client.lookup_batch([("chr1", 100, "A", "T")])
        client.lookup_batch([("chr1", 100, "A", "T")])

        # Second call should hit cache, not the API
        assert call_count[0] == 1


class TestConflictingInterpretations:
    """Multiple submitters with different assertions."""

    def test_picks_highest_review_status(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            url = str(request.url)
            if "esearch" in url:
                return httpx.Response(
                    200, json={"esearchresult": {"idlist": ["111", "222"]}}
                )
            if "esummary" in url:
                return httpx.Response(
                    200,
                    json={
                        "result": {
                            "111": {
                                "clinical_significance": {
                                    "description": "Uncertain significance",
                                    "review_status": "criteria provided, single submitter",
                                }
                            },
                            "222": {
                                "clinical_significance": {
                                    "description": "Pathogenic",
                                    "review_status": "reviewed by expert panel",
                                }
                            },
                        }
                    },
                )
            return httpx.Response(404)

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr13", 32914437, "T", "A")])

        # Expert panel (rank 3) beats single submitter (rank 1)
        assert results[0] == ClinVarAssertion.PATHOGENIC


class TestAPIKeyInjection:
    """NCBI API key handling."""

    def test_api_key_included_in_params(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        captured_urls: list[str] = []

        def handler(request: httpx.Request) -> httpx.Response:
            captured_urls.append(str(request.url))
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = _build_client(
            rate_limiter,
            circuit_breaker,
            cache,
            handler,
            api_key="test_key_12345",
        )
        client.lookup_batch([("chr1", 100, "A", "T")])

        assert any("api_key=test_key_12345" in url for url in captured_urls)

    def test_tool_and_email_included(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        captured_urls: list[str] = []

        def handler(request: httpx.Request) -> httpx.Response:
            captured_urls.append(str(request.url))
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        client.lookup_batch([("chr1", 100, "A", "T")])

        assert any("tool=vartriage" in url for url in captured_urls)
        assert any("email=" in url for url in captured_urls)


class TestErrorHandling:
    """Graceful degradation on API failures."""

    def test_esearch_failure_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(500, text="Server Error")

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr1", 100, "A", "T")])

        assert results[0] is None

    def test_esummary_failure_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            url = str(request.url)
            if "esearch" in url:
                return httpx.Response(
                    200, json={"esearchresult": {"idlist": ["12345"]}}
                )
            # esummary fails
            return httpx.Response(500, text="Internal Error")

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr1", 100, "A", "T")])

        assert results[0] is None

    def test_malformed_json_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, text="not json at all")

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch([("chr1", 100, "A", "T")])

        assert results[0] is None


class TestCaching:
    """Response caching behavior."""

    def test_cached_result_skips_api_call(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            url = str(request.url)
            if "esearch" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": ["555"]}})
            return httpx.Response(
                200,
                json={
                    "result": {
                        "555": {
                            "clinical_significance": {
                                "description": "Likely pathogenic",
                                "review_status": "criteria provided, single submitter",
                            }
                        }
                    }
                },
            )

        client = _build_client(rate_limiter, circuit_breaker, cache, handler)

        # First call hits API
        result1 = client.lookup_batch([("chr7", 5000, "C", "T")])
        assert result1[0] == ClinVarAssertion.LIKELY_PATHOGENIC

        first_call_count = call_count[0]

        # Second call should use cache
        result2 = client.lookup_batch([("chr7", 5000, "C", "T")])
        assert result2[0] == ClinVarAssertion.LIKELY_PATHOGENIC
        assert call_count[0] == first_call_count  # no new API calls
