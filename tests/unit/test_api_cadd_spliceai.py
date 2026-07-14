"""Unit tests for the CADD and SpliceAI API clients."""

from __future__ import annotations

from pathlib import Path

import pytest

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.cadd_client import CADDClient
from vartriage.api.spliceai_client import SpliceAIClient
from vartriage.models.variant import FunctionalConsequence

# --- Shared fixtures ---


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
    return ResponseCache(db_path=tmp_path / "score_test.db", default_ttl_days=7)


# --- CADD Client Tests ---


def _build_cadd_client(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
    handler: object,
) -> CADDClient:
    client = CADDClient(
        rate_limiter=rate_limiter,
        cache=cache,
        circuit_breaker=circuit_breaker,
        max_retries=2,
    )
    client._http._client = httpx.Client(
        base_url="https://cadd.gs.washington.edu/api/v1.0",
        transport=httpx.MockTransport(handler),
        timeout=httpx.Timeout(5.0),
    )
    return client


class TestCADDSuccessfulLookup:
    """CADD score retrieval happy path."""

    def test_extracts_phred_score_for_matching_allele(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        # CADD API returns list-of-lists: header + data rows
        cadd_response = [
            ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
            ["22", 17818804, "G", "A", 4.5, 24.3],
            ["22", 17818804, "G", "C", 3.2, 18.7],
            ["22", 17818804, "G", "T", 2.1, 14.2],
        ]

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json=cadd_response)

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr22", 17818804, "G", "A")

        assert score == pytest.approx(24.3)

    def test_returns_none_when_allele_not_in_response(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        cadd_response = [
            ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
            ["1", 100, "A", "G", 1.0, 10.0],
            ["1", 100, "A", "C", 0.5, 5.0],
        ]

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json=cadd_response)

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        # Query for A>T which isn't in the response
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None

    def test_batch_lookup_returns_positional_results(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            url = str(request.url)
            if "100" in url:
                return httpx.Response(
                    200,
                    json=[
                        ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
                        ["1", 100, "A", "T", 3.0, 20.0],
                    ],
                )
            return httpx.Response(
                200,
                json=[
                    ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
                    ["1", 200, "C", "G", 2.0, 15.0],
                ],
            )

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch(
            [
                ("chr1", 100, "A", "T"),
                ("chr1", 200, "C", "G"),
            ]
        )

        assert results[0] == pytest.approx(20.0)
        assert results[1] == pytest.approx(15.0)


class TestCADDCaching:
    """CADD response caching."""

    def test_second_lookup_uses_cache(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(
                200,
                json=[
                    ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
                    ["22", 500, "G", "A", 5.0, 30.0],
                ],
            )

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        client.lookup_score("chr22", 500, "G", "A")
        client.lookup_score("chr22", 500, "G", "A")

        assert call_count[0] == 1

    def test_cache_stores_none_for_miss(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            # No matching allele in response
            return httpx.Response(
                200,
                json=[
                    ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
                ],
            )

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        client.lookup_score("chr1", 100, "A", "T")
        client.lookup_score("chr1", 100, "A", "T")

        # None result cached, no second API call
        assert call_count[0] == 1


class TestCADDErrorHandling:
    """CADD client error resilience."""

    def test_server_error_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(500, text="Internal Server Error")

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None

    def test_malformed_response_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, text="not json")

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None

    def test_empty_response_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json=[])

        client = _build_cadd_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None


# --- SpliceAI Client Tests ---


def _build_spliceai_client(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
    handler: object,
) -> SpliceAIClient:
    client = SpliceAIClient(
        rate_limiter=rate_limiter,
        cache=cache,
        circuit_breaker=circuit_breaker,
        max_retries=2,
    )
    client._http._client = httpx.Client(
        base_url="https://spliceailookup.broadinstitute.org",
        transport=httpx.MockTransport(handler),
        timeout=httpx.Timeout(5.0),
    )
    return client


class TestSpliceAISuccessfulLookup:
    """SpliceAI score retrieval happy path."""

    def test_extracts_max_delta_score(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        spliceai_response = {
            "scores": [
                {
                    "DS_AG": 0.01,
                    "DS_AL": 0.85,
                    "DS_DG": 0.02,
                    "DS_DL": 0.03,
                }
            ]
        }

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json=spliceai_response)

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr22", 100, "G", "A")

        assert score == pytest.approx(0.85)

    def test_max_across_multiple_score_entries(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        spliceai_response = {
            "scores": [
                {"DS_AG": 0.1, "DS_AL": 0.2, "DS_DG": 0.3, "DS_DL": 0.4},
                {"DS_AG": 0.9, "DS_AL": 0.1, "DS_DG": 0.05, "DS_DL": 0.02},
            ]
        }

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json=spliceai_response)

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 50, "A", "T")

        assert score == pytest.approx(0.9)

    def test_all_zeros_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        spliceai_response = {
            "scores": [
                {"DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.0, "DS_DL": 0.0},
            ]
        }

        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json=spliceai_response)

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 50, "A", "T")

        assert score is None


class TestSpliceAISmartFiltering:
    """Consequence-based query filtering to save rate limit tokens."""

    def test_skips_intergenic_variants(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(200, json={"scores": []})

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score(
            "chr1",
            100,
            "A",
            "T",
            consequence=FunctionalConsequence.INTERGENIC,
        )

        assert score is None
        assert call_count[0] == 0  # No API call made

    def test_skips_frameshift_variants(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(200, json={"scores": []})

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score(
            "chr1",
            100,
            "A",
            "T",
            consequence=FunctionalConsequence.FRAMESHIFT,
        )

        assert score is None
        assert call_count[0] == 0

    def test_queries_missense_variants(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(
                200,
                json={
                    "scores": [{"DS_AG": 0.5, "DS_AL": 0.0, "DS_DG": 0.0, "DS_DL": 0.0}]
                },
            )

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score(
            "chr1",
            100,
            "A",
            "T",
            consequence=FunctionalConsequence.MISSENSE,
        )

        assert score == pytest.approx(0.5)
        assert call_count[0] == 1  # API call was made

    def test_queries_splice_site_variants(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(
                200,
                json={
                    "scores": [{"DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.7, "DS_DL": 0.0}]
                },
            )

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score(
            "chr1",
            100,
            "A",
            "T",
            consequence=FunctionalConsequence.SPLICE_SITE,
        )

        assert score == pytest.approx(0.7)
        assert call_count[0] == 1

    def test_no_consequence_filter_queries_unconditionally(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(200, json={"scores": []})

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        client.lookup_score("chr1", 100, "A", "T", consequence=None)

        assert call_count[0] == 1

    def test_batch_lookup_applies_filtering(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(
                200,
                json={
                    "scores": [{"DS_AG": 0.6, "DS_AL": 0.0, "DS_DG": 0.0, "DS_DL": 0.0}]
                },
            )

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        results = client.lookup_batch(
            [("chr1", 100, "A", "T"), ("chr2", 200, "G", "C"), ("chr3", 300, "T", "A")],
            consequences=[
                FunctionalConsequence.MISSENSE,  # queried
                FunctionalConsequence.INTERGENIC,  # skipped
                FunctionalConsequence.SPLICE_SITE,  # queried
            ],
        )

        assert results[0] == pytest.approx(0.6)
        assert results[1] is None
        assert results[2] == pytest.approx(0.6)
        assert call_count[0] == 2  # Only 2 API calls, not 3


class TestSpliceAICaching:
    """SpliceAI response caching."""

    def test_cached_score_avoids_api_call(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(
                200,
                json={
                    "scores": [{"DS_AG": 0.4, "DS_AL": 0.0, "DS_DG": 0.0, "DS_DL": 0.0}]
                },
            )

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        client.lookup_score("chr1", 100, "A", "T")
        client.lookup_score("chr1", 100, "A", "T")

        assert call_count[0] == 1


class TestSpliceAIErrorHandling:
    """SpliceAI error resilience."""

    def test_server_error_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(500, text="Error")

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None

    def test_empty_scores_list_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json={"scores": []})

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None

    def test_malformed_response_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, text="not json")

        client = _build_spliceai_client(rate_limiter, circuit_breaker, cache, handler)
        score = client.lookup_score("chr1", 100, "A", "T")

        assert score is None
