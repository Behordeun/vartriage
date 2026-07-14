"""Unit tests for the VEP client with mocked HTTP responses."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.vep_client import VEPClient
from vartriage.models.variant import FunctionalConsequence

# --- Fixtures ---


@pytest.fixture
def rate_limiter() -> RateLimiter:
    return RateLimiter(tokens_per_second=1000.0, burst=100, service_name="vep")


@pytest.fixture
def circuit_breaker() -> CircuitBreaker:
    return CircuitBreaker(
        failure_threshold=5, recovery_timeout=60.0, service_name="vep"
    )


@pytest.fixture
def cache(tmp_path: Path) -> ResponseCache:
    return ResponseCache(db_path=tmp_path / "vep_test.db", default_ttl_days=1)


def _make_vep_response_entry(
    input_str: str,
    consequence: str = "missense_variant",
    gene_symbol: str = "BRCA1",
    cadd_phred: float | None = 24.5,
    gnomad_af: float | None = 0.002,
    canonical: int = 1,
) -> dict:
    """Build a single VEP response entry for testing."""
    entry: dict = {
        "input": input_str,
        "most_severe_consequence": consequence,
        "transcript_consequences": [
            {
                "gene_symbol": gene_symbol,
                "consequence_terms": [consequence],
                "transcript_id": "ENST00000000001",
                "canonical": canonical,
                "hgvsc": "ENST00000000001.1:c.123A>G",
                "hgvsp": "ENSP00000000001.1:p.Lys41Arg",
            }
        ],
    }
    if cadd_phred is not None:
        entry["transcript_consequences"][0]["cadd_phred"] = cadd_phred

    if gnomad_af is not None:
        entry["colocated_variants"] = [
            {"frequencies": {"A": {"gnomade_af": gnomad_af}}}
        ]

    return entry


def _build_vep_client(
    rate_limiter: RateLimiter,
    circuit_breaker: CircuitBreaker,
    cache: ResponseCache,
    responses: list[tuple[int, list | dict | str]],
) -> VEPClient:
    """Build a VEPClient with mocked HTTP transport."""
    call_count = [0]

    def handler(request: httpx.Request) -> httpx.Response:
        idx = min(call_count[0], len(responses) - 1)
        call_count[0] += 1
        status, body = responses[idx]
        if isinstance(body, (list, dict)):
            return httpx.Response(
                status,
                json=body,
                headers={"x-ensembl-release": "112"},
            )
        return httpx.Response(status, text=body)

    client = VEPClient(
        rate_limiter=rate_limiter,
        cache=cache,
        circuit_breaker=circuit_breaker,
        genome_build="grch38",
        batch_size=200,
        max_retries=2,
    )
    # Replace internal httpx client with mocked transport
    client._http._client = httpx.Client(
        base_url="https://rest.ensembl.org",
        transport=httpx.MockTransport(handler),
        timeout=httpx.Timeout(5.0),
    )
    return client


# --- Tests ---


class TestSuccessfulAnnotation:
    """Happy path: VEP returns valid annotations."""

    def test_single_variant_annotated(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        vep_response = [
            _make_vep_response_entry(
                "22 17818804 17818804 G/A +",
                consequence="missense_variant",
                gene_symbol="MICAL3",
                cadd_phred=24.5,
                gnomad_af=0.00215,
            )
        ]
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, vep_response)]
        )

        results = client.annotate_batch([("chr22", 17818804, "G", "A")])

        assert len(results) == 1
        ann = results[0]
        assert ann is not None
        assert ann.consequence == FunctionalConsequence.MISSENSE
        assert ann.gene_name == "MICAL3"
        assert ann.cadd_phred == 24.5
        assert ann.allele_frequency == pytest.approx(0.00215)

    def test_multiple_variants_in_single_batch(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        vep_response = [
            _make_vep_response_entry(
                "22 17818804 17818804 G/A +", gene_symbol="MICAL3"
            ),
            _make_vep_response_entry(
                "22 19029764 19029765 -/C +",
                consequence="frameshift_variant",
                gene_symbol="DGCR5",
            ),
        ]
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, vep_response)]
        )

        results = client.annotate_batch(
            [
                ("chr22", 17818804, "G", "A"),
                ("chr22", 19029764, "A", "AC"),
            ]
        )

        assert len(results) == 2
        assert results[0] is not None
        assert results[0].gene_name == "MICAL3"
        assert results[1] is not None
        assert results[1].consequence == FunctionalConsequence.FRAMESHIFT
        assert results[1].gene_name == "DGCR5"

    def test_captures_ensembl_release_version(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        vep_response = [_make_vep_response_entry("1 100 100 A/T +")]
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, vep_response)]
        )
        client.annotate_batch([("chr1", 100, "A", "T")])

        assert client.ensembl_release == "112"


class TestFrequencyExtraction:
    """gnomAD frequency parsing from colocated_variants."""

    def test_extracts_exome_frequency_by_default(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        entry = _make_vep_response_entry("1 100 100 A/T +", gnomad_af=0.05)
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, [entry])]
        )

        results = client.annotate_batch([("chr1", 100, "A", "T")])
        assert results[0] is not None
        assert results[0].allele_frequency == pytest.approx(0.05)

    def test_no_colocated_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        entry = _make_vep_response_entry("1 100 100 A/T +", gnomad_af=None)
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, [entry])]
        )

        results = client.annotate_batch([("chr1", 100, "A", "T")])
        assert results[0] is not None
        assert results[0].allele_frequency is None


class TestCADDExtraction:
    """CADD Phred score from VEP plugin data."""

    def test_extracts_cadd_phred(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        entry = _make_vep_response_entry("1 200 200 C/T +", cadd_phred=32.0)
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, [entry])]
        )

        results = client.annotate_batch([("1", 200, "C", "T")])
        assert results[0] is not None
        assert results[0].cadd_phred == 32.0

    def test_missing_cadd_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        entry = _make_vep_response_entry("1 200 200 C/T +", cadd_phred=None)
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, [entry])]
        )

        results = client.annotate_batch([("1", 200, "C", "T")])
        assert results[0] is not None
        assert results[0].cadd_phred is None


class TestCaching:
    """Response caching behavior."""

    def test_second_call_uses_cache(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        entry = _make_vep_response_entry("22 100 100 G/A +", gene_symbol="TP53")
        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            return httpx.Response(
                200, json=[entry], headers={"x-ensembl-release": "112"}
            )

        client = VEPClient(
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
        )
        client._http._client = httpx.Client(
            base_url="https://rest.ensembl.org",
            transport=httpx.MockTransport(handler),
            timeout=httpx.Timeout(5.0),
        )

        # First call hits API
        client.annotate_batch([("chr22", 100, "G", "A")])
        assert call_count[0] == 1

        # Second call should use cache (no new HTTP call)
        results = client.annotate_batch([("chr22", 100, "G", "A")])
        assert call_count[0] == 1
        assert results[0] is not None
        assert results[0].gene_name == "TP53"


class TestErrorHandling:
    """Graceful degradation on API failures."""

    def test_server_error_returns_none_for_failed_variants(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        client = _build_vep_client(
            rate_limiter,
            circuit_breaker,
            cache,
            [(500, "Internal Server Error")] * 3,
        )

        results = client.annotate_batch([("chr1", 100, "A", "T")])
        assert len(results) == 1
        assert results[0] is None

    def test_malformed_json_returns_none(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        # Return valid HTTP 200 but with non-list JSON
        client = _build_vep_client(
            rate_limiter,
            circuit_breaker,
            cache,
            [(200, {"error": "something unexpected"})],
        )

        results = client.annotate_batch([("chr1", 100, "A", "T")])
        assert results[0] is None

    def test_missing_input_field_in_response_skips_variant(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        # VEP response entry without 'input' field
        bad_entry = {"most_severe_consequence": "missense_variant"}
        client = _build_vep_client(
            rate_limiter,
            circuit_breaker,
            cache,
            [(200, [bad_entry])],
        )

        results = client.annotate_batch([("chr1", 100, "A", "T")])
        assert results[0] is None


class TestBatchFailureRecovery:
    """Progressive recovery: split batch, then individual queries."""

    def test_recovery_splits_batch_on_full_failure(
        self,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        # First call (full batch) fails, second and third (halves) succeed
        entry1 = _make_vep_response_entry("22 100 100 G/A +", gene_symbol="GENE1")
        entry2 = _make_vep_response_entry("22 200 200 C/T +", gene_symbol="GENE2")

        call_count = [0]

        def handler(request: httpx.Request) -> httpx.Response:
            call_count[0] += 1
            if call_count[0] <= 2:
                # First 2 calls fail (original batch retry)
                return httpx.Response(500, text="fail")
            # Recovery calls succeed with appropriate entries
            body = request.content.decode()
            if "100" in body:
                return httpx.Response(
                    200, json=[entry1], headers={"x-ensembl-release": "112"}
                )
            return httpx.Response(
                200, json=[entry2], headers={"x-ensembl-release": "112"}
            )

        client = VEPClient(
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=CircuitBreaker(failure_threshold=20, service_name="vep"),
            max_retries=2,
        )
        client._http._client = httpx.Client(
            base_url="https://rest.ensembl.org",
            transport=httpx.MockTransport(handler),
            timeout=httpx.Timeout(5.0),
        )

        results = client.annotate_batch(
            [
                ("chr22", 100, "G", "A"),
                ("chr22", 200, "C", "T"),
            ]
        )

        # At least one variant should have been recovered
        recovered = [r for r in results if r is not None]
        assert len(recovered) >= 1


class TestConsequenceMapping:
    """VEP response consequence correctly maps to our enum."""

    @pytest.mark.parametrize(
        "vep_term,expected",
        [
            ("missense_variant", FunctionalConsequence.MISSENSE),
            ("frameshift_variant", FunctionalConsequence.FRAMESHIFT),
            ("stop_gained", FunctionalConsequence.NONSENSE),
            ("splice_donor_variant", FunctionalConsequence.SPLICE_SITE),
            ("synonymous_variant", FunctionalConsequence.SYNONYMOUS),
            ("intergenic_variant", FunctionalConsequence.INTERGENIC),
            ("inframe_insertion", FunctionalConsequence.IN_FRAME_INSERTION),
            ("inframe_deletion", FunctionalConsequence.IN_FRAME_DELETION),
        ],
    )
    def test_consequence_mapping(
        self,
        vep_term: str,
        expected: FunctionalConsequence,
        rate_limiter: RateLimiter,
        circuit_breaker: CircuitBreaker,
        cache: ResponseCache,
    ) -> None:
        entry = _make_vep_response_entry("1 100 100 A/T +", consequence=vep_term)
        client = _build_vep_client(
            rate_limiter, circuit_breaker, cache, [(200, [entry])]
        )

        results = client.annotate_batch([("1", 100, "A", "T")])
        assert results[0] is not None
        assert results[0].consequence == expected
