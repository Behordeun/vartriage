"""CADD REST API client for Phred-scaled deleteriousness scores.

Queries the CADD web API at cadd.gs.washington.edu for pre-computed
scores. Used as a fallback when the VEP CADD plugin doesn't return
a score for a given variant (not all positions have pre-computed data).

The API returns all possible substitutions at a position, so we filter
the response by matching ref/alt allele.
"""

from __future__ import annotations

import logging
from typing import Optional

from vartriage.api._base import APIClientError, BaseAPIClient
from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._notation import _strip_chr_prefix
from vartriage.api._rate_limiter import DailyLimitExhausted, RateLimiter

logger = logging.getLogger(__name__)

_CADD_BASE_URL = "https://cadd.gs.washington.edu/api/v1.0"


class CADDClient:
    """CADD REST API client for Phred score lookups.

    Parameters
    ----------
    rate_limiter
        Token bucket for CADD rate limiting (2 req/sec default).
    cache
        Response cache for deduplication.
    circuit_breaker
        Circuit breaker for CADD endpoint.
    genome_build
        Target assembly: "grch38" or "grch37".
    max_retries
        Retry attempts for transient failures.
    timeout
        (connect, read) timeouts in seconds.
    user_agent
        User-Agent header value.
    proxy_url
        HTTP/HTTPS proxy URL.
    """

    def __init__(
        self,
        rate_limiter: RateLimiter,
        cache: ResponseCache,
        circuit_breaker: CircuitBreaker,
        genome_build: str = "grch38",
        max_retries: int = 3,
        timeout: tuple[float, float] = (10.0, 30.0),
        user_agent: str = "vartriage/0.7.0 (https://github.com/Behordeun/vartriage)",
        proxy_url: Optional[str] = None,
    ) -> None:
        self._genome_build = genome_build
        self._cache = cache

        # CADD API path includes the build
        build_path = "GRCh38" if genome_build == "grch38" else "GRCh37"
        self._build_path = build_path

        self._http = BaseAPIClient(
            base_url=_CADD_BASE_URL,
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
            service_name="cadd",
            timeout=timeout,
            max_retries=max_retries,
            user_agent=user_agent,
            proxy_url=proxy_url,
        )

    def lookup_score(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[float]:
        """Look up CADD Phred score for a single variant.

        Parameters
        ----------
        chrom
            Chromosome (e.g., "chr22", "22").
        pos
            1-based genomic position.
        ref
            Reference allele.
        alt
            Alternate allele.

        Returns
        -------
        Optional[float]
            CADD Phred score, or None if not found or API error.
        """
        cache_key = ResponseCache.build_key(
            "cadd", self._genome_build, chrom, pos, ref, alt
        )
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached.get("phred")

        chrom_clean = _strip_chr_prefix(chrom)

        # CADD API: /v1.0/GRCh38/chrom:pos
        path = f"/{self._build_path}/{chrom_clean}:{pos}"

        try:
            response = self._http.request("GET", path)
        except (APIClientError, CircuitBreakerOpen, DailyLimitExhausted) as exc:
            logger.warning("CADD lookup failed for %s:%d: %s", chrom, pos, exc)
            return None

        phred = self._parse_response(response, ref, alt)

        # Cache regardless of hit/miss to avoid re-querying
        self._cache.put(
            key=cache_key,
            value={"phred": phred},
            source="cadd",
            genome_build=self._genome_build,
        )

        return phred

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[float]]:
        """Look up CADD scores for multiple variants.

        Queries each variant individually (CADD API doesn't support
        true batch lookups). Results cached per-variant.

        Parameters
        ----------
        variants
            List of (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[Optional[float]]
            CADD Phred scores in input order. None where unavailable.
        """
        return [
            self.lookup_score(chrom, pos, ref, alt) for chrom, pos, ref, alt in variants
        ]

    def _parse_response(self, response: object, ref: str, alt: str) -> Optional[float]:
        """Extract CADD Phred score matching the specific ref/alt pair.

        The CADD API returns a list-of-lists where the first row is headers
        and subsequent rows are scores for all possible substitutions at
        the queried position. We filter by matching Ref and Alt columns.

        Response format:
            [["Chrom", "Pos", "Ref", "Alt", ..., "PHRED"], [...], ...]
        """
        try:
            data = response.json()  # type: ignore[attr-defined]
        except (ValueError, AttributeError):
            return None

        if not isinstance(data, list) or len(data) < 2:
            return None

        # First row is the header
        header = data[0]
        if not isinstance(header, list):
            return None

        # Find column indices
        try:
            ref_idx = header.index("Ref")
            alt_idx = header.index("Alt")
            phred_idx = header.index("PHRED")
        except ValueError:
            # Fallback column positions if header names differ
            # Standard CADD output: Chrom(0), Pos(1), Ref(2), Alt(3), ..., PHRED(-1)
            if len(header) >= 6:
                ref_idx = 2
                alt_idx = 3
                phred_idx = len(header) - 1
            else:
                return None

        # Search data rows for matching ref/alt
        for row in data[1:]:
            if not isinstance(row, list) or len(row) <= phred_idx:
                continue
            if row[ref_idx] == ref and row[alt_idx] == alt:
                try:
                    return float(row[phred_idx])
                except (ValueError, TypeError):
                    continue

        return None

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self._http.close()
