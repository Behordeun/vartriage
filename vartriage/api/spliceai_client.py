"""SpliceAI Lookup API client for splice-disruption delta scores.

Queries the Broad Institute SpliceAI Lookup service for pre-computed
delta scores. Extremely rate-limited (5 requests/minute per their docs),
so only splice-relevant variants are queried.

Smart filtering: only variants with MISSENSE or SPLICE_SITE consequence
are sent to SpliceAI. All others return None without a network call.
"""

from __future__ import annotations

import logging
from typing import Optional

from vartriage.api._base import APIClientError, BaseAPIClient
from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._notation import _strip_chr_prefix
from vartriage.api._rate_limiter import DailyLimitExhausted, RateLimiter
from vartriage.models.variant import FunctionalConsequence

logger = logging.getLogger(__name__)

_SPLICEAI_BASE_URL = "https://spliceailookup.broadinstitute.org"

# Consequences where SpliceAI is relevant
_SPLICE_RELEVANT_CONSEQUENCES: frozenset[FunctionalConsequence] = frozenset(
    {
        FunctionalConsequence.MISSENSE,
        FunctionalConsequence.SPLICE_SITE,
        FunctionalConsequence.SYNONYMOUS,
        FunctionalConsequence.IN_FRAME_INSERTION,
        FunctionalConsequence.IN_FRAME_DELETION,
    }
)


class SpliceAIClient:
    """Broad Institute SpliceAI Lookup API client.

    Parameters
    ----------
    rate_limiter
        Token bucket for SpliceAI rate limiting (0.08 req/sec = 5/min).
    cache
        Response cache for deduplication.
    circuit_breaker
        Circuit breaker for SpliceAI endpoint.
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
        self._hg_version = "38" if genome_build == "grch38" else "37"

        self._http = BaseAPIClient(
            base_url=_SPLICEAI_BASE_URL,
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
            service_name="spliceai",
            timeout=timeout,
            max_retries=max_retries,
            user_agent=user_agent,
            proxy_url=proxy_url,
        )

    def lookup_score(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        consequence: Optional[FunctionalConsequence] = None,
    ) -> Optional[float]:
        """Look up SpliceAI max delta score for a single variant.

        Skips the query entirely if the variant's consequence is not
        splice-relevant (saves rate limit tokens for variants where
        splice disruption is biologically implausible).

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
        consequence
            Pre-computed consequence for smart filtering. If None,
            the query proceeds unconditionally.

        Returns
        -------
        Optional[float]
            Maximum SpliceAI delta score (0.0-1.0), or None if not
            found, filtered out, or API error.
        """
        # Smart filter: skip irrelevant consequences
        if consequence is not None and consequence not in _SPLICE_RELEVANT_CONSEQUENCES:
            return None

        cache_key = ResponseCache.build_key(
            "spliceai", self._genome_build, chrom, pos, ref, alt
        )
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached.get("max_delta")

        chrom_clean = _strip_chr_prefix(chrom)

        # SpliceAI Lookup API: /spliceai/?hg=38&variant=chr-pos-ref-alt
        # Note: SpliceAI expects the chr prefix in the variant param
        variant_str = f"{chrom_clean}-{pos}-{ref}-{alt}"
        params = {
            "hg": self._hg_version,
            "variant": variant_str,
        }

        try:
            response = self._http.request("GET", "/spliceai/", params=params)
        except (APIClientError, CircuitBreakerOpen, DailyLimitExhausted) as exc:
            logger.warning("SpliceAI lookup failed for %s:%d: %s", chrom, pos, exc)
            return None

        max_delta = self._parse_response(response)

        # Cache the result
        self._cache.put(
            key=cache_key,
            value={"max_delta": max_delta},
            source="spliceai",
            genome_build=self._genome_build,
        )

        return max_delta

    def lookup_batch(
        self,
        variants: list[tuple[str, int, str, str]],
        consequences: Optional[list[Optional[FunctionalConsequence]]] = None,
    ) -> list[Optional[float]]:
        """Look up SpliceAI scores for multiple variants with smart filtering.

        Parameters
        ----------
        variants
            List of (chrom, pos, ref, alt) tuples.
        consequences
            Optional list of consequences (same length as variants) for
            smart filtering. None means query all variants.

        Returns
        -------
        list[Optional[float]]
            Max delta scores in input order. None for filtered/unavailable.
        """
        results: list[Optional[float]] = []

        for i, (chrom, pos, ref, alt) in enumerate(variants):
            consequence = consequences[i] if consequences else None
            score = self.lookup_score(chrom, pos, ref, alt, consequence)
            results.append(score)

        return results

    def _parse_response(self, response: object) -> Optional[float]:
        """Extract max delta score from SpliceAI Lookup response.

        Response format (JSON):
        {
            "scores": [
                {
                    "DS_AG": 0.01,  # Delta score acceptor gain
                    "DS_AL": 0.02,  # Delta score acceptor loss
                    "DS_DG": 0.05,  # Delta score donor gain
                    "DS_DL": 0.03,  # Delta score donor loss
                    ...
                }
            ]
        }

        We take the maximum across all four delta scores.
        """
        try:
            data = response.json()  # type: ignore[attr-defined]
        except (ValueError, AttributeError):
            return None

        if not isinstance(data, dict):
            return None

        scores_list = data.get("scores", [])
        if not scores_list:
            return None

        max_delta: float = 0.0

        for score_entry in scores_list:
            if not isinstance(score_entry, dict):
                continue
            for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
                val = score_entry.get(key)
                if val is not None:
                    try:
                        float_val = float(val)
                        if float_val > max_delta:
                            max_delta = float_val
                    except (ValueError, TypeError):
                        continue

        return max_delta if max_delta > 0.0 else None

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self._http.close()
