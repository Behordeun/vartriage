"""NCBI ClinVar E-utilities client.

Queries ClinVar for clinical significance assertions using the
esearch + esummary API pattern. Supports NCBI API key for higher
rate limits (10 req/sec vs 3 without key).

Query strategy:
  1. esearch: find ClinVar variation IDs matching (chrom, pos, ref, alt)
  2. esummary: fetch clinical significance for matched variation IDs

This avoids the heavier efetch XML parsing while still getting the
fields we need (clinical significance, review status).
"""

from __future__ import annotations

import logging
from typing import Optional

from vartriage.api._base import APIClientError, BaseAPIClient
from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._rate_limiter import DailyLimitExhausted, RateLimiter
from vartriage.models.variant import ClinVarAssertion

logger = logging.getLogger(__name__)

_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ClinVar clinical significance strings to our enum.
# ClinVar uses varied casing and separators; normalize to lowercase for matching.
_SIGNIFICANCE_MAP: dict[str, ClinVarAssertion] = {
    "pathogenic": ClinVarAssertion.PATHOGENIC,
    "likely pathogenic": ClinVarAssertion.LIKELY_PATHOGENIC,
    "pathogenic/likely pathogenic": ClinVarAssertion.LIKELY_PATHOGENIC,
    "uncertain significance": ClinVarAssertion.VUS,
    "likely benign": ClinVarAssertion.LIKELY_BENIGN,
    "benign": ClinVarAssertion.BENIGN,
    "benign/likely benign": ClinVarAssertion.BENIGN,
}

# Review status hierarchy (higher = more authoritative)
_REVIEW_STATUS_RANK: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
}


class ClinVarClient:
    """NCBI ClinVar E-utilities client for clinical significance lookups.

    Parameters
    ----------
    rate_limiter
        Token bucket for NCBI rate limiting.
    cache
        Response cache for deduplication.
    circuit_breaker
        Circuit breaker for NCBI endpoint.
    ncbi_api_key
        NCBI API key for higher rate limits. None uses the anonymous tier.
    genome_build
        Target assembly for coordinate queries.
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
        ncbi_api_key: Optional[str] = None,
        genome_build: str = "grch38",
        max_retries: int = 3,
        timeout: tuple[float, float] = (10.0, 30.0),
        user_agent: str = "vartriage/0.7.0 (https://github.com/Behordeun/vartriage)",
        proxy_url: Optional[str] = None,
    ) -> None:
        self._api_key = ncbi_api_key
        self._genome_build = genome_build
        self._cache = cache

        self._http = BaseAPIClient(
            base_url=_EUTILS_BASE,
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
            service_name="clinvar",
            timeout=timeout,
            max_retries=max_retries,
            user_agent=user_agent,
            proxy_url=proxy_url,
        )

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[ClinVarAssertion]]:
        """Look up ClinVar clinical significance for a batch of variants.

        Queries each variant individually (ClinVar esearch doesn't support
        true batch coordinate lookups). Results are cached per-variant.

        Parameters
        ----------
        variants
            List of (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[Optional[ClinVarAssertion]]
            Clinical significance in input order. None when no ClinVar
            entry exists or the query failed.
        """
        results: list[Optional[ClinVarAssertion]] = []

        for chrom, pos, ref, alt in variants:
            assertion = self._lookup_single(chrom, pos, ref, alt)
            results.append(assertion)

        return results

    def _lookup_single(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Optional[ClinVarAssertion]:
        """Query ClinVar for a single variant, checking cache first."""
        cache_key = ResponseCache.build_key(
            "clinvar", self._genome_build, chrom, pos, ref, alt
        )
        cached = self._cache.get(cache_key)
        if cached is not None:
            return self._parse_cached(cached)

        # esearch to find variation IDs
        variation_ids = self._esearch(chrom, pos, ref, alt)
        if not variation_ids:
            # Cache the miss to avoid repeated lookups
            self._cache.put(
                key=cache_key,
                value={"clinical_significance": None, "review_status": None},
                source="clinvar",
                genome_build=self._genome_build,
            )
            return None

        # esummary to get clinical significance
        assertion, raw_response = self._esummary(variation_ids)

        # Cache the result
        self._cache.put(
            key=cache_key,
            value=raw_response,
            source="clinvar",
            genome_build=self._genome_build,
            source_version=str(raw_response.get("last_evaluated") or ""),
            ttl_days=30,  # ClinVar updates monthly
        )

        return assertion

    def _esearch(self, chrom: str, pos: int, ref: str, alt: str) -> list[str]:
        """Search ClinVar for variation IDs matching coordinates.

        Query format: chrom[Chromosome] AND pos[Base Position] AND
        ref/alt alleles via variant name search.
        """
        # Normalize chromosome: strip 'chr' prefix, ClinVar uses numeric
        chrom_clean = chrom.replace("chr", "").replace("Chr", "")

        # Build the search query
        # ClinVar search supports: "CHR:POS REF>ALT" notation
        query = (
            f"{chrom_clean}[Chromosome] AND {pos}[Base Position for Assembly GRCh38]"
        )

        params: dict[str, str] = {
            "db": "clinvar",
            "term": query,
            "retmode": "json",
            "retmax": "5",
            "tool": "vartriage",
            "email": "vartriage@github.com",
        }
        if self._api_key:
            params["api_key"] = self._api_key

        try:
            response = self._http.request("GET", "/esearch.fcgi", params=params)
        except (APIClientError, CircuitBreakerOpen, DailyLimitExhausted) as exc:
            logger.warning("ClinVar esearch failed: %s", exc)
            return []

        try:
            data = response.json()
        except (ValueError, AttributeError):
            return []

        esearch_result = data.get("esearchresult", {})
        id_list: list[str] = esearch_result.get("idlist", [])

        return id_list[:5]  # Limit to top 5 matches

    def _esummary(
        self, variation_ids: list[str]
    ) -> tuple[Optional[ClinVarAssertion], dict[str, object]]:
        """Fetch clinical significance via esummary for variation IDs.

        When multiple IDs are returned, selects the one with the highest
        review status.
        """
        ids_str = ",".join(variation_ids)

        params: dict[str, str] = {
            "db": "clinvar",
            "id": ids_str,
            "retmode": "json",
            "tool": "vartriage",
            "email": "vartriage@github.com",
        }
        if self._api_key:
            params["api_key"] = self._api_key

        try:
            response = self._http.request("GET", "/esummary.fcgi", params=params)
        except (APIClientError, CircuitBreakerOpen, DailyLimitExhausted) as exc:
            logger.warning("ClinVar esummary failed: %s", exc)
            return None, {"clinical_significance": None, "review_status": None}

        try:
            data = response.json()
        except (ValueError, AttributeError):
            return None, {"clinical_significance": None, "review_status": None}

        result = data.get("result", {})

        # Pick the entry with highest review status
        best_assertion: Optional[ClinVarAssertion] = None
        best_rank = -1
        best_significance: Optional[str] = None
        best_review: Optional[str] = None
        last_evaluated: Optional[str] = None

        for var_id in variation_ids:
            entry = result.get(var_id, {})
            if not isinstance(entry, dict):
                continue

            clinical_sig = entry.get("clinical_significance", {})
            if isinstance(clinical_sig, dict):
                description = clinical_sig.get("description", "")
                review_status = clinical_sig.get("review_status", "")
                last_evaluated = clinical_sig.get("last_evaluated", last_evaluated)
            else:
                # Older API format: flat string
                description = str(clinical_sig) if clinical_sig else ""
                review_status = entry.get("review_status", "")

            rank = _REVIEW_STATUS_RANK.get(str(review_status).lower(), 0)
            assertion = _map_significance(str(description))

            if rank > best_rank or (
                rank == best_rank and assertion is not None and best_assertion is None
            ):
                best_rank = rank
                best_assertion = assertion
                best_significance = str(description) if description else None
                best_review = str(review_status) if review_status else None

        raw_response: dict[str, object] = {
            "clinical_significance": best_significance,
            "review_status": best_review,
            "last_evaluated": last_evaluated,
        }

        return best_assertion, raw_response

    def _parse_cached(self, cached: dict[str, object]) -> Optional[ClinVarAssertion]:
        """Re-parse a cached ClinVar response into ClinVarAssertion."""
        sig = cached.get("clinical_significance")
        if sig is None:
            return None
        return _map_significance(str(sig))

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self._http.close()


def _map_significance(description: str) -> Optional[ClinVarAssertion]:
    """Map a ClinVar clinical significance string to our enum.

    Handles varied casing, compound assertions (e.g., "Pathogenic/Likely pathogenic"),
    and conflicting interpretations by taking the most pathogenic assertion.
    """
    if not description:
        return None

    normalized = description.lower().strip()

    # Direct match
    result = _SIGNIFICANCE_MAP.get(normalized)
    if result is not None:
        return result

    # Check for compound or partial matches
    # "conflicting interpretations of pathogenicity" -> VUS
    if "conflicting" in normalized:
        return ClinVarAssertion.VUS

    # Check if any known term appears as a substring
    for key, assertion in _SIGNIFICANCE_MAP.items():
        if key in normalized:
            return assertion

    return None
