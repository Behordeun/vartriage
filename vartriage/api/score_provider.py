"""API-based pathogenicity score provider.

Provides CADD and SpliceAI scores via remote APIs, matching the same
lookup interface used by the local ScoreLoader. Implements the CADD
score source hierarchy: VEP plugin data is used first (already fetched
during annotation), standalone CADD API is a fallback for positions
where VEP lacked pre-computed scores.

REVEL has no public API. PP3 evidence from REVEL requires a local file.
"""

from __future__ import annotations

import logging
from typing import Optional

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.cadd_client import CADDClient
from vartriage.api.config import APIConfig
from vartriage.api.spliceai_client import SpliceAIClient
from vartriage.api.vep_client import VEPAnnotation
from vartriage.models.variant import FunctionalConsequence
from vartriage.prioritization.score_loader import CoordinateKey

logger = logging.getLogger(__name__)


class APIScoreProvider:
    """Fetches CADD and SpliceAI scores via remote APIs.

    Implements the CADD score hierarchy:
    1. VEP plugin CADD (already retrieved during annotation, no extra call)
    2. Standalone CADD API (fallback for VEP-missing positions)

    SpliceAI queries are filtered to splice-relevant consequences
    to respect the strict rate limit (5 requests/minute).

    REVEL is not available via API. The provider returns empty dicts
    for REVEL lookups and logs a warning on first call.

    Parameters
    ----------
    config
        API configuration with rate limits and timeouts.
    cache
        Shared response cache (same instance as annotation engine).
    """

    def __init__(
        self, config: APIConfig, cache: Optional[ResponseCache] = None
    ) -> None:
        self._config = config
        self._revel_warned = False

        if cache is None:
            cache = ResponseCache(
                db_path=config.cache_path,
                default_ttl_days=config.cache_ttl_days,
            )
        self._cache = cache

        cadd_limiter = RateLimiter(
            tokens_per_second=config.cadd_rate_limit,
            service_name="cadd",
        )
        cadd_breaker = CircuitBreaker(service_name="cadd")

        spliceai_limiter = RateLimiter(
            tokens_per_second=config.spliceai_rate_limit,
            service_name="spliceai",
        )
        spliceai_breaker = CircuitBreaker(service_name="spliceai")

        self._cadd = CADDClient(
            rate_limiter=cadd_limiter,
            cache=self._cache,
            circuit_breaker=cadd_breaker,
            genome_build=config.genome_build,
            max_retries=config.max_retries,
            timeout=(config.connect_timeout, config.read_timeout),
            proxy_url=config.proxy_url,
        )

        self._spliceai = SpliceAIClient(
            rate_limiter=spliceai_limiter,
            cache=self._cache,
            circuit_breaker=spliceai_breaker,
            genome_build=config.genome_build,
            max_retries=config.max_retries,
            timeout=(config.connect_timeout, config.read_timeout),
            proxy_url=config.proxy_url,
        )

    def lookup_cadd_batch(
        self,
        keys: list[CoordinateKey],
        vep_annotations: Optional[list[Optional[VEPAnnotation]]] = None,
    ) -> list[Optional[float]]:
        """Look up CADD Phred scores with VEP-first hierarchy.

        For each variant:
        1. If VEP annotation includes cadd_phred, use it (no network call)
        2. Otherwise, query the standalone CADD API

        Parameters
        ----------
        keys
            Variant coordinate keys (chrom, pos, ref, alt).
        vep_annotations
            VEP results from the annotation step (same length as keys).
            If provided, CADD scores from VEP plugin are used first.

        Returns
        -------
        list[Optional[float]]
            CADD Phred scores in input order. None where unavailable.
        """
        results: list[Optional[float]] = []
        cadd_api_needed: list[tuple[int, CoordinateKey]] = []

        for i, key in enumerate(keys):
            # Check VEP plugin data first
            vep_cadd: Optional[float] = None
            if vep_annotations and i < len(vep_annotations):
                vep_ann = vep_annotations[i]
                if vep_ann is not None and vep_ann.cadd_phred is not None:
                    vep_cadd = vep_ann.cadd_phred

            if vep_cadd is not None:
                results.append(vep_cadd)
            else:
                results.append(None)  # placeholder
                cadd_api_needed.append((i, key))

        # Batch-fetch missing CADD scores from standalone API
        if cadd_api_needed:
            api_keys = [k for _, k in cadd_api_needed]
            api_scores = self._cadd.lookup_batch(api_keys)

            for idx, (original_idx, _) in enumerate(cadd_api_needed):
                if idx < len(api_scores) and api_scores[idx] is not None:
                    results[original_idx] = api_scores[idx]

        return results

    def lookup_spliceai_batch(
        self,
        keys: list[CoordinateKey],
        consequences: Optional[list[Optional[FunctionalConsequence]]] = None,
    ) -> list[Optional[float]]:
        """Look up SpliceAI delta scores with consequence-based filtering.

        Only queries splice-relevant variants (MISSENSE, SPLICE_SITE,
        SYNONYMOUS, IN_FRAME_*) to conserve the strict rate limit.

        Parameters
        ----------
        keys
            Variant coordinate keys (chrom, pos, ref, alt).
        consequences
            Pre-computed consequences for smart filtering. None queries all.

        Returns
        -------
        list[Optional[float]]
            Max SpliceAI delta scores in input order. None where unavailable.
        """
        variants_for_api = [(c, p, r, a) for c, p, r, a in keys]
        return self._spliceai.lookup_batch(variants_for_api, consequences)

    def lookup_revel_batch(self, keys: list[CoordinateKey]) -> list[Optional[float]]:
        """REVEL has no public API. Returns None for all variants.

        Logs a warning on first call to inform the user that PP3
        evidence requires a local REVEL file (--revel-scores).
        """
        if not self._revel_warned:
            logger.warning(
                "REVEL scores unavailable in API mode (no public API exists). "
                "PP3 evidence from REVEL requires --revel-scores with a local file."
            )
            self._revel_warned = True

        return [None] * len(keys)

    def close(self) -> None:
        """Release HTTP connections."""
        self._cadd.close()
        self._spliceai.close()
