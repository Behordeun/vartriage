"""gnomAD GraphQL API client for population allele frequency lookups.

Queries the gnomAD public API at gnomad.broadinstitute.org for
per-population allele frequencies. Used as a fallback when VEP's
colocated_variants field doesn't return gnomAD frequencies, or when
local gnomAD files are unavailable.

The API accepts one variant at a time (no batch support). Responses
are cached in the local SQLite database to avoid re-querying.
"""

from __future__ import annotations

import logging
from typing import Any

from vartriage import __version__
from vartriage.api._base import APIClientError, BaseAPIClient
from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._notation import _strip_chr_prefix
from vartriage.api._rate_limiter import DailyLimitExhausted, RateLimiter
from vartriage.models.variant import PopulationFrequencies

logger = logging.getLogger(__name__)

_GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

# GraphQL query for variant allele frequencies
_VARIANT_QUERY = """
query VariantFrequency($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    genome {
      ac
      an
      af
      populations {
        id
        ac
        an
      }
    }
    exome {
      ac
      an
      af
      populations {
        id
        ac
        an
      }
    }
  }
}
"""

# gnomAD population ID mapping to our PopulationFrequencies fields
_POP_MAP: dict[str, str] = {
    "afr": "afr",
    "amr": "amr",
    "ami": "amr",
    "asj": "asj",
    "eas": "eas",
    "fin": "fin",
    "nfe": "nfe",
    "sas": "sas",
    "mid": "nfe",  # Middle Eastern grouped with NFE in older versions
}


_VALID_DATASETS: frozenset[str] = frozenset(
    {
        "gnomad_r4",
        "gnomad_r3",
        "gnomad_r2_1",
    }
)

_VALID_SOURCES: frozenset[str] = frozenset(
    {
        "combined",
        "exome",
        "genome",
    }
)


class GnomADClient:
    """gnomAD GraphQL API client for allele frequency lookups.

    Parameters
    ----------
    rate_limiter
        Token bucket for gnomAD rate limiting.
    cache
        Response cache for deduplication.
    circuit_breaker
        Circuit breaker for gnomAD endpoint.
    dataset
        gnomAD dataset version: "gnomad_r4", "gnomad_r3", or "gnomad_r2_1".
    prefer_source
        Which data source to prefer: "exome", "genome", or "combined".
        "combined" picks whichever source has higher allele number (more samples).
    max_retries
        Retry attempts for transient failures.
    timeout
        (connect, read) timeouts in seconds.
    """

    def __init__(
        self,
        rate_limiter: RateLimiter,
        cache: ResponseCache,
        circuit_breaker: CircuitBreaker,
        dataset: str = "gnomad_r4",
        prefer_source: str = "combined",
        max_retries: int = 3,
        timeout: tuple[float, float] = (10.0, 30.0),
    ) -> None:
        if dataset not in _VALID_DATASETS:
            raise ValueError(
                f"Invalid dataset '{dataset}'. "
                f"Must be one of: {sorted(_VALID_DATASETS)}"
            )
        if prefer_source not in _VALID_SOURCES:
            raise ValueError(
                f"Invalid prefer_source '{prefer_source}'. "
                f"Must be one of: {sorted(_VALID_SOURCES)}"
            )

        self._dataset = dataset
        self._prefer_source = prefer_source
        self._cache = cache

        self._http = BaseAPIClient(
            base_url=_GNOMAD_API_URL,
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
            service_name="gnomad",
            timeout=timeout,
            max_retries=max_retries,
            user_agent=f"vartriage/{__version__} (https://github.com/Behordeun/vartriage)",
        )

    def lookup_frequency(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> PopulationFrequencies | None:
        """Look up population allele frequencies for a single variant.

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
        Optional[PopulationFrequencies]
            Per-population frequencies, or None if not found.
        """
        chrom_clean = _strip_chr_prefix(chrom)
        variant_id = f"{chrom_clean}-{pos}-{ref}-{alt}"

        # Check cache (use normalized chrom to avoid duplicate entries for chr1 vs 1)
        cache_key = ResponseCache.build_key(
            "gnomad", self._dataset, chrom_clean, pos, ref, alt
        )
        cached = self._cache.get(cache_key)
        if cached is not None:
            return self._parse_cached(cached)

        # Query API
        try:
            response = self._http.request(
                "POST",
                "",
                json_body={
                    "query": _VARIANT_QUERY,
                    "variables": {
                        "variantId": variant_id,
                        "dataset": self._dataset,
                    },
                },
                headers={"Content-Type": "application/json"},
            )
        except (APIClientError, CircuitBreakerOpen, DailyLimitExhausted) as exc:
            logger.debug("gnomAD lookup failed for %s: %s", variant_id, exc)
            return None

        # Parse response
        try:
            data = response.json()
        except (ValueError, AttributeError):
            return None

        # GraphQL can return errors with data.variant == null; don't cache
        # these as genuine misses since they may be transient failures
        errors = data.get("errors")
        if errors:
            logger.debug("gnomAD GraphQL errors for %s: %s", variant_id, errors)
            return None

        variant_data = data.get("data", {}).get("variant")
        if variant_data is None:
            # Variant genuinely not found in gnomAD (no errors) — cache the miss
            self._cache.put(
                key=cache_key,
                value={"not_found": True},
                source="gnomad",
                genome_build=self._dataset,
            )
            return None

        result = self._parse_variant_data(variant_data)

        # Cache the result
        self._cache.put(
            key=cache_key,
            value=self._frequencies_to_dict(result) if result else {"not_found": True},
            source="gnomad",
            genome_build=self._dataset,
        )

        return result

    def lookup_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[PopulationFrequencies | None]:
        """Look up frequencies for multiple variants (sequential, cached).

        gnomAD API doesn't support batch queries, so this iterates
        through variants individually. The cache prevents re-querying
        on subsequent runs.
        """
        return [
            self.lookup_frequency(chrom, pos, ref, alt)
            for chrom, pos, ref, alt in variants
        ]

    def _parse_variant_data(
        self, variant_data: dict[str, Any]
    ) -> PopulationFrequencies | None:
        """Extract population frequencies from gnomAD response."""
        exome = variant_data.get("exome")
        genome = variant_data.get("genome")

        # Pick the source with more data
        source = self._select_source(exome, genome)
        if source is None:
            return None

        global_af = source.get("af")
        populations = source.get("populations", [])

        pop_freqs: dict[str, float | None] = {
            "afr": None,
            "amr": None,
            "asj": None,
            "eas": None,
            "fin": None,
            "nfe": None,
            "sas": None,
        }

        for pop in populations:
            pop_id = pop.get("id", "").lower()
            an = pop.get("an", 0)
            ac = pop.get("ac", 0)

            mapped = _POP_MAP.get(pop_id)
            if mapped is None or an == 0:
                continue

            af = ac / an
            # Keep the higher frequency if multiple gnomAD pops map to same field
            if pop_freqs[mapped] is None or af > pop_freqs[mapped]:
                pop_freqs[mapped] = af

        return PopulationFrequencies(
            global_af=global_af,
            afr=pop_freqs["afr"],
            amr=pop_freqs["amr"],
            asj=pop_freqs["asj"],
            eas=pop_freqs["eas"],
            fin=pop_freqs["fin"],
            nfe=pop_freqs["nfe"],
            sas=pop_freqs["sas"],
        )

    def _select_source(
        self, exome: dict[str, Any] | None, genome: dict[str, Any] | None
    ) -> dict[str, Any] | None:
        """Pick the gnomAD data source based on preference."""
        if self._prefer_source == "exome":
            return exome or genome
        if self._prefer_source == "genome":
            return genome or exome

        # "combined": pick whichever has more alleles (higher AN)
        exome_an = exome.get("an", 0) if exome else 0
        genome_an = genome.get("an", 0) if genome else 0

        if exome_an >= genome_an and exome:
            return exome
        if genome:
            return genome
        return exome

    def _parse_cached(self, cached: dict[str, Any]) -> PopulationFrequencies | None:
        """Reconstruct PopulationFrequencies from cached dict."""
        if cached.get("not_found"):
            return None
        return PopulationFrequencies(
            global_af=cached.get("global_af"),
            afr=cached.get("afr"),
            amr=cached.get("amr"),
            asj=cached.get("asj"),
            eas=cached.get("eas"),
            fin=cached.get("fin"),
            nfe=cached.get("nfe"),
            sas=cached.get("sas"),
        )

    @staticmethod
    def _frequencies_to_dict(freq: PopulationFrequencies) -> dict[str, float | None]:
        """Serialize PopulationFrequencies for cache storage."""
        return {
            "global_af": freq.global_af,
            "afr": freq.afr,
            "amr": freq.amr,
            "asj": freq.asj,
            "eas": freq.eas,
            "fin": freq.fin,
            "nfe": freq.nfe,
            "sas": freq.sas,
        }

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self._http.close()
