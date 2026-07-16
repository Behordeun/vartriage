"""Ensembl VEP REST API client.

Queries the Variant Effect Predictor for consequence, gene name,
CADD score (plugin), and gnomAD allele frequency (colocated_variants)
in batches of up to 200 variants per POST request.

This is the primary annotation source in API mode. A single VEP call
returns data that would require three separate local files (GTF, gnomAD,
CADD) to replicate.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Optional

from vartriage.api._base import APIClientError, BaseAPIClient
from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker, CircuitBreakerOpen
from vartriage.api._consequence_map import (map_vep_most_severe,
                                            most_severe_consequence)
from vartriage.api._notation import vcf_to_vep_notation
from vartriage.api._rate_limiter import DailyLimitExhausted, RateLimiter
from vartriage.models.variant import FunctionalConsequence

logger = logging.getLogger(__name__)

_VEP_GRCH38_URL = "https://rest.ensembl.org"
_VEP_GRCH37_URL = "https://grch37.rest.ensembl.org"

# Response header containing Ensembl release version
_RELEASE_HEADER = "x-ensembl-release"

# VEP batch size ceiling (Ensembl enforces this server-side)
MAX_VEP_BATCH_SIZE = 200


@dataclass(frozen=True)
class VEPAnnotation:
    """Parsed annotation from a single VEP response entry."""

    consequence: FunctionalConsequence
    gene_name: Optional[str]
    allele_frequency: Optional[float]
    cadd_phred: Optional[float]
    transcript_id: Optional[str]
    hgvsc: Optional[str]
    hgvsp: Optional[str]
    ensembl_release: Optional[str]


class VEPClient:
    """Ensembl VEP REST API client with batch annotation.

    Parameters
    ----------
    rate_limiter
        Token bucket for VEP request throttling (15 req/sec default).
    cache
        Response cache for deduplication.
    circuit_breaker
        Circuit breaker for Ensembl endpoint.
    genome_build
        Target assembly: "grch38" or "grch37".
    batch_size
        Variants per POST request (max 200).
    max_retries
        Retry attempts passed to BaseAPIClient.
    timeout
        (connect, read) timeouts in seconds.
    user_agent
        User-Agent header value.
    proxy_url
        HTTP/HTTPS proxy URL.
    preferred_frequency_source
        Which gnomAD frequency field to extract: "gnomad_exome" or "gnomad_genome".
    """

    def __init__(
        self,
        rate_limiter: RateLimiter,
        cache: ResponseCache,
        circuit_breaker: CircuitBreaker,
        genome_build: str = "grch38",
        batch_size: int = 200,
        max_retries: int = 3,
        timeout: tuple[float, float] = (10.0, 30.0),
        user_agent: str = "vartriage/0.7.0 (https://github.com/Behordeun/vartriage)",
        proxy_url: Optional[str] = None,
        preferred_frequency_source: str = "gnomad_exome",
    ) -> None:
        base_url = _VEP_GRCH38_URL if genome_build == "grch38" else _VEP_GRCH37_URL
        self._batch_size = min(batch_size, MAX_VEP_BATCH_SIZE)
        self._genome_build = genome_build
        self._cache = cache
        self._preferred_freq = preferred_frequency_source
        self._ensembl_release: Optional[str] = None

        self._http = BaseAPIClient(
            base_url=base_url,
            rate_limiter=rate_limiter,
            cache=cache,
            circuit_breaker=circuit_breaker,
            service_name="vep",
            timeout=timeout,
            max_retries=max_retries,
            user_agent=user_agent,
            proxy_url=proxy_url,
        )

    @property
    def ensembl_release(self) -> Optional[str]:
        """Ensembl release version from the most recent response."""
        return self._ensembl_release

    def annotate_batch(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[VEPAnnotation]]:
        """Annotate a batch of variants via VEP POST endpoint.

        Splits into sub-batches of self._batch_size (max 200) and sends
        sequential requests. On batch failure, applies progressive recovery:
        split in half, retry sub-batches, then fall back to individual queries.

        Parameters
        ----------
        variants
            List of (chrom, pos, ref, alt) tuples.

        Returns
        -------
        list[Optional[VEPAnnotation]]
            Annotations in input order. None for variants that failed
            annotation (cache miss + API error after all recovery attempts).
        """
        results: list[Optional[VEPAnnotation]] = [None] * len(variants)

        for start in range(0, len(variants), self._batch_size):
            chunk = variants[start : start + self._batch_size]
            chunk_results = self._annotate_with_recovery(chunk)
            for i, annotation in enumerate(chunk_results):
                results[start + i] = annotation

        return results

    def _annotate_with_recovery(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[VEPAnnotation]]:
        """Annotate a chunk with progressive failure recovery.

        Recovery strategy:
        1. Try full chunk
        2. If failed: split in half, try each half
        3. If half fails: try individual variants from that half
        """
        result = self._annotate_chunk(variants)

        # Check if we got any results (all None means total failure)
        if any(r is not None for r in result):
            return result

        # Full chunk failed. If chunk is small enough, no point splitting further.
        if len(variants) <= 1:
            return result

        logger.info(
            "VEP batch failed for %d variants, splitting for recovery",
            len(variants),
        )

        # Split in half and retry
        mid = len(variants) // 2
        left_variants = variants[:mid]
        right_variants = variants[mid:]

        left_results = self._annotate_chunk(left_variants)
        right_results = self._annotate_chunk(right_variants)

        # If a half still failed entirely, try individual queries
        if not any(r is not None for r in left_results):
            left_results = self._annotate_individually(left_variants)

        if not any(r is not None for r in right_results):
            right_results = self._annotate_individually(right_variants)

        return left_results + right_results

    def _annotate_individually(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[VEPAnnotation]]:
        """Last-resort: query each variant individually.

        Slow but maximizes the number of successfully annotated variants
        when batch requests are failing (server-side issue with specific
        variants can poison an entire batch).
        """
        results: list[Optional[VEPAnnotation]] = []
        for variant in variants:
            single_result = self._annotate_chunk([variant])
            results.append(single_result[0] if single_result else None)
        return results

    def _annotate_chunk(
        self, variants: list[tuple[str, int, str, str]]
    ) -> list[Optional[VEPAnnotation]]:
        """Send a single VEP POST request for a chunk (<=200 variants).

        Checks the cache first per-variant. Only sends uncached variants
        to the API, then merges cached + fresh results in input order.
        """
        # Separate cached from uncached
        results: list[Optional[VEPAnnotation]] = [None] * len(variants)
        uncached_indices: list[int] = []
        uncached_variants: list[tuple[str, int, str, str]] = []

        for i, (chrom, pos, ref, alt) in enumerate(variants):
            cache_key = ResponseCache.build_key(
                "vep", self._genome_build, chrom, pos, ref, alt
            )
            cached = self._cache.get(cache_key)
            if cached is not None:
                results[i] = self._parse_cached_entry(cached)
            else:
                uncached_indices.append(i)
                uncached_variants.append((chrom, pos, ref, alt))

        if not uncached_variants:
            return results

        # Build VEP notation for uncached variants
        vep_notations = [
            vcf_to_vep_notation(chrom, pos, ref, alt)
            for chrom, pos, ref, alt in uncached_variants
        ]

        # POST to VEP
        try:
            response = self._http.request(
                "POST",
                "/vep/human/region",
                json_body={
                    "variants": vep_notations,
                    "CADD": 1,
                    "gnomADe": 1,
                    "hgvs": 1,
                },
                headers={"Content-Type": "application/json"},
            )
        except (APIClientError, CircuitBreakerOpen, DailyLimitExhausted) as exc:
            logger.warning("VEP batch request failed: %s", exc)
            return results

        # Extract Ensembl release version from response headers
        release = response.headers.get(_RELEASE_HEADER)
        if release:
            self._ensembl_release = release

        # Parse response
        try:
            response_data = response.json()
        except (ValueError, AttributeError):
            logger.warning("VEP returned non-JSON response")
            return results

        if not isinstance(response_data, list):
            logger.warning("VEP response is not a list: %s", type(response_data))
            return results

        # VEP returns results keyed by input string. Build a lookup.
        parsed_by_input: dict[str, Any] = {}
        for entry in response_data:
            input_str = entry.get("input")
            if input_str:
                parsed_by_input[input_str] = entry

        # Match responses back to input variants
        self._store_results(
            parsed_by_input, vep_notations, uncached_indices,
            uncached_variants, results, release,
        )
        return results

    def _store_results(
        self,
        parsed_by_input: dict[str, Any],
        vep_notations: list[str],
        uncached_indices: list[int],
        uncached_variants: list[tuple[str, int, str, str]],
        results: list[Optional[VEPAnnotation]],
        release: Optional[str],
    ) -> None:
        """Parse VEP responses, populate results, and cache entries."""
        for idx_in_uncached, notation in enumerate(vep_notations):
            entry = parsed_by_input.get(notation)
            if entry is None:
                continue

            results[uncached_indices[idx_in_uncached]] = self._parse_vep_entry(entry)

            chrom, pos, ref, alt = uncached_variants[idx_in_uncached]
            cache_key = ResponseCache.build_key(
                "vep", self._genome_build, chrom, pos, ref, alt
            )
            self._cache.put(
                key=cache_key,
                value=entry,
                source="vep",
                genome_build=self._genome_build,
                source_version=f"Ensembl {release}" if release else None,
            )

    def _parse_vep_entry(self, entry: dict[str, Any]) -> VEPAnnotation:
        """Extract structured annotation from a single VEP response object."""
        # Consequence: use top-level most_severe_consequence
        most_severe_str = entry.get("most_severe_consequence")
        consequence = map_vep_most_severe(most_severe_str)

        # Gene name + transcript details from transcript_consequences
        gene_name: Optional[str] = None
        transcript_id: Optional[str] = None
        hgvsc: Optional[str] = None
        hgvsp: Optional[str] = None
        cadd_phred: Optional[float] = None

        transcript_consequences = entry.get("transcript_consequences", [])
        if transcript_consequences:
            # Pick the canonical transcript, or the first one
            canonical = self._pick_canonical_transcript(transcript_consequences)
            gene_name = canonical.get("gene_symbol")
            transcript_id = canonical.get("transcript_id")
            hgvsc = canonical.get("hgvsc")
            hgvsp = canonical.get("hgvsp")

            # CADD from plugin data (may appear at transcript level)
            cadd_raw = canonical.get("cadd_phred")
            if cadd_raw is not None:
                try:
                    cadd_phred = float(cadd_raw)
                except (ValueError, TypeError):
                    pass

        # gnomAD frequency from colocated_variants
        allele_frequency = self._extract_frequency(entry)

        return VEPAnnotation(
            consequence=consequence,
            gene_name=gene_name,
            allele_frequency=allele_frequency,
            cadd_phred=cadd_phred,
            transcript_id=transcript_id,
            hgvsc=hgvsc,
            hgvsp=hgvsp,
            ensembl_release=self._ensembl_release,
        )

    def _parse_cached_entry(self, cached: dict[str, Any]) -> VEPAnnotation:
        """Re-parse a cached VEP response dict into VEPAnnotation."""
        return self._parse_vep_entry(cached)

    def _pick_canonical_transcript(
        self, transcripts: list[dict[str, Any]]
    ) -> dict[str, Any]:
        """Select the canonical transcript from VEP transcript_consequences.

        VEP marks one transcript per gene as canonical (canonical=1).
        Falls back to the first transcript with a gene_symbol if no
        canonical flag is set.
        """
        for tc in transcripts:
            if tc.get("canonical") == 1:
                return tc

        # Fallback: first transcript with a gene symbol
        for tc in transcripts:
            if tc.get("gene_symbol"):
                return tc

        return transcripts[0] if transcripts else {}

    def _extract_frequency(self, entry: dict[str, Any]) -> Optional[float]:
        """Extract gnomAD allele frequency from colocated_variants.

        Priority:
        1. gnomAD exomes (gnomade_af) when preferred_frequency_source is gnomad_exome
        2. gnomAD genomes (gnomadg_af) as fallback
        3. None if no frequency data
        """
        colocated = entry.get("colocated_variants", [])
        if not colocated:
            return None

        primary_key = (
            "gnomade_af" if self._preferred_freq == "gnomad_exome" else "gnomadg_af"
        )
        fallback_key = (
            "gnomadg_af" if self._preferred_freq == "gnomad_exome" else "gnomade_af"
        )

        for variant in colocated:
            frequencies = variant.get("frequencies", {})
            for _allele, sources in frequencies.items():
                freq = self._freq_from_sources(sources, primary_key, fallback_key)
                if freq is not None:
                    return freq

        return None

    @staticmethod
    def _freq_from_sources(
        sources: object, primary_key: str, fallback_key: str
    ) -> Optional[float]:
        """Try to extract a float frequency from a sources dict."""
        if not isinstance(sources, dict):
            return None
        for key in (primary_key, fallback_key):
            if key in sources:
                try:
                    return float(sources[key])
                except (ValueError, TypeError):
                    continue
        return None

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self._http.close()
