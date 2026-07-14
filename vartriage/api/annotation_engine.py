"""API-based annotation engine composing VEP and ClinVar into AnnotatedVariants.

Same interface as the local AnnotationEngine: accepts an iterator of
raw Variants and yields AnnotatedVariants enriched with consequence,
gene name, allele frequency, and ClinVar assertions. All data comes
from remote APIs rather than local files.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from typing import Iterator, Optional

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.clinvar_client import ClinVarClient
from vartriage.api.config import APIConfig
from vartriage.api.vep_client import VEPAnnotation, VEPClient
from vartriage.models.variant import (AnnotatedVariant, ClinVarAssertion,
                                      FunctionalConsequence, Variant)
from vartriage.models.warnings import MissingDataWarning

logger = logging.getLogger(__name__)


class APIAnnotationEngine:
    """Annotates variants via Ensembl VEP and ClinVar E-utilities.

    Processes variants in batches sized for VEP (default 200). Within
    each batch, VEP provides consequence + gene + frequency, while
    ClinVar provides clinical significance. ClinVar lookups run
    concurrently with VEP response parsing via a thread pool.

    Parameters
    ----------
    config
        API configuration with rate limits, timeouts, and build info.
    """

    def __init__(self, config: APIConfig) -> None:
        self._config = config
        self._batch_size = config.vep_batch_size
        self._warnings: list[MissingDataWarning] = []

        # Shared cache and per-service resilience components
        self._cache = ResponseCache(
            db_path=config.cache_path,
            default_ttl_days=config.cache_ttl_days,
        )

        vep_limiter = RateLimiter(
            tokens_per_second=config.vep_rate_limit,
            daily_limit=config.vep_daily_limit,
            service_name="vep",
        )
        vep_breaker = CircuitBreaker(service_name="vep")

        clinvar_limiter = RateLimiter(
            tokens_per_second=config.clinvar_rate_limit,
            service_name="clinvar",
        )
        clinvar_breaker = CircuitBreaker(service_name="clinvar")

        self._vep = VEPClient(
            rate_limiter=vep_limiter,
            cache=self._cache,
            circuit_breaker=vep_breaker,
            genome_build=config.genome_build,
            batch_size=config.vep_batch_size,
            max_retries=config.max_retries,
            timeout=(config.connect_timeout, config.read_timeout),
            proxy_url=config.proxy_url,
        )

        self._clinvar = ClinVarClient(
            rate_limiter=clinvar_limiter,
            cache=self._cache,
            circuit_breaker=clinvar_breaker,
            ncbi_api_key=config.ncbi_api_key,
            genome_build=config.genome_build,
            max_retries=config.max_retries,
            timeout=(config.connect_timeout, config.read_timeout),
            proxy_url=config.proxy_url,
        )

    @property
    def warnings(self) -> list[MissingDataWarning]:
        """Accumulated missing-data warnings for the limitations section."""
        return self._warnings

    def annotate(self, variants: Iterator[Variant]) -> Iterator[AnnotatedVariant]:
        """Annotate variants via VEP + ClinVar API calls.

        Consumes the input iterator in batches, queries VEP for
        consequence/frequency and ClinVar for significance, then
        yields composed AnnotatedVariant records.

        Parameters
        ----------
        variants
            Input stream of raw variants from VCF parsing.

        Yields
        ------
        AnnotatedVariant
            Variants enriched with consequence, gene name, frequency,
            and ClinVar assertion from remote APIs.
        """
        while True:
            batch = list(islice(variants, self._batch_size))
            if not batch:
                break
            yield from self._annotate_batch(batch)

    def _annotate_batch(self, batch: list[Variant]) -> list[AnnotatedVariant]:
        """Process a single batch through VEP + ClinVar."""
        keys = [(v.chrom, v.pos, v.ref, v.alt) for v in batch]

        # VEP: consequence + gene + frequency (+ optional CADD)
        vep_results = self._vep.annotate_batch(keys)

        # ClinVar: clinical significance (concurrent with processing)
        with ThreadPoolExecutor(max_workers=1) as executor:
            clinvar_future = executor.submit(self._clinvar.lookup_batch, keys)
            clinvar_results = clinvar_future.result()

        # Compose AnnotatedVariant records
        results: list[AnnotatedVariant] = []
        for i, variant in enumerate(batch):
            vep_ann: Optional[VEPAnnotation] = (
                vep_results[i] if i < len(vep_results) else None
            )
            clinvar_assertion: Optional[ClinVarAssertion] = (
                clinvar_results[i] if i < len(clinvar_results) else None
            )

            consequence = FunctionalConsequence.INTERGENIC
            gene_name: Optional[str] = None
            allele_frequency: Optional[float] = None

            if vep_ann is not None:
                consequence = vep_ann.consequence
                gene_name = vep_ann.gene_name
                allele_frequency = vep_ann.allele_frequency

            frequency_unknown = allele_frequency is None
            clinvar_unknown = clinvar_assertion is None

            # Track missing data for the limitations section
            if vep_ann is None:
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="VEP",
                        reason="api_error",
                    )
                )

            if frequency_unknown and vep_ann is not None:
                # VEP responded but had no gnomAD data for this position
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="gnomAD",
                        reason="not_found",
                    )
                )

            if clinvar_unknown:
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="ClinVar",
                        reason="not_found",
                    )
                )

            results.append(
                AnnotatedVariant(
                    variant=variant,
                    consequence=consequence,
                    allele_frequency=allele_frequency,
                    clinvar_assertion=clinvar_assertion,
                    frequency_unknown=frequency_unknown,
                    clinvar_unknown=clinvar_unknown,
                    gene_name=gene_name,
                )
            )

        return results

    def close(self) -> None:
        """Release HTTP connections and close the cache."""
        self._vep.close()
        self._clinvar.close()
        self._cache.close()
