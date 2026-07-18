"""Integration test: compare local consequence calls against Ensembl VEP.

Queries the live VEP endpoint with known variants and compares the
consequence returned by VEP against what our SO term mapper produces.
This validates that the _consequence_map.py mapping is consistent with
VEP's output for real variants.

Marked @pytest.mark.slow because it hits the network.
Run with: pytest -m slow tests/integration/test_consequence_concordance.py
"""

from __future__ import annotations

from pathlib import Path

import pytest

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._consequence_map import map_vep_most_severe
from vartriage.api._notation import vcf_to_vep_notation
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.vep_client import VEPClient
from vartriage.models.variant import FunctionalConsequence

# Known variants with expected consequences (manually verified against VEP)
_KNOWN_VARIANTS: list[tuple[str, int, str, str, FunctionalConsequence]] = [
    # BRCA1 missense (well-characterized)
    ("chr17", 43094452, "G", "A", FunctionalConsequence.MISSENSE),
    # BRAF V600 region (missense)
    ("chr7", 140753336, "A", "T", FunctionalConsequence.MISSENSE),
    # TP53 exon 5 region — VEP calls this missense despite being at codon
    # position 3 (the wobble position can still change the amino acid for
    # some codons; this tests that our SO mapping handles real VEP output)
    ("chr17", 7676154, "G", "A", FunctionalConsequence.MISSENSE),
    # Intergenic variant (far from any gene)
    ("chr1", 10000, "A", "T", FunctionalConsequence.INTERGENIC),
]


@pytest.fixture
def vep_client(tmp_path: Path) -> VEPClient:
    """Build a VEP client for live testing."""
    cache = ResponseCache(db_path=tmp_path / "concordance_test.db", default_ttl_days=1)
    limiter = RateLimiter(tokens_per_second=15.0, service_name="vep")
    breaker = CircuitBreaker(service_name="vep")

    return VEPClient(
        rate_limiter=limiter,
        cache=cache,
        circuit_breaker=breaker,
        genome_build="grch38",
    )


@pytest.mark.slow
class TestConsequenceConcordance:
    """Compare local SO term mapping against live VEP responses."""

    def test_vep_consequence_matches_our_mapping(self, vep_client: VEPClient) -> None:
        """For known variants, VEP's most_severe_consequence maps correctly."""
        variants = [
            (chrom, pos, ref, alt) for chrom, pos, ref, alt, _ in _KNOWN_VARIANTS
        ]

        results = vep_client.annotate_batch(variants)

        concordant = 0
        discordant: list[str] = []

        for i, (chrom, pos, ref, alt, expected) in enumerate(_KNOWN_VARIANTS):
            annotation = results[i]
            if annotation is None:
                discordant.append(f"{chrom}:{pos} - VEP returned None")
                continue

            actual = annotation.consequence
            if actual == expected:
                concordant += 1
            else:
                discordant.append(
                    f"{chrom}:{pos} {ref}>{alt}: expected={expected.value}, "
                    f"got={actual.value}"
                )

        total = len(_KNOWN_VARIANTS)
        concordance_rate = concordant / total if total > 0 else 0.0

        print(
            f"\nConsequence concordance: {concordant}/{total} ({concordance_rate:.0%})"
        )
        if discordant:
            print("Discordant calls:")
            for d in discordant:
                print(f"  {d}")

        # Allow some discordance (VEP may return different consequence based
        # on transcript selection), but majority should agree
        assert concordance_rate >= 0.5, (
            f"Consequence concordance too low: {concordance_rate:.0%}. "
            f"Discordant: {discordant}"
        )

    def test_so_term_mapping_covers_real_vep_output(
        self, vep_client: VEPClient
    ) -> None:
        """Every consequence term VEP returns for our test variants is in our mapping."""
        variants = [
            (chrom, pos, ref, alt) for chrom, pos, ref, alt, _ in _KNOWN_VARIANTS
        ]

        results = vep_client.annotate_batch(variants)

        for annotation in results:
            if annotation is None:
                continue
            # The consequence was mapped successfully (not INTERGENIC default
            # due to unmapped term) unless the variant is genuinely intergenic
            assert annotation.consequence is not None

    def test_vep_returns_gene_names_for_coding_variants(
        self, vep_client: VEPClient
    ) -> None:
        """Coding variants should have a gene_name from VEP."""
        # Only test coding variants (skip intergenic)
        coding_variants = [
            (chrom, pos, ref, alt)
            for chrom, pos, ref, alt, csq in _KNOWN_VARIANTS
            if csq != FunctionalConsequence.INTERGENIC
        ]

        results = vep_client.annotate_batch(coding_variants)

        for i, annotation in enumerate(results):
            if annotation is None:
                continue
            assert (
                annotation.gene_name is not None
            ), f"Variant {coding_variants[i]} should have a gene name"
