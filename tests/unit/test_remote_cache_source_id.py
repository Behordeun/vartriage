"""Remote cache source ids are qualified by build, version, and dataset.

A cached remote score is keyed on its source id plus coordinates. If the
source id does not encode the reference build and dataset, a GRCh37 run and
a GRCh38 run (or exomes vs genomes) share cache rows and read each other's
values. The source-id resolver derives a distinct id per preset so that
cannot happen, and gives a raw URL its own stable id.
"""

from __future__ import annotations

from vartriage.remote.presets import resolve_cache_source_id


class TestCacheSourceIdIsolation:
    def test_different_builds_get_different_ids(self) -> None:
        grch38 = resolve_cache_source_id("cadd", "cadd-v1.7-grch38")
        grch37 = resolve_cache_source_id("cadd", "cadd-v1.7-grch37")
        assert grch38 != grch37

    def test_exomes_and_genomes_get_different_ids(self) -> None:
        exomes = resolve_cache_source_id("gnomad", "gnomad-exomes-v4-grch38")
        genomes = resolve_cache_source_id("gnomad", "gnomad-genomes-v4-grch38")
        assert exomes != genomes

    def test_snv_and_indel_cadd_get_different_ids(self) -> None:
        snv = resolve_cache_source_id("cadd", "cadd-v1.7-grch38")
        indel = resolve_cache_source_id("cadd", "cadd-v1.7-indels-grch38")
        assert snv != indel

    def test_id_is_stable_for_same_preset(self) -> None:
        a = resolve_cache_source_id("gnomad", "gnomad-genomes-v4-grch38")
        b = resolve_cache_source_id("gnomad", "gnomad-genomes-v4-grch38")
        assert a == b

    def test_id_carries_the_base_source(self) -> None:
        # The id remains attributable to its scoring system for cache stats.
        assert resolve_cache_source_id("cadd", "cadd-v1.7-grch38").startswith("cadd")
        assert resolve_cache_source_id("gnomad", "gnomad-genomes-v4-grch38").startswith(
            "gnomad"
        )

    def test_raw_url_gets_a_stable_distinct_id(self) -> None:
        url_a = "https://example.org/a/GRCh38/scores.tsv.gz"
        url_b = "https://example.org/b/GRCh37/scores.tsv.gz"
        id_a = resolve_cache_source_id("cadd", url_a)
        id_b = resolve_cache_source_id("cadd", url_b)
        assert id_a != id_b
        assert id_a == resolve_cache_source_id("cadd", url_a)
