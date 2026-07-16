"""Unit tests for VEP SO term to FunctionalConsequence mapping."""

from __future__ import annotations

import pytest

from vartriage.api._consequence_map import (map_so_term, map_vep_most_severe,
                                            most_severe_consequence)
from vartriage.models.variant import FunctionalConsequence


class TestMapSoTerm:
    """Single-term mapping."""

    def test_missense_variant(self) -> None:
        consequence, rank = map_so_term("missense_variant")
        assert consequence == FunctionalConsequence.MISSENSE
        assert rank == 7

    def test_frameshift_variant(self) -> None:
        consequence, rank = map_so_term("frameshift_variant")
        assert consequence == FunctionalConsequence.FRAMESHIFT
        assert rank == 3

    def test_stop_gained(self) -> None:
        consequence, rank = map_so_term("stop_gained")
        assert consequence == FunctionalConsequence.NONSENSE
        assert rank == 2

    def test_splice_donor_variant(self) -> None:
        consequence, rank = map_so_term("splice_donor_variant")
        assert consequence == FunctionalConsequence.SPLICE_SITE
        assert rank == 1

    def test_splice_acceptor_variant(self) -> None:
        consequence, rank = map_so_term("splice_acceptor_variant")
        assert consequence == FunctionalConsequence.SPLICE_SITE
        assert rank == 1

    def test_splice_region_variant(self) -> None:
        consequence, rank = map_so_term("splice_region_variant")
        assert consequence == FunctionalConsequence.SPLICE_SITE
        assert rank == 8

    def test_synonymous_variant(self) -> None:
        consequence, rank = map_so_term("synonymous_variant")
        assert consequence == FunctionalConsequence.SYNONYMOUS
        assert rank == 10

    def test_intron_variant(self) -> None:
        consequence, rank = map_so_term("intron_variant")
        assert consequence == FunctionalConsequence.SYNONYMOUS
        assert rank == 13

    def test_intergenic_variant(self) -> None:
        consequence, rank = map_so_term("intergenic_variant")
        assert consequence == FunctionalConsequence.INTERGENIC
        assert rank == 14

    def test_upstream_gene_variant(self) -> None:
        consequence, rank = map_so_term("upstream_gene_variant")
        assert consequence == FunctionalConsequence.INTERGENIC
        assert rank == 14

    def test_inframe_insertion(self) -> None:
        consequence, rank = map_so_term("inframe_insertion")
        assert consequence == FunctionalConsequence.IN_FRAME_INSERTION
        assert rank == 5

    def test_inframe_deletion(self) -> None:
        consequence, rank = map_so_term("inframe_deletion")
        assert consequence == FunctionalConsequence.IN_FRAME_DELETION
        assert rank == 6

    def test_stop_lost(self) -> None:
        consequence, rank = map_so_term("stop_lost")
        assert consequence == FunctionalConsequence.MISSENSE
        assert rank == 4

    def test_utr_variants(self) -> None:
        consequence, _ = map_so_term("5_prime_UTR_variant")
        assert consequence == FunctionalConsequence.SYNONYMOUS
        consequence, _ = map_so_term("3_prime_UTR_variant")
        assert consequence == FunctionalConsequence.SYNONYMOUS

    def test_unknown_term_defaults_to_intergenic(self) -> None:
        consequence, rank = map_so_term("totally_unknown_variant_type")
        assert consequence == FunctionalConsequence.INTERGENIC
        assert rank == 99


class TestMostSevereConsequence:
    """Selecting most severe from multiple terms."""

    def test_empty_list_returns_intergenic(self) -> None:
        assert most_severe_consequence([]) == FunctionalConsequence.INTERGENIC

    def test_single_term(self) -> None:
        result = most_severe_consequence(["missense_variant"])
        assert result == FunctionalConsequence.MISSENSE

    def test_selects_most_severe_from_multiple(self) -> None:
        terms = ["synonymous_variant", "missense_variant", "intron_variant"]
        result = most_severe_consequence(terms)
        assert result == FunctionalConsequence.MISSENSE

    def test_frameshift_beats_missense(self) -> None:
        terms = ["missense_variant", "frameshift_variant"]
        result = most_severe_consequence(terms)
        assert result == FunctionalConsequence.FRAMESHIFT

    def test_splice_donor_beats_frameshift(self) -> None:
        terms = ["frameshift_variant", "splice_donor_variant", "intron_variant"]
        result = most_severe_consequence(terms)
        assert result == FunctionalConsequence.SPLICE_SITE

    def test_stop_gained_beats_splice_region(self) -> None:
        terms = ["splice_region_variant", "stop_gained"]
        result = most_severe_consequence(terms)
        assert result == FunctionalConsequence.NONSENSE

    def test_unknown_term_among_known_uses_known(self) -> None:
        terms = ["unknown_thing", "synonymous_variant"]
        result = most_severe_consequence(terms)
        assert result == FunctionalConsequence.SYNONYMOUS


class TestMapVepMostSevere:
    """Mapping VEP's pre-computed most_severe_consequence field."""

    def test_maps_missense(self) -> None:
        assert map_vep_most_severe("missense_variant") == FunctionalConsequence.MISSENSE

    def test_maps_frameshift(self) -> None:
        assert (
            map_vep_most_severe("frameshift_variant")
            == FunctionalConsequence.FRAMESHIFT
        )

    def test_none_returns_intergenic(self) -> None:
        assert map_vep_most_severe(None) == FunctionalConsequence.INTERGENIC

    def test_maps_splice_acceptor(self) -> None:
        assert (
            map_vep_most_severe("splice_acceptor_variant")
            == FunctionalConsequence.SPLICE_SITE
        )


class TestSeverityOrdering:
    """Verify the severity ranking is internally consistent."""

    def test_splice_more_severe_than_nonsense(self) -> None:
        _, splice_rank = map_so_term("splice_donor_variant")
        _, nonsense_rank = map_so_term("stop_gained")
        assert splice_rank < nonsense_rank

    def test_nonsense_more_severe_than_frameshift(self) -> None:
        _, nonsense_rank = map_so_term("stop_gained")
        _, frame_rank = map_so_term("frameshift_variant")
        assert nonsense_rank < frame_rank

    def test_frameshift_more_severe_than_missense(self) -> None:
        _, frame_rank = map_so_term("frameshift_variant")
        _, missense_rank = map_so_term("missense_variant")
        assert frame_rank < missense_rank

    def test_missense_more_severe_than_synonymous(self) -> None:
        _, missense_rank = map_so_term("missense_variant")
        _, syn_rank = map_so_term("synonymous_variant")
        assert missense_rank < syn_rank

    def test_synonymous_more_severe_than_intergenic(self) -> None:
        _, syn_rank = map_so_term("synonymous_variant")
        _, inter_rank = map_so_term("intergenic_variant")
        assert syn_rank < inter_rank
