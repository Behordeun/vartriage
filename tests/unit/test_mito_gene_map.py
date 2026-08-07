"""Unit tests for the mitochondrial gene map."""

from __future__ import annotations

import pytest

from vartriage.mito.gene_map import MtGeneMap


@pytest.fixture
def gene_map() -> MtGeneMap:
    """Load the bundled gene map."""
    return MtGeneMap()


class TestMtGeneMapLoading:
    """Verify the gene map loads correctly from bundled data."""

    def test_loads_at_least_37_entries(self, gene_map: MtGeneMap):
        # 37 genes + 2 D-loop entries = 39 minimum
        assert len(gene_map.entries) >= 39

    def test_first_entry_is_dloop(self, gene_map: MtGeneMap):
        first = gene_map.entries[0]
        assert first.gene_name == "MT-Dloop"
        assert first.start == 1


class TestMtGeneMapQuery:
    """Verify position-based gene lookups."""

    def test_position_in_nd1_returns_protein_coding(self, gene_map: MtGeneMap):
        # MT-ND1: 3307-4262
        ctx = gene_map.query(3500)
        assert ctx.gene_name == "MT-ND1"
        assert ctx.gene_type == "protein_coding"
        assert ctx.is_in_coding_or_trna is True

    def test_position_in_tl1_returns_trna(self, gene_map: MtGeneMap):
        # MT-TL1: 3230-3304
        ctx = gene_map.query(3243)
        assert ctx.gene_name == "MT-TL1"
        assert ctx.gene_type == "tRNA"
        assert ctx.is_in_coding_or_trna is True

    def test_position_in_rnr1_returns_rrna(self, gene_map: MtGeneMap):
        # MT-RNR1: 648-1601
        ctx = gene_map.query(1000)
        assert ctx.gene_name == "MT-RNR1"
        assert ctx.gene_type == "rRNA"
        assert ctx.is_in_coding_or_trna is False

    def test_position_in_dloop_returns_control_region(self, gene_map: MtGeneMap):
        # MT-Dloop: 16024-16569
        ctx = gene_map.query(16100)
        assert ctx.gene_name == "MT-Dloop"
        assert ctx.gene_type == "control_region"
        assert ctx.is_in_coding_or_trna is False

    def test_intergenic_position(self, gene_map: MtGeneMap):
        # Between MT-TS1 (7446-7514) and MT-TD (7518-7585)
        # Positions 7515-7517 are intergenic
        ctx = gene_map.query(7516)
        assert ctx.gene_type == "intergenic"
        assert ctx.gene_name is None
        assert ctx.is_in_coding_or_trna is False

    def test_position_at_gene_boundary_start(self, gene_map: MtGeneMap):
        # MT-ND1 starts at 3307
        ctx = gene_map.query(3307)
        assert ctx.gene_name == "MT-ND1"

    def test_position_at_gene_boundary_end(self, gene_map: MtGeneMap):
        # MT-ND1 ends at 4262
        ctx = gene_map.query(4262)
        assert ctx.gene_name == "MT-ND1"
