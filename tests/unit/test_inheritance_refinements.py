"""Tests for v0.19.0 NMD-escape geometry and inheritance-filter refinements."""

from __future__ import annotations

from vartriage.annotation.transcript_index import TranscriptCDSIndex
from vartriage.filter.inheritance_filter import InheritanceFilter
from vartriage.models.config import InheritanceConfig
from vartriage.models.variant import Variant


def _index(exons: list[tuple[int, int]], strand: str) -> TranscriptCDSIndex:
    index = TranscriptCDSIndex()
    for start, end in exons:
        index.add_cds_exon("t1", "GENE", "chr1", start, end, strand, 0)
    index.finalize()
    return index


class TestNMDEscapeZone:
    def test_plus_strand_last_exon_and_junction(self) -> None:
        zone = _index([(0, 100), (200, 300)], "+").escape_zone("GENE", "chr1", 250)
        assert zone is not None
        assert not zone.is_single_exon
        assert zone.last_exon_start == 200
        assert zone.last_exon_end == 300
        assert zone.last_junction_pos == 200
        assert zone.is_in_last_exon(250)
        assert not zone.is_in_last_exon(50)

    def test_minus_strand_last_exon_is_lowest_coordinate(self) -> None:
        zone = _index([(0, 100), (200, 300)], "-").escape_zone("GENE", "chr1", 50)
        assert zone is not None
        assert zone.last_exon_start == 0
        assert zone.last_exon_end == 100
        assert zone.last_junction_pos == 99

    def test_single_exon_flag(self) -> None:
        zone = _index([(0, 300)], "+").escape_zone("GENE", "chr1", 150)
        assert zone is not None
        assert zone.is_single_exon

    def test_near_last_junction_margin(self) -> None:
        zone = _index([(0, 100), (200, 300)], "+").escape_zone("GENE", "chr1", 250)
        assert zone is not None
        assert zone.is_near_last_junction(230, margin=50)
        assert not zone.is_near_last_junction(260, margin=50)

    def test_wrong_gene_returns_none(self) -> None:
        assert _index([(0, 100)], "+").escape_zone("OTHER", "chr1", 50) is None

    def test_position_outside_cds_returns_none(self) -> None:
        assert _index([(0, 100)], "+").escape_zone("GENE", "chr1", 5000) is None


def _config() -> InheritanceConfig:
    return InheritanceConfig(
        proband="child",
        mother="mom",
        father="dad",
        patterns=["de_novo", "dominant", "recessive", "x_linked"],
    )


def _variant(
    proband_gt: tuple,
    mother_gt: tuple,
    father_gt: tuple,
    gene: str | None = None,
) -> Variant:
    info: dict = {
        "_pysam_samples": {
            "child": {"GT": proband_gt, "GQ": 99},
            "mom": {"GT": mother_gt, "GQ": 80},
            "dad": {"GT": father_gt, "GQ": 75},
        }
    }
    if gene is not None:
        info["gene"] = gene
    return Variant("chr1", 100, None, "A", "T", 30.0, "PASS", info)


class TestConsanguinity:
    def test_standard_both_parents_het(self) -> None:
        filt = InheritanceFilter(_config(), ["child", "mom", "dad"])
        assert filt._classify_recessive("1/1", "0/1", "0/1")

    def test_one_parent_het_other_unavailable(self) -> None:
        filt = InheritanceFilter(_config(), ["child", "mom", "dad"])
        assert filt._classify_recessive("1/1", "0/1", "./.")

    def test_one_parent_het_other_hom_alt(self) -> None:
        filt = InheritanceFilter(_config(), ["child", "mom", "dad"])
        assert filt._classify_recessive("1/1", "0/1", "1/1")

    def test_rejected_when_other_parent_confirmed_hom_ref(self) -> None:
        filt = InheritanceFilter(_config(), ["child", "mom", "dad"])
        assert not filt._classify_recessive("1/1", "0/1", "0/0")

    def test_rejected_when_neither_parent_het(self) -> None:
        filt = InheritanceFilter(_config(), ["child", "mom", "dad"])
        assert not filt._classify_recessive("1/1", "0/0", "0/0")


class TestGeneModeFiltering:
    def test_ad_gene_drops_recessive_call(self) -> None:
        filt = InheritanceFilter(
            _config(),
            ["child", "mom", "dad"],
            gene_inheritance_modes={"BRCA1": "AD"},
        )
        out = list(filt.apply(iter([_variant((1, 1), (0, 1), (0, 1), gene="BRCA1")])))
        assert out[0].info["inheritance_pattern"] == []

    def test_ar_gene_keeps_recessive_call(self) -> None:
        filt = InheritanceFilter(
            _config(),
            ["child", "mom", "dad"],
            gene_inheritance_modes={"CFTR": "AR"},
        )
        out = list(filt.apply(iter([_variant((1, 1), (0, 1), (0, 1), gene="CFTR")])))
        assert "recessive" in out[0].info["inheritance_pattern"]

    def test_ad_gene_keeps_dominant_call(self) -> None:
        filt = InheritanceFilter(
            _config(),
            ["child", "mom", "dad"],
            gene_inheritance_modes={"BRCA1": "AD"},
        )
        out = list(filt.apply(iter([_variant((0, 1), (0, 1), (0, 0), gene="BRCA1")])))
        assert "dominant" in out[0].info["inheritance_pattern"]

    def test_unknown_gene_leaves_patterns_untouched(self) -> None:
        filt = InheritanceFilter(
            _config(),
            ["child", "mom", "dad"],
            gene_inheritance_modes={"BRCA1": "AD"},
        )
        out = list(filt.apply(iter([_variant((1, 1), (0, 1), (0, 1), gene="XYZ")])))
        assert "recessive" in out[0].info["inheritance_pattern"]

    def test_no_mode_map_leaves_patterns_untouched(self) -> None:
        filt = InheritanceFilter(_config(), ["child", "mom", "dad"])
        out = list(filt.apply(iter([_variant((1, 1), (0, 1), (0, 1), gene="BRCA1")])))
        assert "recessive" in out[0].info["inheritance_pattern"]
