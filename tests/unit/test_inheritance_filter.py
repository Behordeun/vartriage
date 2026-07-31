"""Unit tests for trio-based Mendelian inheritance pattern classification."""

from __future__ import annotations

import pytest

from vartriage.filter.inheritance_filter import InheritanceFilter
from vartriage.models.config import InheritanceConfig
from vartriage.models.variant import Variant


def _config(patterns: list[str] | None = None) -> InheritanceConfig:
    return InheritanceConfig(
        proband="child",
        mother="mom",
        father="dad",
        patterns=patterns or ["de_novo", "dominant", "recessive", "x_linked"],
    )


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    proband_gt: tuple = (0, 1),
    mother_gt: tuple = (0, 0),
    father_gt: tuple = (0, 0),
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
    return Variant(
        chrom=chrom, pos=pos, id=None, ref="A", alt="T",
        qual=30.0, filter_status="PASS", info=info,
    )


class TestConstructor:
    def test_rejects_missing_proband(self) -> None:
        config = _config()
        with pytest.raises(ValueError, match="Proband"):
            InheritanceFilter(config, ["mom", "dad"])

    def test_rejects_missing_mother(self) -> None:
        config = _config()
        with pytest.raises(ValueError, match="Mother"):
            InheritanceFilter(config, ["child", "dad"])

    def test_rejects_missing_father(self) -> None:
        config = _config()
        with pytest.raises(ValueError, match="Father"):
            InheritanceFilter(config, ["child", "mom"])

    def test_accepts_valid_samples(self) -> None:
        config = _config()
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        assert filt._proband == "child"


class TestGenotypeHelpers:
    def test_format_gt_diploid(self) -> None:
        assert InheritanceFilter._format_gt((0, 1)) == "0/1"

    def test_format_gt_missing(self) -> None:
        assert InheritanceFilter._format_gt((None, None)) == "./."

    def test_has_alt_allele_het(self) -> None:
        assert InheritanceFilter._has_alt_allele("0/1") is True

    def test_has_alt_allele_hom_ref(self) -> None:
        assert InheritanceFilter._has_alt_allele("0/0") is False

    def test_has_alt_allele_hom_alt(self) -> None:
        assert InheritanceFilter._has_alt_allele("1/1") is True

    def test_has_alt_allele_missing(self) -> None:
        assert InheritanceFilter._has_alt_allele("./.") is False

    def test_is_het_true(self) -> None:
        assert InheritanceFilter._is_het("0/1") is True

    def test_is_het_hom_alt_false(self) -> None:
        assert InheritanceFilter._is_het("1/1") is False

    def test_is_het_hom_ref_false(self) -> None:
        assert InheritanceFilter._is_het("0/0") is False

    def test_is_het_missing_false(self) -> None:
        assert InheritanceFilter._is_het("./.") is False

    def test_is_hom_alt_true(self) -> None:
        assert InheritanceFilter._is_hom_alt("1/1") is True

    def test_is_hom_alt_het_false(self) -> None:
        assert InheritanceFilter._is_hom_alt("0/1") is False

    def test_is_hom_ref_true(self) -> None:
        assert InheritanceFilter._is_hom_ref("0/0") is True

    def test_is_hom_ref_het_false(self) -> None:
        assert InheritanceFilter._is_hom_ref("0/1") is False

    def test_phased_separator(self) -> None:
        assert InheritanceFilter._is_het("0|1") is True
        assert InheritanceFilter._is_hom_alt("1|1") is True


class TestDeNovo:
    def test_de_novo_detected(self) -> None:
        config = _config(["de_novo"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 0), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert len(results) == 1
        assert "de_novo" in results[0].info["inheritance_pattern"]

    def test_not_de_novo_when_mother_carries(self) -> None:
        config = _config(["de_novo"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert results[0].info["inheritance_pattern"] == []


class TestDominant:
    def test_dominant_from_mother(self) -> None:
        config = _config(["dominant"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert "dominant" in results[0].info["inheritance_pattern"]

    def test_dominant_from_father(self) -> None:
        config = _config(["dominant"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 0), father_gt=(0, 1))
        results = list(filt.apply(iter([v])))
        assert "dominant" in results[0].info["inheritance_pattern"]

    def test_not_dominant_both_parents_het(self) -> None:
        config = _config(["dominant"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 1))
        results = list(filt.apply(iter([v])))
        assert "dominant" not in results[0].info["inheritance_pattern"]


class TestRecessive:
    def test_recessive_detected(self) -> None:
        config = _config(["recessive"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(1, 1), mother_gt=(0, 1), father_gt=(0, 1))
        results = list(filt.apply(iter([v])))
        assert "recessive" in results[0].info["inheritance_pattern"]

    def test_not_recessive_proband_het(self) -> None:
        config = _config(["recessive"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 1))
        results = list(filt.apply(iter([v])))
        assert "recessive" not in results[0].info["inheritance_pattern"]


class TestXLinked:
    def test_x_linked_detected(self) -> None:
        config = _config(["x_linked"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(chrom="chrX", proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert "x_linked" in results[0].info["inheritance_pattern"]

    def test_not_x_linked_on_autosome(self) -> None:
        config = _config(["x_linked"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(chrom="chr1", proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert "x_linked" not in results[0].info["inheritance_pattern"]


class TestCompoundHet:
    def test_compound_het_trans_pair(self) -> None:
        config = _config(["compound_het"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        # Variant 1: alt from mother (mother het, father hom-ref)
        v1 = _make_variant(pos=100, proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0), gene="BRCA1")
        # Variant 2: alt from father (father het, mother hom-ref)
        v2 = _make_variant(pos=200, proband_gt=(0, 1), mother_gt=(0, 0), father_gt=(0, 1), gene="BRCA1")

        results = list(filt.apply(iter([v1, v2])))
        assert len(results) == 2
        assert "compound_het" in results[0].info["inheritance_pattern"]
        assert "compound_het" in results[1].info["inheritance_pattern"]

    def test_not_compound_het_same_parent(self) -> None:
        config = _config(["compound_het"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        # Both variants from mother — cis, not trans
        v1 = _make_variant(pos=100, proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0), gene="TP53")
        v2 = _make_variant(pos=200, proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0), gene="TP53")

        results = list(filt.apply(iter([v1, v2])))
        assert "compound_het" not in results[0].info["inheritance_pattern"]
        assert "compound_het" not in results[1].info["inheritance_pattern"]

    def test_gene_boundary_triggers_flush(self) -> None:
        config = _config(["compound_het"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v1 = _make_variant(pos=100, proband_gt=(0, 1), mother_gt=(0, 1), father_gt=(0, 0), gene="GENE_A")
        v2 = _make_variant(pos=200, proband_gt=(0, 1), mother_gt=(0, 0), father_gt=(0, 1), gene="GENE_B")

        results = list(filt.apply(iter([v1, v2])))
        assert len(results) == 2
        # Different genes, so no compound het
        assert "compound_het" not in results[0].info["inheritance_pattern"]
        assert "compound_het" not in results[1].info["inheritance_pattern"]


class TestOutputEnrichment:
    def test_output_strips_pysam_samples(self) -> None:
        config = _config(["de_novo"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 0), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert "_pysam_samples" not in results[0].info

    def test_output_contains_sample_info(self) -> None:
        config = _config(["de_novo"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 1), mother_gt=(0, 0), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert results[0].info["sample_gt"] == "0/1"
        assert results[0].info["sample_name"] == "child"
        assert results[0].info["sample_gq"] == 99

    def test_skips_variant_when_proband_has_no_alt(self) -> None:
        config = _config(["de_novo"])
        filt = InheritanceFilter(config, ["child", "mom", "dad"])
        v = _make_variant(proband_gt=(0, 0), mother_gt=(0, 0), father_gt=(0, 0))
        results = list(filt.apply(iter([v])))
        assert results == []
