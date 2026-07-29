"""Unit tests for structural variant VCF parser."""

from __future__ import annotations

import pytest

from vartriage.structural.models import SVType, StructuralVariant
from vartriage.structural.parser import (
    SVParser,
    _BND_PATTERN,
    _SVTYPE_MAP,
    _SYMBOLIC_ALT_PATTERN,
)


class TestSVTypeMapping:
    """Tests for SV type resolution from INFO and ALT fields."""

    def test_standard_types_mapped(self) -> None:
        assert _SVTYPE_MAP["DEL"] == SVType.DEL
        assert _SVTYPE_MAP["DUP"] == SVType.DUP
        assert _SVTYPE_MAP["INV"] == SVType.INV
        assert _SVTYPE_MAP["INS"] == SVType.INS
        assert _SVTYPE_MAP["BND"] == SVType.BND
        assert _SVTYPE_MAP["CNV"] == SVType.CNV

    def test_tra_maps_to_bnd(self) -> None:
        assert _SVTYPE_MAP["TRA"] == SVType.BND

    def test_dup_tandem_maps_to_dup(self) -> None:
        assert _SVTYPE_MAP["DUP:TANDEM"] == SVType.DUP
        assert _SVTYPE_MAP["DUP:DISPERSED"] == SVType.DUP

    def test_mobile_element_deletions_map_to_del(self) -> None:
        assert _SVTYPE_MAP["DEL:ME"] == SVType.DEL
        assert _SVTYPE_MAP["DEL:ME:ALU"] == SVType.DEL

    def test_mobile_element_insertions_map_to_ins(self) -> None:
        assert _SVTYPE_MAP["INS:ME"] == SVType.INS
        assert _SVTYPE_MAP["INS:ME:ALU"] == SVType.INS


class TestSymbolicAltPattern:
    """Tests for symbolic ALT allele regex."""

    def test_matches_del(self) -> None:
        match = _SYMBOLIC_ALT_PATTERN.match("<DEL>")
        assert match is not None
        assert match.group(1) == "DEL"

    def test_matches_dup(self) -> None:
        match = _SYMBOLIC_ALT_PATTERN.match("<DUP>")
        assert match is not None
        assert match.group(1) == "DUP"

    def test_matches_inv(self) -> None:
        match = _SYMBOLIC_ALT_PATTERN.match("<INV>")
        assert match is not None

    def test_matches_ins(self) -> None:
        match = _SYMBOLIC_ALT_PATTERN.match("<INS>")
        assert match is not None

    def test_matches_cnv(self) -> None:
        match = _SYMBOLIC_ALT_PATTERN.match("<CNV>")
        assert match is not None

    def test_no_match_for_regular_allele(self) -> None:
        assert _SYMBOLIC_ALT_PATTERN.match("ATCG") is None

    def test_no_match_for_empty(self) -> None:
        assert _SYMBOLIC_ALT_PATTERN.match("") is None


class TestBNDPattern:
    """Tests for BND bracket notation regex."""

    def test_forward_bracket_notation(self) -> None:
        # t[p[ format: N[chr2:12345[
        match = _BND_PATTERN.match("N[chr2:12345[")
        assert match is not None

    def test_reverse_bracket_notation(self) -> None:
        # ]p]t format: ]chr2:12345]N
        match = _BND_PATTERN.match("]chr2:12345]N")
        assert match is not None

    def test_alt_forward_bracket(self) -> None:
        match = _BND_PATTERN.match("A]chr5:98765]")
        assert match is not None

    def test_no_match_for_symbolic(self) -> None:
        assert _BND_PATTERN.match("<DEL>") is None

    def test_no_match_for_regular(self) -> None:
        assert _BND_PATTERN.match("ATCG") is None


class TestStructuralVariantModel:
    """Tests for the StructuralVariant dataclass."""

    def test_length_from_svlen(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=1000, end=2000, sv_type=SVType.DEL,
            svlen=-1001,
        )
        assert sv.length == 1001

    def test_length_from_coordinates_when_svlen_none(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=1000, end=5000, sv_type=SVType.DUP,
        )
        assert sv.length == 4001

    def test_is_intrachromosomal_del(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=500, sv_type=SVType.DEL,
        )
        assert sv.is_intrachromosomal is True

    def test_bnd_not_intrachromosomal(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=100, sv_type=SVType.BND,
        )
        assert sv.is_intrachromosomal is False

    def test_frozen_dataclass(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=500, sv_type=SVType.DEL,
        )
        with pytest.raises(AttributeError):
            sv.start = 200  # type: ignore[misc]


class TestSVParserFiltering:
    """Tests for parser size and quality filtering logic."""

    def test_sv_below_min_size_filtered(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=120, sv_type=SVType.DEL,
            svlen=-21,
        )
        # Simulate the filter check
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 0
        parser._min_quality = 0.0
        parser._sv_types = None
        assert parser._passes_filters(sv) is False

    def test_sv_above_min_size_passes(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=200, sv_type=SVType.DEL,
            svlen=-101,
        )
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 0
        parser._min_quality = 0.0
        parser._sv_types = None
        assert parser._passes_filters(sv) is True

    def test_sv_above_max_size_filtered(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=200000, sv_type=SVType.DEL,
        )
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 100000
        parser._min_quality = 0.0
        parser._sv_types = None
        assert parser._passes_filters(sv) is False

    def test_bnd_bypasses_size_filter(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=100, sv_type=SVType.BND,
        )
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 0
        parser._min_quality = 0.0
        parser._sv_types = None
        assert parser._passes_filters(sv) is True

    def test_low_quality_filtered(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=500, sv_type=SVType.DEL,
            qual=10.0,
        )
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 0
        parser._min_quality = 20.0
        parser._sv_types = None
        assert parser._passes_filters(sv) is False

    def test_type_filter_excludes_unwanted(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=500, sv_type=SVType.INV,
        )
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 0
        parser._min_quality = 0.0
        parser._sv_types = {SVType.DEL, SVType.DUP}
        assert parser._passes_filters(sv) is False

    def test_type_filter_passes_wanted(self) -> None:
        sv = StructuralVariant(
            chrom="chr1", start=100, end=500, sv_type=SVType.DEL,
        )
        parser = object.__new__(SVParser)
        parser._min_size = 50
        parser._max_size = 0
        parser._min_quality = 0.0
        parser._sv_types = {SVType.DEL, SVType.DUP}
        assert parser._passes_filters(sv) is True


class TestSVParserFileValidation:
    """Tests for file validation in SVParser."""

    def test_missing_file_raises(self, tmp_path: "pytest.TempPathFactory") -> None:
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            SVParser(tmp_path / "nonexistent.vcf")  # type: ignore[arg-type]
