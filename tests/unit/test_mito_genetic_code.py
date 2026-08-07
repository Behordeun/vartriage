"""Unit tests for the mitochondrial genetic code module."""

from __future__ import annotations

import pytest

from vartriage._internal.genetic_code import translate_codon
from vartriage.mito.genetic_code import (
    MT_CODON_OVERRIDES,
    is_mitochondrial,
    translate_codon_mt,
)


class TestMtCodonOverrides:
    """Verify the 4 mitochondrial codon differences from the standard table."""

    def test_tga_is_trp_in_mito(self):
        assert translate_codon_mt("TGA") == "W"

    def test_tga_is_stop_in_standard(self):
        assert translate_codon("TGA") == "*"

    def test_ata_is_met_in_mito(self):
        assert translate_codon_mt("ATA") == "M"

    def test_ata_is_ile_in_standard(self):
        assert translate_codon("ATA") == "I"

    def test_aga_is_stop_in_mito(self):
        assert translate_codon_mt("AGA") == "*"

    def test_aga_is_arg_in_standard(self):
        assert translate_codon("AGA") == "R"

    def test_agg_is_stop_in_mito(self):
        assert translate_codon_mt("AGG") == "*"

    def test_agg_is_arg_in_standard(self):
        assert translate_codon("AGG") == "R"

    def test_override_dict_has_exactly_four_entries(self):
        assert len(MT_CODON_OVERRIDES) == 4


class TestTranslateCodonMtFallthrough:
    """Non-overridden codons should match standard translation."""

    def test_atg_is_met(self):
        assert translate_codon_mt("ATG") == "M"

    def test_ttt_is_phe(self):
        assert translate_codon_mt("TTT") == "F"

    def test_taa_is_stop(self):
        assert translate_codon_mt("TAA") == "*"

    def test_ggg_is_gly(self):
        assert translate_codon_mt("GGG") == "G"

    def test_lowercase_input_normalized(self):
        assert translate_codon_mt("tga") == "W"

    def test_ambiguous_codon_returns_question_mark(self):
        assert translate_codon_mt("NNN") == "?"

    def test_invalid_codon_returns_question_mark(self):
        assert translate_codon_mt("XYZ") == "?"


class TestIsMitochondrial:
    """Verify chromosome name detection for mtDNA."""

    @pytest.mark.parametrize("chrom", ["chrM", "MT", "M", "CHRM", "mt", "chrm"])
    def test_mitochondrial_names_detected(self, chrom: str):
        assert is_mitochondrial(chrom) is True

    @pytest.mark.parametrize(
        "chrom", ["chr1", "chr22", "chrX", "chrY", "1", "X", "chr10"]
    )
    def test_nuclear_chromosomes_not_mitochondrial(self, chrom: str):
        assert is_mitochondrial(chrom) is False
