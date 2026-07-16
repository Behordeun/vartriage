"""Unit tests for VCF-to-VEP notation conversion."""

from __future__ import annotations

import pytest

from vartriage.api._notation import _strip_chr_prefix, vcf_to_vep_notation


class TestSNVs:
    """Single nucleotide variants."""

    def test_simple_snv(self) -> None:
        result = vcf_to_vep_notation("chr22", 17818804, "G", "A")
        assert result == "22 17818804 17818804 G/A +"

    def test_snv_without_chr_prefix(self) -> None:
        result = vcf_to_vep_notation("1", 12345, "C", "T")
        assert result == "1 12345 12345 C/T +"

    def test_snv_on_x_chromosome(self) -> None:
        result = vcf_to_vep_notation("chrX", 100000, "A", "G")
        assert result == "X 100000 100000 A/G +"


class TestDeletions:
    """Variants where ref is longer than alt (VCF padding base)."""

    def test_single_base_deletion(self) -> None:
        # VCF: pos=100, ref=AT, alt=A (deletes the T at pos 101)
        result = vcf_to_vep_notation("chr1", 100, "AT", "A")
        assert result == "1 101 101 T/- +"

    def test_three_base_deletion(self) -> None:
        # VCF: pos=100, ref=AGCT, alt=A (deletes GCT at pos 101-103)
        result = vcf_to_vep_notation("chr1", 100, "AGCT", "A")
        assert result == "1 101 103 GCT/- +"

    def test_deletion_preserves_end_coordinate(self) -> None:
        # VCF: pos=500, ref=GC, alt=G (single base del at 501)
        result = vcf_to_vep_notation("22", 500, "GC", "G")
        assert result == "22 501 501 C/- +"


class TestInsertions:
    """Variants where alt is longer than ref (VCF padding base)."""

    def test_single_base_insertion(self) -> None:
        # VCF: pos=100, ref=A, alt=AC (inserts C after pos 100)
        result = vcf_to_vep_notation("chr1", 100, "A", "AC")
        assert result == "1 100 101 -/C +"

    def test_three_base_insertion(self) -> None:
        # VCF: pos=200, ref=G, alt=GACT (inserts ACT)
        result = vcf_to_vep_notation("chr5", 200, "G", "GACT")
        assert result == "5 200 201 -/ACT +"


class TestMNVs:
    """Multi-nucleotide variants (same length ref and alt > 1bp)."""

    def test_dinucleotide_substitution(self) -> None:
        result = vcf_to_vep_notation("chr7", 1000, "AG", "CT")
        assert result == "7 1000 1001 AG/CT +"

    def test_trinucleotide_substitution(self) -> None:
        result = vcf_to_vep_notation("chr3", 500, "TAC", "GGG")
        assert result == "3 500 502 TAC/GGG +"


class TestComplexIndels:
    """Both ref and alt > 1bp with different lengths."""

    def test_complex_del_ins(self) -> None:
        # VCF: pos=100, ref=ACTG, alt=AG (padding A, delete CTG, insert G)
        # After stripping padding: delete CTG, insert G at pos 101
        result = vcf_to_vep_notation("chr2", 100, "ACTG", "AG")
        assert result == "2 101 103 CTG/G +"

    def test_complex_shorter_del_longer_ins(self) -> None:
        # VCF: pos=50, ref=AC, alt=ATTTG (padding A, del C, ins TTTG)
        result = vcf_to_vep_notation("chr11", 50, "AC", "ATTTG")
        assert result == "11 51 51 C/TTTG +"


class TestChrPrefixStripping:
    """Chromosome prefix handling."""

    def test_strips_chr_from_autosome(self) -> None:
        assert _strip_chr_prefix("chr1") == "1"
        assert _strip_chr_prefix("chr22") == "22"

    def test_strips_chr_from_sex_chromosomes(self) -> None:
        assert _strip_chr_prefix("chrX") == "X"
        assert _strip_chr_prefix("chrY") == "Y"

    def test_converts_chrM_to_MT(self) -> None:
        assert _strip_chr_prefix("chrM") == "MT"

    def test_no_prefix_passes_through(self) -> None:
        assert _strip_chr_prefix("22") == "22"
        assert _strip_chr_prefix("X") == "X"
        assert _strip_chr_prefix("MT") == "MT"

    def test_handles_lowercase_chr(self) -> None:
        assert _strip_chr_prefix("Chr1") == "1"
        assert _strip_chr_prefix("CHR22") == "22"
