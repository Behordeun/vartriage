"""Unit tests for the VCF report writer."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import pysam
import pytest

from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.reporting.vcf_writer import (
    _build_lookup,
    _inject_info_fields,
    write_vcf,
)


def _make_classified(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    allele_frequency: Optional[float] = 0.001,
    composite_rank: Optional[float] = 0.85,
    classification: ACMGClassification = ACMGClassification.LIKELY_PATHOGENIC,
    evidence_tags: frozenset[EvidenceTag] = frozenset({EvidenceTag.PP3}),
) -> ClassifiedVariant:
    """Build a ClassifiedVariant for testing."""
    variant = Variant(
        chrom=chrom, pos=pos, id=None,
        ref=ref, alt=alt, qual=30.0, filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=allele_frequency,
    )
    scored = ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
    )
    return ClassifiedVariant(
        scored=scored,
        classification=classification,
        evidence_tags=evidence_tags,
    )


def _create_source_vcf(path: Path, records: list[dict]) -> None:
    """Write a minimal source VCF for testing."""
    header = pysam.VariantHeader()
    header.add_sample("SAMPLE")

    chroms = set()
    for rec in records:
        chroms.add(rec["chrom"])

    for chrom in sorted(chroms):
        header.add_line(
            f"##contig=<ID={chrom},length=1000000>"
        )

    header.add_line(
        '##FORMAT=<ID=GT,Number=1,Type=String,'
        'Description="Genotype">'
    )

    with pysam.VariantFile(str(path), "wz", header=header) as out:
        for rec in records:
            alleles = (rec["ref"],) + tuple(rec.get("alts", [rec.get("alt")]))
            new_rec = out.new_record(
                contig=rec["chrom"],
                start=rec["pos"] - 1,
                alleles=alleles,
            )
            new_rec.samples["SAMPLE"]["GT"] = (0, 1)
            out.write(new_rec)

    pysam.tabix_index(str(path), preset="vcf", force=True)


class TestBuildLookup:
    """Tests for _build_lookup."""

    def test_empty_list(self) -> None:
        assert _build_lookup([]) == {}

    def test_single_variant(self) -> None:
        cv = _make_classified(chrom="chr1", pos=100, ref="A", alt="T")
        lookup = _build_lookup([cv])
        assert ("chr1", 100, "A", "T") in lookup
        assert lookup[("chr1", 100, "A", "T")] is cv

    def test_duplicate_keys_last_wins(self) -> None:
        first = _make_classified(
            chrom="chr2", pos=300, ref="C", alt="G",
            classification=ACMGClassification.VUS,
        )
        second = _make_classified(
            chrom="chr2", pos=300, ref="C", alt="G",
            classification=ACMGClassification.PATHOGENIC,
        )
        lookup = _build_lookup([first, second])
        assert lookup[("chr2", 300, "C", "G")] is second


class TestInjectInfoFields:
    """Tests for _inject_info_fields."""

    def test_mandatory_fields_always_present(self, tmp_path: Path) -> None:
        cv = _make_classified(
            allele_frequency=None,
            composite_rank=None,
            evidence_tags=frozenset(),
        )
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        with pysam.VariantFile(str(output)) as vcf:
            for rec in vcf:
                assert "VARTRIAGE_CONSEQUENCE" in rec.info
                assert "VARTRIAGE_ACMG" in rec.info
                assert "VARTRIAGE_AF" not in rec.info
                assert "VARTRIAGE_RANK" not in rec.info
                assert "VARTRIAGE_TAGS" not in rec.info

    def test_optional_fields_present_when_non_null(
        self, tmp_path: Path
    ) -> None:
        cv = _make_classified(
            allele_frequency=0.005,
            composite_rank=0.9,
            evidence_tags=frozenset({EvidenceTag.PVS1, EvidenceTag.PM2}),
        )
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        with pysam.VariantFile(str(output)) as vcf:
            for rec in vcf:
                assert rec.info["VARTRIAGE_AF"] == pytest.approx(0.005)
                assert rec.info["VARTRIAGE_RANK"] == pytest.approx(0.9)
                tags = rec.info["VARTRIAGE_TAGS"]
                assert "PM2" in tags
                assert "PVS1" in tags


class TestWriteVcfUnmatched:
    """Unmatched records pass through without VARTRIAGE_* fields."""

    def test_unmatched_variant_has_no_vartriage_fields(
        self, tmp_path: Path
    ) -> None:
        cv = _make_classified(chrom="chr1", pos=100, ref="A", alt="T")
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
            {"chrom": "chr1", "pos": 200, "ref": "G", "alt": "C"},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        with pysam.VariantFile(str(output)) as vcf:
            records = list(vcf)
            assert len(records) == 2

            matched = records[0]
            assert "VARTRIAGE_ACMG" in matched.info

            unmatched = records[1]
            for key in unmatched.info:
                assert not key.startswith("VARTRIAGE_")


class TestWriteVcfMultiAllelic:
    """Multi-allelic records: only first ALT is matched."""

    def test_first_alt_annotated_second_ignored(
        self, tmp_path: Path
    ) -> None:
        cv = _make_classified(chrom="chr1", pos=100, ref="A", alt="T")
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alts": ["T", "C"]},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        with pysam.VariantFile(str(output)) as vcf:
            for rec in vcf:
                assert "VARTRIAGE_ACMG" in rec.info


class TestWriteVcfAtomicity:
    """Atomic write: both VCF and .tbi produced, or neither."""

    def test_output_and_index_both_exist(self, tmp_path: Path) -> None:
        cv = _make_classified()
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        assert output.exists()
        assert Path(str(output) + ".tbi").exists()

    def test_no_temp_files_remain_on_success(
        self, tmp_path: Path
    ) -> None:
        cv = _make_classified()
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        tmp_files = list(tmp_path.glob("*.tmp*"))
        assert tmp_files == []


class TestWriteVcfRecordCount:
    """All source records are included in output regardless of matching."""

    def test_all_records_preserved(self, tmp_path: Path) -> None:
        cv = _make_classified(chrom="chr1", pos=100, ref="A", alt="T")
        source = tmp_path / "src.vcf.gz"
        _create_source_vcf(source, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
            {"chrom": "chr1", "pos": 200, "ref": "G", "alt": "C"},
            {"chrom": "chr1", "pos": 300, "ref": "T", "alt": "A"},
        ])

        output = tmp_path / "out.vcf.gz"
        write_vcf([cv], source, output)

        with pysam.VariantFile(str(output)) as vcf:
            records = list(vcf)
            assert len(records) == 3
