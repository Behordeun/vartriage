"""Unit tests for ClinVar truth-set construction (validation-benchmarks FR-1)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from scripts.build_truth_set import build_truth_set

# A minimal GRCh38 FASTA with the contigs the fixtures touch. Written with a
# .fai alongside so pysam can open it without a separate index step.
_FASTA = ">chr1\n" + ("ACGT" * 50) + "\n>chr17\n" + ("ACGT" * 50) + "\n"


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    import pysam

    fa = tmp_path / "ref.fa"
    fa.write_text(_FASTA, encoding="utf-8")
    pysam.faidx(str(fa))
    return fa


def _clinvar_vcf(tmp_path: Path, records: list[str]) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="sig">\n'
        '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="rev">\n'
        '##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="conf">\n'
        '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="gene">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    path = tmp_path / "clinvar.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(header)
        for record in records:
            handle.write(record + "\n")
    return path


def _row(chrom, pos, ref, alt, sig, rev, gene="GENEX", conf=None):
    info = f"CLNSIG={sig};CLNREVSTAT={rev};GENEINFO={gene}:1"
    if conf:
        info += f";CLNSIGCONF={conf}"
    return f"{chrom}\t{pos}\t1\t{ref}\t{alt}\t.\t.\t{info}"


def _run(tmp_path: Path, fasta: Path, records: list[str]):
    clinvar = _clinvar_vcf(tmp_path, records)
    out = tmp_path / "truth"
    counters = build_truth_set(clinvar, out, fasta)
    labels = (out / "truth_labels.tsv").read_text(encoding="utf-8").splitlines()[1:]
    provenance = json.loads((out / "truth_set_provenance.json").read_text())
    return counters, labels, provenance


class TestTruthSetFiltering:
    def test_expert_panel_pathogenic_is_positive(self, tmp_path, fasta):
        counters, labels, _ = _run(
            tmp_path,
            fasta,
            [_row("1", 10, "A", "T", "Pathogenic", "reviewed_by_expert_panel")],
        )
        assert counters.positive == 1
        assert counters.negative == 0
        assert labels[0].endswith("positive")
        assert labels[0].startswith("chr1\t")

    def test_practice_guideline_benign_is_negative(self, tmp_path, fasta):
        counters, _, _ = _run(
            tmp_path,
            fasta,
            [_row("1", 20, "A", "T", "Benign", "practice_guideline")],
        )
        assert counters.negative == 1

    def test_single_submitter_is_rejected(self, tmp_path, fasta):
        counters, labels, _ = _run(
            tmp_path,
            fasta,
            [
                _row(
                    "1",
                    30,
                    "A",
                    "T",
                    "Pathogenic",
                    "criteria_provided,_single_submitter",
                )
            ],
        )
        assert counters.positive == 0
        assert counters.wrong_review_status == 1
        assert labels == []

    def test_conflicting_is_excluded(self, tmp_path, fasta):
        counters, _, _ = _run(
            tmp_path,
            fasta,
            [
                _row(
                    "1",
                    40,
                    "A",
                    "T",
                    "Pathogenic",
                    "reviewed_by_expert_panel",
                    conf="Benign(1)",
                )
            ],
        )
        assert counters.positive == 0
        assert counters.conflicting == 1

    def test_vus_significance_is_unmapped(self, tmp_path, fasta):
        counters, _, _ = _run(
            tmp_path,
            fasta,
            [
                _row(
                    "1",
                    50,
                    "A",
                    "T",
                    "Uncertain_significance",
                    "reviewed_by_expert_panel",
                )
            ],
        )
        assert counters.unmapped_significance == 1

    def test_compound_pathogenic_likely_pathogenic_is_positive(self, tmp_path, fasta):
        counters, _, _ = _run(
            tmp_path,
            fasta,
            [
                _row(
                    "1",
                    60,
                    "A",
                    "T",
                    "Pathogenic/Likely_pathogenic",
                    "reviewed_by_expert_panel",
                )
            ],
        )
        assert counters.positive == 1

    def test_multiallelic_is_dropped(self, tmp_path, fasta):
        counters, _, _ = _run(
            tmp_path,
            fasta,
            [_row("1", 70, "A", "T,C", "Pathogenic", "reviewed_by_expert_panel")],
        )
        assert counters.positive == 0
        assert counters.dropped_reasons.get("multiallelic_or_no_alt") == 1

    def test_bare_chrom_is_prefixed(self, tmp_path, fasta):
        _, labels, _ = _run(
            tmp_path,
            fasta,
            [_row("17", 80, "A", "T", "Pathogenic", "reviewed_by_expert_panel")],
        )
        assert labels[0].startswith("chr17\t")

    def test_provenance_pins_release_and_filters(self, tmp_path, fasta):
        _, _, provenance = _run(
            tmp_path,
            fasta,
            [_row("1", 90, "A", "T", "Pathogenic", "reviewed_by_expert_panel")],
        )
        assert provenance["clinvar_source"] == "clinvar.vcf.gz"
        assert provenance["reference"] == "GRCh38"
        assert provenance["filters"]["conflicting_excluded"] is True
        assert provenance["counts"]["truth_set_size"] == 1
