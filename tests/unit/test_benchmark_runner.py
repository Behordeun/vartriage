"""Unit tests for the benchmark runner (validation-benchmarks FR-2)."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.validation.benchmark import BenchmarkRunner, wilson_ci


def _classified(
    chrom: str, pos: int, ref: str, alt: str, tier: ACMGClassification
) -> ClassifiedVariant:
    annotated = AnnotatedVariant(
        variant=Variant(chrom, pos, None, ref, alt, 99.0, "PASS", {}),
        consequence=FunctionalConsequence.MISSENSE,
    )
    return ClassifiedVariant(
        scored=ScoredVariant(annotated=annotated), classification=tier
    )


def _labels_file(tmp_path: Path, rows: list[tuple]) -> Path:
    path = tmp_path / "truth_labels.tsv"
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("chrom\tpos\tref\talt\tgene\tclinvar_class\tlabel\n")
        for chrom, pos, ref, alt, label in rows:
            handle.write(f"{chrom}\t{pos}\t{ref}\t{alt}\tGENE\tsig\t{label}\n")
    return path


class TestWilsonCI:
    def test_known_value_half_successes(self) -> None:
        # 50/100 -> Wilson 95% CI is approximately [0.404, 0.596].
        point, low, high = wilson_ci(50, 100)
        assert point == pytest.approx(0.5)
        assert low == pytest.approx(0.404, abs=0.005)
        assert high == pytest.approx(0.596, abs=0.005)

    def test_perfect_success_upper_bound_is_one(self) -> None:
        point, low, high = wilson_ci(20, 20)
        assert point == pytest.approx(1.0)
        assert high == pytest.approx(1.0)
        assert low < 1.0

    def test_zero_denominator_is_zeros(self) -> None:
        assert wilson_ci(0, 0) == (0.0, 0.0, 0.0)

    def test_bounds_are_clamped_to_unit_interval(self) -> None:
        _, low, high = wilson_ci(1, 1000)
        assert low >= 0.0
        assert high <= 1.0


class TestBenchmarkTallying:
    def test_perfect_concordance(self, tmp_path: Path) -> None:
        labels = _labels_file(
            tmp_path,
            [
                ("chr1", 10, "A", "T", "positive"),
                ("chr1", 20, "C", "G", "negative"),
            ],
        )
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter(
                [
                    _classified("chr1", 10, "A", "T", ACMGClassification.PATHOGENIC),
                    _classified("chr1", 20, "C", "G", ACMGClassification.BENIGN),
                ]
            )
        )
        assert report.counts["tp"] == 1
        assert report.counts["tn"] == 1
        assert report.counts["fp"] == 0
        assert report.counts["fn"] == 0
        assert report.ppv.value == pytest.approx(1.0)
        assert report.sensitivity.value == pytest.approx(1.0)
        assert report.cohen_kappa == pytest.approx(1.0)

    def test_false_positive_lowers_ppv(self, tmp_path: Path) -> None:
        labels = _labels_file(
            tmp_path,
            [
                ("chr1", 10, "A", "T", "positive"),
                ("chr1", 20, "C", "G", "negative"),
            ],
        )
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter(
                [
                    _classified("chr1", 10, "A", "T", ACMGClassification.PATHOGENIC),
                    _classified(
                        "chr1", 20, "C", "G", ACMGClassification.LIKELY_PATHOGENIC
                    ),
                ]
            )
        )
        assert report.counts["tp"] == 1
        assert report.counts["fp"] == 1
        assert report.ppv.value == pytest.approx(0.5)

    def test_vus_on_positive_is_a_miss_not_a_benign_call(self, tmp_path: Path) -> None:
        labels = _labels_file(tmp_path, [("chr1", 10, "A", "T", "positive")])
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter([_classified("chr1", 10, "A", "T", ACMGClassification.VUS)])
        )
        assert report.counts["fn"] == 1
        assert report.counts["vus_on_positive"] == 1
        assert report.sensitivity.value == pytest.approx(0.0)

    def test_vus_on_negative_is_noncall_not_true_negative(self, tmp_path: Path) -> None:
        labels = _labels_file(tmp_path, [("chr1", 20, "C", "G", "negative")])
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter([_classified("chr1", 20, "C", "G", ACMGClassification.VUS)])
        )
        assert report.counts["tn"] == 0
        assert report.counts["fp"] == 0
        assert report.counts["vus_on_negative"] == 1

    def test_unmatched_truth_is_counted(self, tmp_path: Path) -> None:
        labels = _labels_file(
            tmp_path,
            [
                ("chr1", 10, "A", "T", "positive"),
                ("chr1", 99, "A", "T", "positive"),
            ],
        )
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter([_classified("chr1", 10, "A", "T", ACMGClassification.PATHOGENIC)])
        )
        assert report.matched == 1
        assert report.unmatched_truth == 1

    def test_prediction_not_in_truth_is_ignored(self, tmp_path: Path) -> None:
        labels = _labels_file(tmp_path, [("chr1", 10, "A", "T", "positive")])
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter(
                [
                    _classified("chr1", 10, "A", "T", ACMGClassification.PATHOGENIC),
                    _classified("chr9", 50, "G", "C", ACMGClassification.PATHOGENIC),
                ]
            )
        )
        assert report.matched == 1

    def test_confusion_matrix_cell(self, tmp_path: Path) -> None:
        labels = _labels_file(tmp_path, [("chr1", 10, "A", "T", "positive")])
        runner = BenchmarkRunner(labels)
        report = runner.evaluate(
            iter(
                [
                    _classified(
                        "chr1", 10, "A", "T", ACMGClassification.LIKELY_PATHOGENIC
                    )
                ]
            )
        )
        assert report.confusion["Likely_Pathogenic"]["Pathogenic"] == 1
