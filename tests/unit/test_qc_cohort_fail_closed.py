"""QC and cohort fail closed on absent data instead of scoring it PASS.

A sample with no genotype calls or no indels is not a clean sample: it is a
sample the metric could not evaluate, which is a warning, not a pass. And a
cohort config that keeps singletons while excluding doubletons is an
inversion that silently drops real recurrence, so it is rejected.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.models.cohort import CohortConfig
from vartriage.qc.config import QCConfig
from vartriage.qc.metrics import QCMetrics
from vartriage.qc.validator import QCStatus, QCValidator


def _metrics(**overrides) -> QCMetrics:  # noqa: ANN003
    base = dict(
        total_variants=0,
        snv_count=0,
        indel_count=0,
        insertion_count=0,
        deletion_count=0,
        transition_count=0,
        transversion_count=0,
        ti_tv_ratio=2.0,
        het_count=0,
        hom_alt_count=0,
        het_hom_ratio=None,
        ins_del_ratio=0.0,
        per_chrom_counts={},
    )
    base.update(overrides)
    return QCMetrics(**base)


def _check(metrics: QCMetrics, name: str) -> QCStatus:
    report = QCValidator(QCConfig(assay_type="wgs")).validate(metrics)
    for check in report.checks:
        if check.metric_name == name:
            return check.status
    raise AssertionError(f"no check named {name}")


class TestQCNoDataWarns:
    def test_no_genotype_calls_warns_not_passes(self) -> None:
        # metrics.finalize yields het_hom_ratio 0.0 (not None) when there are
        # no hom-alt calls, so the check runs and must warn on zero calls.
        status = _check(
            _metrics(het_count=0, hom_alt_count=0, het_hom_ratio=0.0),
            "Het/Hom Ratio",
        )
        assert status == QCStatus.WARN

    def test_no_indels_warns_not_passes(self) -> None:
        status = _check(_metrics(insertion_count=0, deletion_count=0), "Ins/Del Ratio")
        assert status == QCStatus.WARN


class TestCohortConfigInversion:
    def test_singletons_with_high_min_recurrence_rejected(self) -> None:
        with pytest.raises(ValueError, match="min_recurrence"):
            CohortConfig(
                sample_vcfs=[Path("a.vcf"), Path("b.vcf"), Path("c.vcf")],
                output_path=Path("out"),
                min_recurrence=3,
                include_singletons=True,
            )

    def test_singletons_with_min_recurrence_two_is_allowed(self) -> None:
        cfg = CohortConfig(
            sample_vcfs=[Path("a.vcf"), Path("b.vcf")],
            output_path=Path("out"),
            min_recurrence=2,
            include_singletons=True,
        )
        assert cfg.include_singletons

    def test_high_min_recurrence_without_singletons_is_allowed(self) -> None:
        cfg = CohortConfig(
            sample_vcfs=[Path("a.vcf"), Path("b.vcf"), Path("c.vcf")],
            output_path=Path("out"),
            min_recurrence=3,
            include_singletons=False,
        )
        assert cfg.min_recurrence == 3
