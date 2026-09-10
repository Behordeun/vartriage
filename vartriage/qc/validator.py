"""QC validation: compare computed metrics against assay thresholds."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from vartriage.qc.config import QCConfig
from vartriage.qc.metrics import QCMetrics
from vartriage.qc.thresholds import get_thresholds


class QCStatus(Enum):
    """QC check outcome."""

    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


@dataclass(frozen=True)
class QCCheckResult:
    """Result of a single QC metric check."""

    metric_name: str
    value: float
    expected_min: float
    expected_max: float
    status: QCStatus
    message: str


@dataclass(frozen=True)
class QCReport:
    """Aggregated QC report across all metric checks."""

    metrics: QCMetrics
    checks: list[QCCheckResult]
    overall_status: QCStatus
    assay_type: str


class QCValidator:
    """Compare QCMetrics against assay-specific thresholds.

    Produces a QCReport with per-metric pass/warn/fail status and
    an overall verdict (worst status across all checks).

    Parameters
    ----------
    config : QCConfig
        QC configuration with assay type and optional threshold overrides.
    """

    def __init__(self, config: QCConfig) -> None:
        self._config = config
        self._thresholds = get_thresholds(config.assay_type)

    def validate(self, metrics: QCMetrics) -> QCReport:
        """Run all checks and produce a QCReport.

        Parameters
        ----------
        metrics : QCMetrics
            Computed metrics from a QC pass.

        Returns
        -------
        QCReport
            Frozen report with per-metric results and overall verdict.
        """
        checks: list[QCCheckResult] = []

        checks.append(self._check_ti_tv(metrics))
        checks.append(self._check_variant_count(metrics))
        checks.append(self._check_ins_del(metrics))

        if metrics.het_hom_ratio is not None:
            checks.append(self._check_het_hom(metrics))

        overall = QCStatus.PASS
        for check in checks:
            if check.status == QCStatus.FAIL:
                overall = QCStatus.FAIL
                break
            if check.status == QCStatus.WARN:
                overall = QCStatus.WARN

        return QCReport(
            metrics=metrics,
            checks=checks,
            overall_status=overall,
            assay_type=self._config.assay_type,
        )

    def _check_ti_tv(self, metrics: QCMetrics) -> QCCheckResult:
        """Check Ti/Tv ratio against thresholds."""
        value = metrics.ti_tv_ratio
        thresholds = self._thresholds

        # Allow CLI override
        warn_range = self._config.expected_ti_tv or thresholds.ti_tv_warn
        fail_range = thresholds.ti_tv_fail

        if value < fail_range[0] or value > fail_range[1]:
            status = QCStatus.FAIL
            msg = (
                f"Ti/Tv ratio {value:.2f} is critically out of range "
                f"(expected {fail_range[0]}-{fail_range[1]})"
            )
        elif value < warn_range[0] or value > warn_range[1]:
            status = QCStatus.WARN
            msg = (
                f"Ti/Tv ratio {value:.2f} is borderline "
                f"(expected {warn_range[0]}-{warn_range[1]})"
            )
        else:
            status = QCStatus.PASS
            msg = f"Ti/Tv ratio {value:.2f} within expected range"

        return QCCheckResult(
            metric_name="Ti/Tv Ratio",
            value=value,
            expected_min=warn_range[0],
            expected_max=warn_range[1],
            status=status,
            message=msg,
        )

    def _check_het_hom(self, metrics: QCMetrics) -> QCCheckResult:
        """Check het/hom ratio against thresholds."""
        value = metrics.het_hom_ratio or 0.0
        thresholds = self._thresholds

        warn_range = self._config.expected_het_hom or thresholds.het_hom_warn
        fail_range = thresholds.het_hom_fail

        # No genotype calls for this sample: ratio is undefined, not a signal
        if metrics.het_count == 0 and metrics.hom_alt_count == 0:
            return QCCheckResult(
                metric_name="Het/Hom Ratio",
                value=value,
                expected_min=warn_range[0],
                expected_max=warn_range[1],
                status=QCStatus.PASS,
                message="No genotype calls present; het/hom ratio not evaluated",
            )

        if value < fail_range[0] or value > fail_range[1]:
            status = QCStatus.FAIL
            msg = (
                f"Het/Hom ratio {value:.2f} is critically out of range "
                f"(expected {fail_range[0]}-{fail_range[1]})"
            )
        elif value < warn_range[0] or value > warn_range[1]:
            status = QCStatus.WARN
            msg = (
                f"Het/Hom ratio {value:.2f} is borderline "
                f"(expected {warn_range[0]}-{warn_range[1]})"
            )
        else:
            status = QCStatus.PASS
            msg = f"Het/Hom ratio {value:.2f} within expected range"

        return QCCheckResult(
            metric_name="Het/Hom Ratio",
            value=value,
            expected_min=warn_range[0],
            expected_max=warn_range[1],
            status=status,
            message=msg,
        )

    def _check_variant_count(self, metrics: QCMetrics) -> QCCheckResult:
        """Check total variant count against assay-type expected range."""
        value = float(metrics.total_variants)
        thresholds = self._thresholds
        expected_min = float(thresholds.variant_count_min)
        expected_max = float(thresholds.variant_count_max)

        if expected_min == 0 and expected_max == 999_999_999:
            # Panel mode with no count constraint
            return QCCheckResult(
                metric_name="Total Variants",
                value=value,
                expected_min=expected_min,
                expected_max=expected_max,
                status=QCStatus.PASS,
                message=f"Total variant count {int(value):,} (no range constraint for panel)",
            )

        # FAIL if >2x or <0.5x expected range
        fail_min = expected_min * 0.5
        fail_max = expected_max * 2.0

        if value < fail_min or value > fail_max:
            status = QCStatus.FAIL
            msg = (
                f"Total variant count {int(value):,} is critically "
                f"outside expected range ({int(expected_min):,}-{int(expected_max):,})"
            )
        elif value < expected_min or value > expected_max:
            status = QCStatus.WARN
            msg = (
                f"Total variant count {int(value):,} is outside "
                f"expected range ({int(expected_min):,}-{int(expected_max):,})"
            )
        else:
            status = QCStatus.PASS
            msg = (
                f"Total variant count {int(value):,} within expected range "
                f"({int(expected_min):,}-{int(expected_max):,})"
            )

        return QCCheckResult(
            metric_name="Total Variants",
            value=value,
            expected_min=expected_min,
            expected_max=expected_max,
            status=status,
            message=msg,
        )

    def _check_ins_del(self, metrics: QCMetrics) -> QCCheckResult:
        """Check insertion/deletion ratio against thresholds."""
        value = metrics.ins_del_ratio
        thresholds = self._thresholds
        warn_range = thresholds.ins_del_warn
        fail_range = thresholds.ins_del_fail

        # No indels present: ratio is undefined, not a failure signal
        if metrics.insertion_count == 0 and metrics.deletion_count == 0:
            return QCCheckResult(
                metric_name="Ins/Del Ratio",
                value=value,
                expected_min=warn_range[0],
                expected_max=warn_range[1],
                status=QCStatus.PASS,
                message="No indels present; ins/del ratio not evaluated",
            )

        if value < fail_range[0] or value > fail_range[1]:
            status = QCStatus.FAIL
            msg = (
                f"Ins/Del ratio {value:.2f} is critically out of range "
                f"(expected {fail_range[0]}-{fail_range[1]})"
            )
        elif value < warn_range[0] or value > warn_range[1]:
            status = QCStatus.WARN
            msg = (
                f"Ins/Del ratio {value:.2f} is borderline "
                f"(expected {warn_range[0]}-{warn_range[1]})"
            )
        else:
            status = QCStatus.PASS
            msg = f"Ins/Del ratio {value:.2f} within expected range"

        return QCCheckResult(
            metric_name="Ins/Del Ratio",
            value=value,
            expected_min=warn_range[0],
            expected_max=warn_range[1],
            status=status,
            message=msg,
        )
