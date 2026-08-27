"""QC report formatting: stderr table, JSON, and clinical report section."""

from __future__ import annotations

import json
import sys
from typing import Any

from vartriage.qc.validator import QCCheckResult, QCReport, QCStatus

_STATUS_SYMBOLS = {
    QCStatus.PASS: "\u2713 PASS",
    QCStatus.WARN: "\u26a0 WARN",
    QCStatus.FAIL: "\u2717 FAIL",
}


def format_qc_stderr(report: QCReport) -> str:
    """Format QC report as a human-readable table for stderr output.

    Parameters
    ----------
    report : QCReport
        Validated QC report to format.

    Returns
    -------
    str
        Multi-line formatted table with metric/value/expected/status columns.
    """
    lines: list[str] = []
    lines.append("")
    lines.append("=" * 70)
    lines.append(f"  VCF Quality Control  |  Assay: {report.assay_type.upper()}")
    lines.append("=" * 70)
    lines.append(f"  {'Metric':<20} {'Value':>12} {'Expected':>16} {'Status':>10}")
    lines.append("  " + "-" * 62)

    for check in report.checks:
        value_str = _format_value(check)
        expected_str = f"{check.expected_min:.1f}-{check.expected_max:.1f}"
        # Special formatting for variant count
        if check.metric_name == "Total Variants":
            value_str = f"{int(check.value):,}"
            if check.expected_max >= 999_999_999:
                expected_str = "no constraint"
            else:
                expected_str = (
                    f"{int(check.expected_min):,}-{int(check.expected_max):,}"
                )
        status_str = _STATUS_SYMBOLS[check.status]
        lines.append(
            f"  {check.metric_name:<20} {value_str:>12} {expected_str:>16} {status_str:>10}"
        )

    lines.append("  " + "-" * 62)
    overall_sym = _STATUS_SYMBOLS[report.overall_status]
    lines.append(f"  Overall QC Verdict: {overall_sym}")
    lines.append("=" * 70)
    lines.append("")

    return "\n".join(lines)


def print_qc_stderr(report: QCReport) -> None:
    """Print the QC report table to stderr."""
    sys.stderr.write(format_qc_stderr(report))
    sys.stderr.flush()


def serialize_qc_json(report: QCReport) -> dict[str, Any]:
    """Serialize QCReport to a JSON-compatible dictionary.

    Parameters
    ----------
    report : QCReport
        Validated QC report to serialize.

    Returns
    -------
    dict[str, Any]
        JSON-serializable representation.
    """
    metrics = report.metrics
    return {
        "assay_type": report.assay_type,
        "overall_status": report.overall_status.value,
        "metrics": {
            "total_variants": metrics.total_variants,
            "snv_count": metrics.snv_count,
            "indel_count": metrics.indel_count,
            "insertion_count": metrics.insertion_count,
            "deletion_count": metrics.deletion_count,
            "transition_count": metrics.transition_count,
            "transversion_count": metrics.transversion_count,
            "ti_tv_ratio": round(metrics.ti_tv_ratio, 4),
            "het_count": metrics.het_count,
            "hom_alt_count": metrics.hom_alt_count,
            "het_hom_ratio": (
                round(metrics.het_hom_ratio, 4)
                if metrics.het_hom_ratio is not None
                else None
            ),
            "ins_del_ratio": round(metrics.ins_del_ratio, 4),
            "per_chrom_counts": metrics.per_chrom_counts,
        },
        "checks": [_check_to_dict(c) for c in report.checks],
    }


def write_qc_json(report: QCReport, path: Any) -> None:
    """Write QC report as JSON to a file path.

    Parameters
    ----------
    report : QCReport
        Validated QC report.
    path : Path
        Output file path.
    """
    from pathlib import Path

    resolved = Path(path)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    with open(resolved, "w") as f:
        json.dump(serialize_qc_json(report), f, indent=2)


def format_clinical_qc_section(report: QCReport) -> str:
    """Generate HTML for the Sample Quality Control section of clinical reports.

    Parameters
    ----------
    report : QCReport
        Validated QC report.

    Returns
    -------
    str
        HTML fragment for the clinical report.
    """
    rows: list[str] = []
    for check in report.checks:
        value_str = _format_value(check)
        if check.metric_name == "Total Variants":
            value_str = f"{int(check.value):,}"
            expected = (
                "no constraint"
                if check.expected_max >= 999_999_999
                else f"{int(check.expected_min):,}-{int(check.expected_max):,}"
            )
        else:
            expected = f"{check.expected_min:.1f}-{check.expected_max:.1f}"
        status_sym = _STATUS_SYMBOLS[check.status]
        rows.append(
            f"    <tr><td>{check.metric_name}</td>"
            f"<td>{value_str}</td>"
            f"<td>{expected}</td>"
            f"<td>{status_sym}</td></tr>"
        )

    return (
        '<section id="sample-qc">\n'
        "  <h2>Sample Quality Control</h2>\n"
        "  <table>\n"
        "    <tr><th>Metric</th><th>Value</th><th>Expected</th><th>Status</th></tr>\n"
        + "\n".join(rows)
        + "\n  </table>\n"
        "</section>"
    )


def _format_value(check: QCCheckResult) -> str:
    """Format a check value for display."""
    if check.metric_name == "Total Variants":
        return f"{int(check.value):,}"
    return f"{check.value:.2f}"


def _check_to_dict(check: QCCheckResult) -> dict[str, Any]:
    """Convert a single check result to a dict for JSON serialization."""
    return {
        "metric_name": check.metric_name,
        "value": round(check.value, 4),
        "expected_min": check.expected_min,
        "expected_max": check.expected_max,
        "status": check.status.value,
        "message": check.message,
    }
