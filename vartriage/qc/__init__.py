"""VCF quality control metrics and validation.

Pre-flight data sanity checks computed before annotation. Detects
sample swaps, contamination, and variant caller artifacts via
population-level summary statistics (Ti/Tv, het/hom, variant counts).
"""

from vartriage.qc.config import QCConfig
from vartriage.qc.metrics import QCComputer, QCMetrics
from vartriage.qc.report import (
    format_clinical_qc_section,
    format_qc_stderr,
    qc_check_rows,
    serialize_qc_json,
)
from vartriage.qc.thresholds import AssayThresholds, AssayType
from vartriage.qc.validator import QCCheckResult, QCReport, QCStatus, QCValidator

__all__ = [
    "AssayThresholds",
    "AssayType",
    "QCCheckResult",
    "QCComputer",
    "QCConfig",
    "QCMetrics",
    "QCReport",
    "QCStatus",
    "QCValidator",
    "format_clinical_qc_section",
    "format_qc_stderr",
    "qc_check_rows",
    "serialize_qc_json",
]
