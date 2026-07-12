"""Clinical report generation subpackage.

Produces structured, auditable clinical variant reports from
ClassifiedVariant data. Supports HTML, PDF (via WeasyPrint),
and DOCX (via python-docx) output formats.
"""

from vartriage.reporting.clinical.audit import AuditTrailWriter
from vartriage.reporting.clinical.generator import ClinicalReportGenerator
from vartriage.reporting.clinical.narrative import EvidenceNarrativeBuilder
from vartriage.reporting.clinical.template_engine import ReportTemplateEngine

__all__ = [
    "AuditTrailWriter",
    "ClinicalReportGenerator",
    "EvidenceNarrativeBuilder",
    "ReportTemplateEngine",
]
