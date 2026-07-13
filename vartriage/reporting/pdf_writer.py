"""PDF report renderer using reportlab (optional dependency).

This module implements the PDFRenderer protocol via the reportlab library,
generating clinical variant reports with a title section, variant table,
and page numbers.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import (PageBreak, Paragraph, SimpleDocTemplate,
                                    Spacer, Table, TableStyle)

    HAS_REPORTLAB = True
except ImportError:
    HAS_REPORTLAB = False

from vartriage.models.variant import ClassifiedVariant

_OUTPUT_FIELDS = [
    "Chromosome",
    "Position",
    "Ref Allele",
    "Alt Allele",
    "Functional Consequence",
    "Allele Frequency",
    "Composite Rank",
    "ClinVar Assertion",
    "ACMG Classification",
    "Evidence Tags",
]


def _format_value(value: object) -> str:
    """Format a field value for PDF display, representing None as 'N/A'."""
    if value is None:
        return "N/A"
    return str(value)


def _extract_row(variant: ClassifiedVariant) -> list[str]:
    """Extract a table row from a ClassifiedVariant.

    Parameters
    ----------
    variant : ClassifiedVariant
        The classified variant to extract fields from.

    Returns
    -------
    list[str]
        Row of formatted string values for the PDF table.
    """
    scored = variant.scored
    annotated = scored.annotated
    raw = annotated.variant

    consequence_str = (
        annotated.consequence.value if annotated.consequence is not None else None
    )
    af_str: Optional[str] = (
        f"{annotated.allele_frequency:.8g}"
        if annotated.allele_frequency is not None
        else None
    )
    rank_str: Optional[str] = (
        f"{scored.composite_rank:.4f}" if scored.composite_rank is not None else None
    )
    clinvar_str = (
        annotated.clinvar_assertion.value
        if annotated.clinvar_assertion is not None
        else None
    )
    classification_str = (
        variant.classification.value if variant.classification is not None else None
    )
    tags_str = (
        ", ".join(sorted(tag.value for tag in variant.evidence_tags))
        if variant.evidence_tags
        else None
    )

    return [
        _format_value(raw.chrom),
        _format_value(raw.pos),
        _format_value(raw.ref),
        _format_value(raw.alt),
        _format_value(consequence_str),
        _format_value(af_str),
        _format_value(rank_str),
        _format_value(clinvar_str),
        _format_value(classification_str),
        _format_value(tags_str),
    ]


class ReportlabPDFRenderer:
    """PDF renderer using reportlab to generate clinical variant reports.

    Generates a PDF document containing:
    - Title section with the report generation timestamp
    - Variant table with all output fields
    - Page numbers on each page

    Methods
    -------
    render(variants, output_path)
        Render classified variants to a PDF clinical report.

    Raises
    ------
    ImportError
        If reportlab is not installed.
    """

    def __init__(self) -> None:
        if not HAS_REPORTLAB:
            raise ImportError(
                "PDF output requires reportlab. "
                "Install with: pip install vartriage[pdf]"
            )

    def render(self, variants: list[Any], output_path: Path) -> Path:
        """Render classified variants to a PDF clinical report.

        Parameters
        ----------
        variants : list[Any]
            List of ClassifiedVariant instances to include in the report.
        output_path : Path
            Filesystem path where the PDF should be written.

        Returns
        -------
        Path
            The path to the generated PDF file.

        Raises
        ------
        IOError
            If the file cannot be written to the specified path.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        doc = SimpleDocTemplate(
            str(output_path),
            pagesize=landscape(A4),
            leftMargin=1.5 * cm,
            rightMargin=1.5 * cm,
            topMargin=2 * cm,
            bottomMargin=2 * cm,
        )

        styles = getSampleStyleSheet()
        elements: list[Any] = []

        timestamp = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        title = Paragraph(
            "Variant Prioritization Clinical Report",
            styles["Title"],
        )
        elements.append(title)
        elements.append(Spacer(1, 0.5 * cm))

        subtitle = Paragraph(
            f"Generated: {timestamp}",
            styles["Normal"],
        )
        elements.append(subtitle)
        elements.append(Spacer(1, 1 * cm))

        table_data: list[list[str]] = [_OUTPUT_FIELDS]
        for variant in variants:
            table_data.append(_extract_row(variant))

        col_widths = [
            2.0 * cm,  # Chromosome
            2.0 * cm,  # Position
            1.8 * cm,  # Ref Allele
            1.8 * cm,  # Alt Allele
            3.5 * cm,  # Functional Consequence
            2.5 * cm,  # Allele Frequency
            2.5 * cm,  # Composite Rank
            3.0 * cm,  # ClinVar Assertion
            3.0 * cm,  # ACMG Classification
            3.5 * cm,  # Evidence Tags
        ]

        table = Table(table_data, colWidths=col_widths, repeatRows=1)
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c3e50")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, 0), 7),
                    ("FONTSIZE", (0, 1), (-1, -1), 6),
                    ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
                    ("TOPPADDING", (0, 0), (-1, 0), 8),
                    ("BOTTOMPADDING", (0, 1), (-1, -1), 4),
                    ("TOPPADDING", (0, 1), (-1, -1), 4),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
                    (
                        "ROWBACKGROUNDS",
                        (0, 1),
                        (-1, -1),
                        [colors.white, colors.HexColor("#ecf0f1")],
                    ),
                ]
            )
        )

        elements.append(table)

        doc.build(elements, onFirstPage=_add_page_number, onLaterPages=_add_page_number)

        return output_path


def _add_page_number(canvas: object, doc: object) -> None:
    """Add page number footer to each page.

    Parameters
    ----------
    canvas : reportlab Canvas
        The canvas object for the current page.
    doc : reportlab BaseDocTemplate
        The document template being built.
    """
    page_num = canvas.getPageNumber()  # type: ignore[attr-defined]
    text = f"Page {page_num}"
    canvas.saveState()  # type: ignore[attr-defined]
    canvas.setFont("Helvetica", 8)  # type: ignore[attr-defined]
    canvas.drawCentredString(  # type: ignore[attr-defined]
        doc.pagesize[0] / 2.0,  # type: ignore[attr-defined]
        1.0 * cm,
        text,
    )
    canvas.restoreState()  # type: ignore[attr-defined]
