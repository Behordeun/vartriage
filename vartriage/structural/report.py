"""Structural variant findings section for clinical reports.

Generates the SV-specific portion of clinical reports, including
a findings table, per-SV evidence narrative, and summary statistics.
Designed to integrate with the existing clinical report generator
as an additional section alongside SNV findings.
"""

from __future__ import annotations

import html
from dataclasses import dataclass, field
from typing import Sequence

from vartriage.structural.models import (
    ClassifiedSV,
    SVClassification,
    SVConsequence,
)


@dataclass(frozen=True)
class SVFindingsRow:
    """One row in the SV findings table of the clinical report."""

    sv_type: str
    coordinates: str
    size: str
    consequence: str
    classification: str
    genes: str
    syndrome: str
    evidence_score: float
    pathogenicity_score: str


@dataclass(frozen=True)
class SVSummaryData:
    """Summary statistics for the SV findings section."""

    total_svs: int
    pathogenic_count: int
    likely_pathogenic_count: int
    vus_count: int
    syndromes_identified: list[str]
    hi_genes_affected: list[str]


@dataclass(frozen=True)
class SVReportSection:
    """Complete SV findings section for a clinical report.

    Parameters
    ----------
    findings : list[SVFindingsRow]
        Sorted table rows for all reportable SVs.
    summary : SVSummaryData
        Aggregate metrics for the section header.
    narratives : list[str]
        Per-SV evidence narratives in clinical language.
    """

    findings: list[SVFindingsRow] = field(default_factory=list)
    summary: SVSummaryData = field(
        default_factory=lambda: SVSummaryData(
            total_svs=0,
            pathogenic_count=0,
            likely_pathogenic_count=0,
            vus_count=0,
            syndromes_identified=[],
            hi_genes_affected=[],
        )
    )
    narratives: list[str] = field(default_factory=list)


# Classification sort priority (lower = more severe)
_SV_TIER_ORDER: dict[SVClassification, int] = {
    SVClassification.PATHOGENIC: 0,
    SVClassification.LIKELY_PATHOGENIC: 1,
    SVClassification.VUS: 2,
    SVClassification.LIKELY_BENIGN: 3,
    SVClassification.BENIGN: 4,
}


class SVReportBuilder:
    """Build the SV findings section for a clinical report.

    Takes a list of ClassifiedSV records and produces a structured
    SVReportSection with findings table, summary, and narratives.
    """

    def build(self, sv_results: Sequence[ClassifiedSV]) -> SVReportSection:
        """Assemble the SV report section from classified results.

        Parameters
        ----------
        sv_results : Sequence[ClassifiedSV]
            Classified SVs to include in the report.

        Returns
        -------
        SVReportSection
            Complete section data for template rendering.
        """
        if not sv_results:
            return SVReportSection()

        sorted_svs = sorted(
            sv_results,
            key=lambda s: (
                _SV_TIER_ORDER.get(s.classification, 99),
                -(s.scored.pathogenicity_score or 0.0),
            ),
        )

        findings = [self._build_row(sv) for sv in sorted_svs]
        summary = self._build_summary(sorted_svs)
        narratives = [self._build_narrative(sv) for sv in sorted_svs]

        return SVReportSection(
            findings=findings,
            summary=summary,
            narratives=narratives,
        )

    def _build_row(self, classified: ClassifiedSV) -> SVFindingsRow:
        """Convert one ClassifiedSV to a findings table row."""
        sv = classified.scored.annotated.sv
        annotated = classified.scored.annotated

        genes = ", ".join(
            o.gene_symbol for o in annotated.gene_overlaps[:5]
        )
        if annotated.genes_affected > 5:
            genes += f" (+{annotated.genes_affected - 5} more)"

        size_str = _format_size(sv.length)

        score_str = (
            f"{classified.scored.pathogenicity_score:.3f}"
            if classified.scored.pathogenicity_score is not None
            else "N/A"
        )

        return SVFindingsRow(
            sv_type=sv.sv_type.value,
            coordinates=f"{sv.chrom}:{sv.start}-{sv.end}",
            size=size_str,
            consequence=annotated.consequence.value,
            classification=classified.classification.value,
            genes=genes,
            syndrome=classified.syndrome_name or "",
            evidence_score=classified.evidence_score,
            pathogenicity_score=score_str,
        )

    def _build_summary(self, svs: list[ClassifiedSV]) -> SVSummaryData:
        """Compute aggregate statistics for the section."""
        path_count = sum(
            1 for s in svs if s.classification == SVClassification.PATHOGENIC
        )
        lp_count = sum(
            1 for s in svs
            if s.classification == SVClassification.LIKELY_PATHOGENIC
        )
        vus_count = sum(
            1 for s in svs if s.classification == SVClassification.VUS
        )

        syndromes = [
            s.syndrome_name for s in svs
            if s.syndrome_name is not None
        ]

        hi_genes: list[str] = []
        for sv in svs:
            for overlap in sv.scored.annotated.gene_overlaps:
                if overlap.is_haploinsufficient and overlap.gene_symbol not in hi_genes:
                    hi_genes.append(overlap.gene_symbol)

        return SVSummaryData(
            total_svs=len(svs),
            pathogenic_count=path_count,
            likely_pathogenic_count=lp_count,
            vus_count=vus_count,
            syndromes_identified=syndromes,
            hi_genes_affected=hi_genes,
        )

    def _build_narrative(self, classified: ClassifiedSV) -> str:
        """Generate a clinical evidence narrative for one SV.

        External-sourced strings (gene symbols, syndrome names) are
        HTML-escaped to prevent XSS if rendered in HTML reports.
        """
        sv = classified.scored.annotated.sv
        annotated = classified.scored.annotated
        parts: list[str] = []

        # Opening: describe the SV
        size_str = _format_size(sv.length)
        chrom_safe = html.escape(sv.chrom)
        parts.append(
            f"{sv.sv_type.value} at {chrom_safe}:{sv.start}-{sv.end} "
            f"({size_str})"
        )

        # Gene content
        if annotated.genes_affected == 0:
            parts.append("does not overlap any protein-coding genes.")
        elif annotated.genes_affected == 1:
            gene = annotated.gene_overlaps[0]
            overlap_pct = f"{gene.overlap_fraction * 100:.0f}%"
            gene_name_safe = html.escape(gene.gene_symbol)
            parts.append(
                f"overlaps {overlap_pct} of {gene_name_safe} "
                f"({gene.exons_affected}/{gene.total_exons} exons)."
            )
        else:
            parts.append(
                f"overlaps {annotated.genes_affected} protein-coding genes."
            )

        # Dosage sensitivity
        if annotated.hi_genes_affected > 0:
            hi_names = [
                html.escape(o.gene_symbol) for o in annotated.gene_overlaps
                if o.is_haploinsufficient
            ]
            parts.append(
                f"Haploinsufficient gene(s) affected: {', '.join(hi_names)}."
            )

        # Syndrome match
        if classified.syndrome_name:
            syndrome_safe = html.escape(classified.syndrome_name)
            parts.append(
                f"Matches known pathogenic region: {syndrome_safe}."
            )

        # Population frequency
        if annotated.frequency_unknown:
            parts.append("Not observed in gnomAD-SV (absent from population databases).")
        elif annotated.population_frequency is not None:
            af_pct = f"{annotated.population_frequency * 100:.3f}%"
            parts.append(f"Population frequency: {af_pct}.")

        # Classification
        parts.append(
            f"Classification: {classified.classification.value} "
            f"(evidence score: {classified.evidence_score:.2f})."
        )

        return " ".join(parts)


def _format_size(bp: int) -> str:
    """Format base-pair size into human-readable string."""
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.2f} Mb"
    if bp >= 1_000:
        return f"{bp / 1_000:.1f} kb"
    return f"{bp} bp"
