"""Clinical report generator orchestrating all report sections.

Materializes classified variants, assembles section data, delegates
rendering to the template engine, and writes the audit trail sidecar.
Uses atomic file writing (temp file + rename) so the target path never
contains partial output.
"""

from __future__ import annotations

import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Sequence, Union

from vartriage.models.config import ClinicalReportConfig
from vartriage.models.variant import ACMGClassification, ClassifiedVariant
from vartriage.reporting.clinical.audit import AuditTrailWriter
from vartriage.reporting.clinical.models import (EvidenceCardData,
                                                 ExecutiveSummaryData,
                                                 FindingsRow, HeaderData,
                                                 MethodologyData,
                                                 ReportSections, SignOffData)
from vartriage.reporting.clinical.narrative import EvidenceNarrativeBuilder
from vartriage.reporting.clinical.template_engine import ReportTemplateEngine

# Tier ordering for findings table sort. Lower value = higher priority.
_TIER_ORDER: dict[ACMGClassification, int] = {
    ACMGClassification.PATHOGENIC: 0,
    ACMGClassification.LIKELY_PATHOGENIC: 1,
    ACMGClassification.VUS: 2,
    ACMGClassification.LIKELY_BENIGN: 3,
    ACMGClassification.BENIGN: 4,
}


class ClinicalReportGenerator:
    """Orchestrates clinical report generation across all sections.

    Materializes the variant iterator, sorts variants by clinical
    priority, assembles all report sections, renders via the template
    engine, and writes the audit trail sidecar.

    Parameters
    ----------
    config : ClinicalReportConfig
        Clinical report configuration with patient and format info.
    pipeline_version : str
        Version string of the vartriage pipeline.
    reference_checksums : dict[str, str]
        Mapping of reference file paths to SHA-256 checksums.
        Defaults to an empty dict.
    """

    def __init__(
        self,
        config: ClinicalReportConfig,
        pipeline_version: str,
        reference_checksums: dict[str, str] | None = None,
    ) -> None:
        self._config = config
        self._pipeline_version = pipeline_version
        self._reference_checksums = (
            reference_checksums if reference_checksums is not None else {}
        )
        self._narrative_builder = EvidenceNarrativeBuilder()
        self._template_engine = ReportTemplateEngine()
        self._audit_writer = AuditTrailWriter()

    def generate(
        self,
        variants: Union[
            Iterator[ClassifiedVariant],
            Sequence[ClassifiedVariant],
        ],
        output_path: Path,
    ) -> Path:
        """Generate a clinical report with audit trail.

        Materializes the variant iterator, assembles all sections,
        renders to the configured format, and writes the audit
        sidecar. Uses atomic file writing (temp + rename).

        Parameters
        ----------
        variants : Iterator or Sequence of ClassifiedVariant
            Classified variants to include in the report. May be
            empty.
        output_path : Path
            Destination path for the final report file.

        Returns
        -------
        Path
            The path where the report was written.

        Raises
        ------
        IOError
            On filesystem write failure or rendering errors.
        ImportError
            If a required rendering backend (WeasyPrint or
            python-docx) is not installed.
        """
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Materialize the variant iterator for random access.
        variant_list: list[ClassifiedVariant] = list(variants)

        # Sort for findings table display.
        sorted_variants = self._sort_variants(variant_list)

        # Generate execution timestamp once for consistency.
        timestamp = datetime.now(timezone.utc).isoformat()

        # Assemble all report sections.
        sections = self._assemble_sections(sorted_variants, timestamp)

        # Render and write atomically.
        self._render_atomic(sections, output_path)

        # Write audit trail sidecar. If this fails after the report
        # was written, remove the report to avoid an inconsistent
        # state (report exists without its audit trail).
        try:
            self._audit_writer.write(
                output_path=output_path,
                config=self._config,
                variants=sorted_variants,
                reference_checksums=self._reference_checksums,
                pipeline_version=self._pipeline_version,
                execution_timestamp=timestamp,
            )
        except Exception as exc:
            try:
                output_path.unlink(missing_ok=True)
            except OSError:
                pass
            raise IOError(f"Audit trail write failed: {exc}") from exc

        return output_path

    def _sort_variants(
        self, variants: list[ClassifiedVariant]
    ) -> list[ClassifiedVariant]:
        """Sort variants by tier priority then composite rank.

        Primary sort: classification tier ascending (Pathogenic
        first). Secondary sort: composite_rank descending within
        the same tier. Variants with None composite_rank sort last
        within their tier.
        """

        def sort_key(
            v: ClassifiedVariant,
        ) -> tuple[int, float]:
            tier = _TIER_ORDER.get(v.classification, 99)
            rank = v.scored.composite_rank
            # Negate rank for descending order. None sorts last.
            rank_key = -rank if rank is not None else float("inf")
            return (tier, rank_key)

        return sorted(variants, key=sort_key)

    def _assemble_sections(
        self,
        variants: list[ClassifiedVariant],
        timestamp: str,
    ) -> ReportSections:
        """Build all report section data from sorted variants."""
        header = self._build_header(timestamp)
        executive_summary = self._build_executive_summary(variants)
        findings_table = self._build_findings_table(variants)
        evidence_cards = self._build_evidence_cards(variants)
        limitations = self._collect_limitations(variants)
        methodology = self._build_methodology(timestamp)
        sign_off = SignOffData()

        return ReportSections(
            header=header,
            executive_summary=executive_summary,
            findings_table=findings_table,
            evidence_cards=evidence_cards,
            limitations=limitations,
            methodology=methodology,
            sign_off=sign_off,
        )

    def _build_header(self, timestamp: str) -> HeaderData:
        """Build the report header section."""
        return HeaderData(
            patient_id=self._config.patient_id,
            panel_name=self._config.panel_name,
            analysis_date=timestamp,
            pipeline_version=self._pipeline_version,
        )

    def _build_executive_summary(
        self, variants: list[ClassifiedVariant]
    ) -> ExecutiveSummaryData:
        """Build executive summary with variant counts per tier."""
        total = len(variants)
        # All variants in this list have passed pipeline filters.
        passed = total

        pathogenic = sum(
            1 for v in variants if v.classification == ACMGClassification.PATHOGENIC
        )
        likely_pathogenic = sum(
            1
            for v in variants
            if v.classification == ACMGClassification.LIKELY_PATHOGENIC
        )
        vus = sum(1 for v in variants if v.classification == ACMGClassification.VUS)

        return ExecutiveSummaryData(
            total_variants_analyzed=total,
            variants_passed_filters=passed,
            pathogenic_count=pathogenic,
            likely_pathogenic_count=likely_pathogenic,
            vus_count=vus,
        )

    def _build_findings_table(
        self, variants: list[ClassifiedVariant]
    ) -> list[FindingsRow]:
        """Build findings table rows from sorted variants."""
        rows: list[FindingsRow] = []
        for v in variants:
            annotated = v.scored.annotated
            rows.append(
                FindingsRow(
                    gene_name=annotated.gene_name,
                    consequence=annotated.consequence.value,
                    classification=v.classification.value,
                    composite_rank=v.scored.composite_rank,
                    chromosome=annotated.variant.chrom,
                    position=annotated.variant.pos,
                )
            )
        return rows

    def _build_evidence_cards(
        self, variants: list[ClassifiedVariant]
    ) -> list[EvidenceCardData]:
        """Build one evidence card per variant."""
        cards: list[EvidenceCardData] = []

        for v in variants:
            annotated = v.scored.annotated
            scored = v.scored

            # Format allele frequency if available.
            af_formatted = None
            if annotated.allele_frequency is not None:
                af_formatted = self._narrative_builder.format_allele_frequency(
                    annotated.allele_frequency
                )

            # Format predictor scores.
            scores_formatted: list[str] = []
            for name, value in [
                ("REVEL", scored.revel_score),
                ("CADD", scored.cadd_phred),
                ("SpliceAI", scored.spliceai_score),
            ]:
                if value is not None:
                    scores_formatted.append(
                        self._narrative_builder.format_predictor_score(name, value)
                    )

            # ClinVar assertion.
            clinvar = None
            if annotated.clinvar_assertion is not None:
                clinvar = annotated.clinvar_assertion.value

            # Inheritance pattern from variant info.
            inheritance = annotated.variant.info.get("inheritance_pattern")

            # Evidence tags with explanations.
            tags_explained: list[str] = []
            for tag in sorted(v.evidence_tags, key=lambda t: t.value):
                tags_explained.append(
                    self._narrative_builder.format_evidence_tag(tag, v)
                )

            # Full narrative text.
            narrative = self._narrative_builder.build_narrative(v)

            cards.append(
                EvidenceCardData(
                    gene_name=annotated.gene_name,
                    consequence=annotated.consequence.value,
                    allele_frequency_formatted=af_formatted,
                    predictor_scores_formatted=scores_formatted,
                    clinvar_assertion=clinvar,
                    inheritance_pattern=inheritance,
                    evidence_tags_with_explanations=(tags_explained),
                    narrative=narrative,
                )
            )

        return cards

    def _collect_limitations(self, variants: list[ClassifiedVariant]) -> list[str]:
        """Collect all missing data source names for limitations."""
        all_sources: set[str] = set()
        for v in variants:
            all_sources.update(v.missing_data_sources)

        limitations: list[str] = []
        for source in sorted(all_sources):
            limitations.append(
                f"{source} data was not available during "
                f"classification, which may affect evidence "
                f"tag assignment for affected variants."
            )
        return limitations

    def _build_methodology(self, timestamp: str) -> MethodologyData:
        """Build the methodology section from config metadata."""
        # NOTE: classification_parameters currently only includes
        # report-level settings. Adding scoring thresholds (e.g.
        # max_allele_frequency, PP3 threshold) requires access to
        # PrioritizationConfig which is not available here.
        # Planned for v0.6.0.
        params: dict[str, str] = {
            "report_template": self._config.report_template,
            "output_format": self._config.output_format,
        }

        return MethodologyData(
            pipeline_version=self._pipeline_version,
            reference_files=dict(self._reference_checksums),
            classification_parameters=params,
            analysis_timestamp=timestamp,
        )

    def _render_atomic(self, sections: ReportSections, output_path: Path) -> None:
        """Render the report to a temp file, then atomically rename.

        This ensures the target path never holds partial output.
        Clinical reports may contain PHI (patient identifiers); file
        permissions are restricted to owner-only (0o600) to limit
        exposure on shared filesystems.
        """
        fmt = self._config.output_format
        tmp_fd = None
        tmp_path: Path | None = None

        # Resolve output_path to prevent path traversal
        output_path = output_path.resolve()
        if not output_path.parent.exists():
            raise IOError(f"Parent directory does not exist: {output_path.parent}")

        try:
            tmp_fd, tmp_name = tempfile.mkstemp(
                dir=output_path.parent,
                prefix=".clinical_report_",
                suffix=".tmp",
            )
            os.close(tmp_fd)
            tmp_fd = None
            tmp_path = Path(tmp_name)

            # Restrict permissions before writing PHI content
            os.chmod(tmp_path, 0o600)

            if fmt == "clinical-html":
                html_content = self._template_engine.render_html(sections)
                tmp_path.write_text(html_content, encoding="utf-8")
            elif fmt == "clinical-pdf":
                self._template_engine.render_pdf(sections, tmp_path)
            elif fmt == "clinical-docx":
                self._template_engine.render_docx(sections, tmp_path)
            else:
                raise IOError(f"Unsupported clinical format: {fmt}")

            os.replace(str(tmp_path), str(output_path))
            # Ensure final file also has restricted permissions
            os.chmod(output_path, 0o600)
            tmp_path = None

        except (IOError, ImportError):
            raise
        except Exception as exc:
            raise IOError(f"Failed to generate clinical report: {exc}") from exc
        finally:
            if tmp_path is not None:
                try:
                    tmp_path.unlink(missing_ok=True)
                except OSError:
                    pass
