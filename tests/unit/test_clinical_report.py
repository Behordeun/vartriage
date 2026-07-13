"""Tests for the clinical report generation module.

Covers: models, narrative builder, template engine, audit trail writer,
and the clinical report generator orchestration.
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path
from typing import Any

import pytest

from vartriage.models.config import ClinicalReportConfig
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    ClinVarAssertion,
    EvidenceTag,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.reporting.clinical.audit import AuditTrailWriter
from vartriage.reporting.clinical.generator import ClinicalReportGenerator
from vartriage.reporting.clinical.models import (
    EvidenceCardData,
    ExecutiveSummaryData,
    FindingsRow,
    HeaderData,
    MethodologyData,
    ReportSections,
    SignOffData,
)
from vartriage.reporting.clinical.narrative import EvidenceNarrativeBuilder
from vartriage.reporting.clinical.template_engine import ReportTemplateEngine


# --- Fixtures ---


def _make_variant(
    chrom: str = "chr17",
    pos: int = 43094464,
    ref: str = "ACTG",
    alt: str = "A",
    consequence: FunctionalConsequence = FunctionalConsequence.FRAMESHIFT,
    af: float | None = None,
    clinvar: ClinVarAssertion | None = ClinVarAssertion.PATHOGENIC,
    cadd: float | None = 38.0,
    revel: float | None = 0.94,
    spliceai: float | None = None,
    composite_rank: float | None = 0.97,
    tags: frozenset[EvidenceTag] = frozenset(
        {EvidenceTag.PVS1, EvidenceTag.PM2, EvidenceTag.PP3, EvidenceTag.PP5}
    ),
    classification: ACMGClassification = ACMGClassification.PATHOGENIC,
    gene_name: str | None = "BRCA1",
) -> ClassifiedVariant:
    raw = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=99.0,
        filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=raw,
        consequence=consequence,
        allele_frequency=af,
        clinvar_assertion=clinvar,
        frequency_unknown=af is None,
        clinvar_unknown=clinvar is None,
        gene_name=gene_name,
    )
    scored = ScoredVariant(
        annotated=annotated,
        cadd_phred=cadd,
        cadd_normalized=cadd / 99.0 if cadd else None,
        revel_score=revel,
        spliceai_score=spliceai,
        composite_rank=composite_rank,
    )
    return ClassifiedVariant(
        scored=scored,
        evidence_tags=tags,
        classification=classification,
        missing_data_sources=frozenset(),
    )


@pytest.fixture
def pathogenic_variant() -> ClassifiedVariant:
    return _make_variant()


@pytest.fixture
def vus_variant() -> ClassifiedVariant:
    return _make_variant(
        chrom="chr11",
        pos=108289632,
        ref="C",
        alt="T",
        consequence=FunctionalConsequence.MISSENSE,
        af=0.00089,
        clinvar=ClinVarAssertion.VUS,
        cadd=24.1,
        revel=0.58,
        composite_rank=0.62,
        tags=frozenset({EvidenceTag.PP3}),
        classification=ACMGClassification.VUS,
        gene_name="ATM",
    )


@pytest.fixture
def sample_config() -> ClinicalReportConfig:
    return ClinicalReportConfig(
        patient_id="TEST-001",
        panel_name="Test Panel v1",
        output_format="clinical-html",
    )


@pytest.fixture
def sample_sections() -> ReportSections:
    return ReportSections(
        header=HeaderData(
            patient_id="TEST-001",
            panel_name="Test Panel v1",
            analysis_date="2025-07-13",
            pipeline_version="0.5.0",
        ),
        executive_summary=ExecutiveSummaryData(
            total_variants_analyzed=1000,
            variants_passed_filters=50,
            pathogenic_count=1,
            likely_pathogenic_count=2,
            vus_count=10,
        ),
        findings_table=[
            FindingsRow(
                gene_name="BRCA1",
                consequence="Frameshift",
                classification="Pathogenic",
                composite_rank=0.97,
                chromosome="chr17",
                position=43094464,
            ),
        ],
        evidence_cards=[
            EvidenceCardData(
                gene_name="BRCA1",
                consequence="Frameshift",
                allele_frequency_formatted="Absent from gnomAD",
                predictor_scores_formatted=["CADD 38.0", "REVEL 0.94"],
                clinvar_assertion="Pathogenic",
                inheritance_pattern=None,
                evidence_tags_with_explanations=[
                    "PVS1: null variant in gene with LOF mechanism",
                ],
                narrative="This frameshift variant introduces a premature stop codon.",
            ),
        ],
        limitations=["REVEL unavailable for 2 non-coding variants"],
        methodology=MethodologyData(
            pipeline_version="0.5.0",
            reference_files={"refs/gencode.gtf": "sha256:abc123"},
            classification_parameters={"max_af": "0.0001"},
            analysis_timestamp="2025-07-13T10:00:00Z",
        ),
        sign_off=SignOffData(),
    )


# --- ClinicalReportConfig tests ---


class TestClinicalReportConfig:
    def test_valid_config(self) -> None:
        cfg = ClinicalReportConfig(
            patient_id="PAT-001",
            panel_name="Panel",
            output_format="clinical-html",
        )
        assert cfg.patient_id == "PAT-001"
        assert cfg.report_template == "standard"

    def test_empty_patient_id_raises(self) -> None:
        with pytest.raises(ValueError, match="patient_id"):
            ClinicalReportConfig(
                patient_id="",
                panel_name="Panel",
                output_format="clinical-html",
            )

    def test_whitespace_panel_name_raises(self) -> None:
        with pytest.raises(ValueError, match="panel_name"):
            ClinicalReportConfig(
                patient_id="PAT-001",
                panel_name="   ",
                output_format="clinical-html",
            )


# --- Models tests ---


class TestClinicalModels:
    def test_header_data_frozen(self) -> None:
        h = HeaderData(
            patient_id="P1",
            panel_name="Panel",
            analysis_date="2025-01-01",
            pipeline_version="0.5.0",
        )
        with pytest.raises(Exception):
            h.patient_id = "P2"  # type: ignore[misc]

    def test_findings_row_none_gene(self) -> None:
        row = FindingsRow(
            gene_name=None,
            consequence="Intergenic",
            classification="VUS",
            composite_rank=None,
            chromosome="chr1",
            position=100,
        )
        assert row.gene_name is None

    def test_sign_off_defaults(self) -> None:
        s = SignOffData()
        assert s.reviewer_name_placeholder == "[Reviewer Name]"
        assert s.review_date_placeholder == "[Review Date]"


# --- Narrative builder tests ---


class TestEvidenceNarrativeBuilder:
    def test_build_narrative_pathogenic(
        self, pathogenic_variant: ClassifiedVariant
    ) -> None:
        builder = EvidenceNarrativeBuilder()
        narrative = builder.build_narrative(pathogenic_variant)
        assert len(narrative) > 50
        assert "BRCA1" in narrative or "frameshift" in narrative.lower()

    def test_build_narrative_vus(self, vus_variant: ClassifiedVariant) -> None:
        builder = EvidenceNarrativeBuilder()
        narrative = builder.build_narrative(vus_variant)
        assert len(narrative) > 20

    def test_format_allele_frequency(self) -> None:
        builder = EvidenceNarrativeBuilder()
        formatted = builder.format_allele_frequency(0.001)
        assert "0.001" in formatted or "1,000" in formatted or "1000" in formatted

    def test_format_allele_frequency_very_rare(self) -> None:
        builder = EvidenceNarrativeBuilder()
        formatted = builder.format_allele_frequency(0.000007)
        assert "0.000007" in formatted or "143" in formatted

    def test_format_predictor_score(self) -> None:
        builder = EvidenceNarrativeBuilder()
        formatted = builder.format_predictor_score("CADD", 38.0)
        assert "CADD" in formatted
        assert "38" in formatted

    def test_format_evidence_tag_pvs1(
        self, pathogenic_variant: ClassifiedVariant
    ) -> None:
        builder = EvidenceNarrativeBuilder()
        text = builder.format_evidence_tag(EvidenceTag.PVS1, pathogenic_variant)
        assert "PVS1" in text

    def test_format_evidence_tag_pm2(
        self, pathogenic_variant: ClassifiedVariant
    ) -> None:
        builder = EvidenceNarrativeBuilder()
        text = builder.format_evidence_tag(EvidenceTag.PM2, pathogenic_variant)
        assert "PM2" in text


# --- Template engine tests ---


class TestReportTemplateEngine:
    def test_render_html_produces_valid_structure(
        self, sample_sections: ReportSections
    ) -> None:
        engine = ReportTemplateEngine()
        html = engine.render_html(sample_sections)
        assert "<!DOCTYPE html>" in html
        assert "Clinical Variant Report" in html
        assert "TEST-001" in html
        assert "section-header" in html
        assert "section-executive-summary" in html
        assert "section-findings-table" in html
        assert "section-evidence-cards" in html
        assert "section-limitations" in html
        assert "section-methodology" in html
        assert "section-sign-off" in html

    def test_render_html_escapes_patient_id(self) -> None:
        sections = ReportSections(
            header=HeaderData(
                patient_id='<script>alert("xss")</script>',
                panel_name="Panel",
                analysis_date="2025-01-01",
                pipeline_version="0.5.0",
            ),
            executive_summary=ExecutiveSummaryData(
                total_variants_analyzed=0,
                variants_passed_filters=0,
                pathogenic_count=0,
                likely_pathogenic_count=0,
                vus_count=0,
            ),
            findings_table=[],
            evidence_cards=[],
            limitations=[],
            methodology=MethodologyData(
                pipeline_version="0.5.0",
                reference_files={},
                classification_parameters={},
                analysis_timestamp="2025-01-01T00:00:00Z",
            ),
        )
        engine = ReportTemplateEngine()
        html = engine.render_html(sections)
        assert "<script>" not in html
        assert "&lt;script&gt;" in html

    def test_render_html_sections_ordered(
        self, sample_sections: ReportSections
    ) -> None:
        engine = ReportTemplateEngine()
        html = engine.render_html(sample_sections)
        header_pos = html.index("section-header")
        summary_pos = html.index("section-executive-summary")
        findings_pos = html.index("section-findings-table")
        cards_pos = html.index("section-evidence-cards")
        assert header_pos < summary_pos < findings_pos < cards_pos

    def test_empty_findings_produces_valid_html(self) -> None:
        sections = ReportSections(
            header=HeaderData(
                patient_id="P1",
                panel_name="Panel",
                analysis_date="2025-01-01",
                pipeline_version="0.5.0",
            ),
            executive_summary=ExecutiveSummaryData(
                total_variants_analyzed=100,
                variants_passed_filters=0,
                pathogenic_count=0,
                likely_pathogenic_count=0,
                vus_count=0,
            ),
            findings_table=[],
            evidence_cards=[],
            limitations=[],
            methodology=MethodologyData(
                pipeline_version="0.5.0",
                reference_files={},
                classification_parameters={},
                analysis_timestamp="2025-01-01T00:00:00Z",
            ),
        )
        engine = ReportTemplateEngine()
        html = engine.render_html(sections)
        assert "<!DOCTYPE html>" in html
        assert "section-findings-table" in html


# --- Audit trail tests ---


class TestAuditTrailWriter:
    def test_write_creates_sidecar_file(
        self, tmp_path: Path, sample_config: ClinicalReportConfig
    ) -> None:
        report_path = tmp_path / "report.html"
        report_path.write_text("<html></html>")

        writer = AuditTrailWriter()
        sidecar = writer.write(
            output_path=report_path,
            config=sample_config,
            variants=[],
            reference_checksums={"refs/test.tsv": "sha256:abc"},
            pipeline_version="0.5.0",
            execution_timestamp="2025-07-13T10:00:00Z",
        )

        assert sidecar.exists()
        assert sidecar.name == "report.html.audit.json"

        data = json.loads(sidecar.read_text())
        assert "run_manifest" in data
        assert "decision_log" in data
        assert data["run_manifest"]["patient_id"] == "TEST-001"
        assert data["run_manifest"]["pipeline_version"] == "0.5.0"

    def test_write_restricted_permissions(
        self, tmp_path: Path, sample_config: ClinicalReportConfig
    ) -> None:
        report_path = tmp_path / "report.html"
        report_path.write_text("<html></html>")

        writer = AuditTrailWriter()
        sidecar = writer.write(
            output_path=report_path,
            config=sample_config,
            variants=[],
            reference_checksums={},
            pipeline_version="0.5.0",
            execution_timestamp="2025-07-13T10:00:00Z",
        )

        mode = stat.S_IMODE(os.stat(sidecar).st_mode)
        assert mode == 0o600

    def test_write_with_variants(
        self,
        tmp_path: Path,
        sample_config: ClinicalReportConfig,
        pathogenic_variant: ClassifiedVariant,
    ) -> None:
        report_path = tmp_path / "report.html"
        report_path.write_text("<html></html>")

        writer = AuditTrailWriter()
        sidecar = writer.write(
            output_path=report_path,
            config=sample_config,
            variants=[pathogenic_variant],
            reference_checksums={},
            pipeline_version="0.5.0",
            execution_timestamp="2025-07-13T10:00:00Z",
        )

        data = json.loads(sidecar.read_text())
        assert len(data["decision_log"]) == 1
        entry = data["decision_log"][0]
        assert entry["chromosome"] == "chr17"
        assert entry["position"] == 43094464
        assert entry["classification"] == "Pathogenic"
        assert "PVS1" in entry["evidence_tags_assigned"]

    def test_compute_file_checksum(self, tmp_path: Path) -> None:
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello world")

        writer = AuditTrailWriter()
        checksum = writer.compute_file_checksum(test_file)
        assert checksum.startswith("sha256:")
        assert len(checksum) == len("sha256:") + 64


# --- Generator integration tests ---


class TestClinicalReportGenerator:
    def test_generate_html_report(
        self,
        tmp_path: Path,
        sample_config: ClinicalReportConfig,
        pathogenic_variant: ClassifiedVariant,
        vus_variant: ClassifiedVariant,
    ) -> None:
        output_path = tmp_path / "test_report.html"
        generator = ClinicalReportGenerator(
            config=sample_config,
            pipeline_version="0.5.0",
            reference_checksums={"refs/test.tsv": "sha256:abc123"},
        )
        generator.generate(
            variants=[pathogenic_variant, vus_variant],
            output_path=output_path,
        )

        assert output_path.exists()
        html = output_path.read_text()
        assert "Clinical Variant Report" in html
        assert "TEST-001" in html
        assert "BRCA1" in html

        # Check audit sidecar was created
        sidecar = Path(str(output_path) + ".audit.json")
        assert sidecar.exists()

    def test_generate_restricts_file_permissions(
        self,
        tmp_path: Path,
        sample_config: ClinicalReportConfig,
        pathogenic_variant: ClassifiedVariant,
    ) -> None:
        output_path = tmp_path / "test_report.html"
        generator = ClinicalReportGenerator(
            config=sample_config,
            pipeline_version="0.5.0",
        )
        generator.generate(
            variants=[pathogenic_variant],
            output_path=output_path,
        )

        mode = stat.S_IMODE(os.stat(output_path).st_mode)
        assert mode == 0o600

    def test_generate_empty_variants(
        self, tmp_path: Path, sample_config: ClinicalReportConfig
    ) -> None:
        output_path = tmp_path / "empty_report.html"
        generator = ClinicalReportGenerator(
            config=sample_config,
            pipeline_version="0.5.0",
        )
        generator.generate(
            variants=[],
            output_path=output_path,
        )

        assert output_path.exists()
        html = output_path.read_text()
        assert "Pathogenic:</strong> 0" in html
        assert "Likely Pathogenic:</strong> 0" in html

    def test_generate_sorts_by_classification_tier(
        self,
        tmp_path: Path,
        sample_config: ClinicalReportConfig,
        pathogenic_variant: ClassifiedVariant,
        vus_variant: ClassifiedVariant,
    ) -> None:
        output_path = tmp_path / "sorted_report.html"
        generator = ClinicalReportGenerator(
            config=sample_config,
            pipeline_version="0.5.0",
        )
        # Pass VUS first, Pathogenic second
        generator.generate(
            variants=[vus_variant, pathogenic_variant],
            output_path=output_path,
        )

        html = output_path.read_text()
        # Pathogenic should appear before VUS in findings table
        pathogenic_pos = html.index("Pathogenic</td>")
        vus_pos = html.index("VUS</td>")
        assert pathogenic_pos < vus_pos
