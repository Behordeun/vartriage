"""Tests for PDF report writer and fallback behavior."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from vartriage.models.variant import (ACMGClassification, AnnotatedVariant,
                                      ClassifiedVariant, ClinVarAssertion,
                                      EvidenceTag, FunctionalConsequence,
                                      ScoredVariant, Variant)
from vartriage.reporting.pdf_fallback import PDFFallbackRenderer


def _make_classified_variant(
    chrom: str = "chr1",
    pos: int = 12345,
    ref: str = "A",
    alt: str = "T",
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    allele_frequency: float | None = 0.0005,
    composite_rank: float | None = 0.85,
    clinvar_assertion: ClinVarAssertion | None = ClinVarAssertion.PATHOGENIC,
    classification: ACMGClassification = ACMGClassification.LIKELY_PATHOGENIC,
    evidence_tags: frozenset[EvidenceTag] | None = None,
) -> ClassifiedVariant:
    """Helper to build a ClassifiedVariant for testing."""
    if evidence_tags is None:
        evidence_tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PM2})

    raw = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
    )
    annotated = AnnotatedVariant(
        variant=raw,
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
    )
    scored = ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
    )
    return ClassifiedVariant(
        scored=scored,
        evidence_tags=evidence_tags,
        classification=classification,
    )


class TestPDFFallbackRenderer:
    """Tests for the fallback renderer that raises ImportError."""

    def test_render_raises_import_error(self) -> None:
        renderer = PDFFallbackRenderer()
        with pytest.raises(ImportError, match="PDF output requires reportlab"):
            renderer.render([], Path("/tmp/test.pdf"))

    def test_error_message_includes_install_instructions(self) -> None:
        renderer = PDFFallbackRenderer()
        with pytest.raises(ImportError, match="pip install vartriage\\[pdf\\]"):
            renderer.render([], Path("/tmp/test.pdf"))


class TestReportlabPDFRenderer:
    """Tests for actual PDF generation using reportlab (skipped if not installed)."""

    @pytest.fixture(autouse=True)
    def _check_reportlab(self) -> None:
        pytest.importorskip("reportlab")

    def test_renders_empty_variant_list(self, tmp_path: Path) -> None:
        from vartriage.reporting.pdf_writer import ReportlabPDFRenderer

        renderer = ReportlabPDFRenderer()
        output = tmp_path / "empty_report.pdf"
        result = renderer.render([], output)

        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_renders_single_variant(self, tmp_path: Path) -> None:
        from vartriage.reporting.pdf_writer import ReportlabPDFRenderer

        renderer = ReportlabPDFRenderer()
        variant = _make_classified_variant()
        output = tmp_path / "single_variant.pdf"
        result = renderer.render([variant], output)

        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_renders_multiple_variants(self, tmp_path: Path) -> None:
        from vartriage.reporting.pdf_writer import ReportlabPDFRenderer

        renderer = ReportlabPDFRenderer()
        variants = [
            _make_classified_variant(chrom="chr1", pos=100),
            _make_classified_variant(chrom="chr2", pos=200, allele_frequency=None),
            _make_classified_variant(
                chrom="chr3",
                pos=300,
                composite_rank=None,
                clinvar_assertion=None,
                evidence_tags=frozenset(),
                classification=ACMGClassification.VUS,
            ),
        ]
        output = tmp_path / "multi_variant.pdf"
        result = renderer.render(variants, output)

        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_null_values_rendered_as_na(self, tmp_path: Path) -> None:
        from vartriage.reporting.pdf_writer import (ReportlabPDFRenderer,
                                                    _extract_row)

        variant = _make_classified_variant(
            allele_frequency=None,
            composite_rank=None,
            clinvar_assertion=None,
            evidence_tags=frozenset(),
        )
        row = _extract_row(variant)

        assert row[6] == "N/A"  # Allele Frequency
        assert row[8] == "N/A"  # Composite Rank
        assert row[9] == "N/A"  # ClinVar Assertion
        assert row[11] == "N/A"  # Evidence Tags

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        from vartriage.reporting.pdf_writer import ReportlabPDFRenderer

        renderer = ReportlabPDFRenderer()
        output = tmp_path / "nested" / "dir" / "report.pdf"
        result = renderer.render([], output)

        assert result == output
        assert output.exists()

    def test_pdf_contains_expected_content(self, tmp_path: Path) -> None:
        """Verify PDF file has valid PDF header bytes."""
        from vartriage.reporting.pdf_writer import ReportlabPDFRenderer

        renderer = ReportlabPDFRenderer()
        variant = _make_classified_variant()
        output = tmp_path / "check_content.pdf"
        renderer.render([variant], output)

        content = output.read_bytes()
        assert content[:5] == b"%PDF-"

    def test_renderer_raises_without_reportlab(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Verify ReportlabPDFRenderer raises at construction if reportlab missing."""
        import vartriage.reporting.pdf_writer as pdf_module

        monkeypatch.setattr(pdf_module, "HAS_REPORTLAB", False)
        with pytest.raises(ImportError, match="PDF output requires reportlab"):
            pdf_module.ReportlabPDFRenderer()
