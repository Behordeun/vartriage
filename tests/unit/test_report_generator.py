"""Smoke tests for the ReportGenerator orchestrator."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from vartriage.models.config import ReportConfig
from vartriage.models.variant import (ACMGClassification, AnnotatedVariant,
                                      ClassifiedVariant, ClinVarAssertion,
                                      EvidenceTag, FunctionalConsequence,
                                      ScoredVariant, Variant)
from vartriage.reporting.generator import ReportGenerator


def _make_variant(
    chrom: str = "chr1",
    pos: int = 12345,
    ref: str = "A",
    alt: str = "T",
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    allele_frequency: float | None = 0.001,
    composite_rank: float | None = 0.85,
    clinvar: ClinVarAssertion | None = ClinVarAssertion.PATHOGENIC,
    classification: ACMGClassification = ACMGClassification.LIKELY_PATHOGENIC,
    tags: frozenset[EvidenceTag] | None = None,
) -> ClassifiedVariant:
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
        clinvar_assertion=clinvar,
    )
    scored = ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
    )
    return ClassifiedVariant(
        scored=scored,
        evidence_tags=tags or frozenset(),
        classification=classification,
    )


class TestReportGeneratorJSON:
    """Tests for JSON output format routing."""

    def test_generates_json_with_variants(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        variants = [_make_variant(), _make_variant(pos=67890)]
        output = tmp_path / "report.json"

        result = gen.generate(variants, output)

        assert result == output
        assert output.exists()
        data = json.loads(output.read_text(encoding="utf-8"))
        assert len(data) == 2
        assert data[0]["chromosome"] == "chr1"
        assert data[0]["position"] == 12345
        assert data[1]["position"] == 67890

    def test_generates_json_empty_list(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        output = tmp_path / "empty.json"

        result = gen.generate([], output)

        assert result == output
        assert output.exists()
        data = json.loads(output.read_text(encoding="utf-8"))
        assert data == []

    def test_no_partial_file_on_write_failure(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from vartriage.reporting import generator

        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        target = tmp_path / "report.json"

        def _failing_write(*args, **kwargs):
            raise OSError("Simulated disk full")

        monkeypatch.setattr(generator, "write_json", _failing_write)

        with pytest.raises(IOError, match="Simulated disk full"):
            gen.generate([_make_variant()], target)

        assert not target.exists()
        tmp_files = list(tmp_path.glob(".report_*"))
        assert tmp_files == []


class TestReportGeneratorCSV:
    """Tests for CSV output format routing."""

    def test_generates_csv_with_variants(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)
        variants = [_make_variant(), _make_variant(pos=99999)]
        output = tmp_path / "report.csv"

        result = gen.generate(variants, output)

        assert result == output
        assert output.exists()
        lines = output.read_text(encoding="utf-8").strip().split("\n")
        assert len(lines) == 3  # header + 2 data rows
        assert lines[0].startswith("chromosome,")

    def test_generates_csv_empty_list(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)
        output = tmp_path / "empty.csv"

        result = gen.generate([], output)

        assert result == output
        lines = output.read_text(encoding="utf-8").strip().split("\n")
        assert len(lines) == 1  # header only


class TestReportGeneratorPDF:
    """Tests for PDF output format routing."""

    def test_generates_pdf_when_reportlab_available(self, tmp_path: Path) -> None:
        pytest.importorskip("reportlab")
        config = ReportConfig(output_format="pdf")
        gen = ReportGenerator(config)
        variants = [_make_variant()]
        output = tmp_path / "report.pdf"

        result = gen.generate(variants, output)

        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0

    def test_generates_pdf_empty_list(self, tmp_path: Path) -> None:
        pytest.importorskip("reportlab")
        config = ReportConfig(output_format="pdf")
        gen = ReportGenerator(config)
        output = tmp_path / "empty.pdf"

        result = gen.generate([], output)

        assert result == output
        assert output.exists()


class TestReportGeneratorAtomicity:
    """Tests for atomic write behavior."""

    def test_atomic_overwrite_existing(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        output = tmp_path / "report.json"
        output.write_text("old content")

        gen.generate([_make_variant()], output)

        data = json.loads(output.read_text(encoding="utf-8"))
        assert len(data) == 1
        assert data[0]["chromosome"] == "chr1"

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)
        output = tmp_path / "deep" / "nested" / "dir" / "report.csv"

        result = gen.generate([_make_variant()], output)

        assert result == output
        assert output.exists()

    def test_no_temp_file_left_on_success(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        output = tmp_path / "report.json"

        gen.generate([_make_variant()], output)

        tmp_files = list(tmp_path.glob(".report_*"))
        assert tmp_files == []
