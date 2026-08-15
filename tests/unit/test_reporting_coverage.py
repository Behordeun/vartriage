"""Unit tests for reporting/generator.py and reporting/csv_writer.py uncovered paths."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vartriage.models.config import ClinicalReportConfig, ReportConfig
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
from vartriage.reporting.csv_writer import (
    _format_disease_associations,
    _format_field,
    _get_constraint_field,
    _get_gene_context_field,
    _variant_to_row,
    write_csv,
)
from vartriage.reporting.generator import ReportGenerator


def _make_classified(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    gene_name: str | None = None,
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    allele_frequency: float | None = 0.001,
    composite_rank: float | None = 0.7,
    revel_score: float | None = None,
    clinvar: ClinVarAssertion | None = None,
    classification: ACMGClassification = ACMGClassification.VUS,
    evidence_tags: frozenset[EvidenceTag] | None = None,
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
        gene_name=gene_name,
    )
    scored = ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
        revel_score=revel_score,
    )
    return ClassifiedVariant(
        scored=scored,
        evidence_tags=evidence_tags or frozenset(),
        classification=classification,
    )


class TestFormatField:
    """_format_field converts values to CSV-safe strings."""

    def test_none_becomes_empty(self) -> None:
        assert _format_field(None) == ""

    def test_int_converted(self) -> None:
        assert _format_field(42) == "42"

    def test_float_converted(self) -> None:
        assert _format_field(0.123) == "0.123"

    def test_string_passthrough(self) -> None:
        assert _format_field("hello") == "hello"

    def test_bool_converted(self) -> None:
        assert _format_field(True) == "True"


class TestFormatDiseaseAssociations:
    """_format_disease_associations handles gene_context presence and absence."""

    def test_returns_none_when_no_gene_context(self) -> None:
        obj = MagicMock(spec=[])
        del obj.gene_context
        assert _format_disease_associations(obj) is None

    def test_returns_none_when_gene_context_is_none(self) -> None:
        obj = MagicMock()
        obj.gene_context = None
        assert _format_disease_associations(obj) is None

    def test_returns_none_when_no_disease_associations(self) -> None:
        ctx = MagicMock()
        ctx.disease_associations = []
        obj = MagicMock()
        obj.gene_context = ctx
        assert _format_disease_associations(obj) is None

    def test_formats_disease_with_mim_and_inheritance(self) -> None:
        assoc = MagicMock()
        assoc.disease_name = "Breast cancer"
        assoc.mim_number = "113705"
        assoc.inheritance_mode = "AD"

        ctx = MagicMock()
        ctx.disease_associations = [assoc]

        obj = MagicMock()
        obj.gene_context = ctx

        result = _format_disease_associations(obj)
        assert "Breast cancer" in result
        assert "[MIM:113705]" in result
        assert "(AD)" in result

    def test_formats_disease_without_mim_and_mode(self) -> None:
        assoc = MagicMock()
        assoc.disease_name = "Some disease"
        assoc.mim_number = None
        assoc.inheritance_mode = None

        ctx = MagicMock()
        ctx.disease_associations = [assoc]

        obj = MagicMock()
        obj.gene_context = ctx

        result = _format_disease_associations(obj)
        assert result == "Some disease"

    def test_multiple_associations_semicolon_separated(self) -> None:
        assoc1 = MagicMock()
        assoc1.disease_name = "Disease A"
        assoc1.mim_number = None
        assoc1.inheritance_mode = None

        assoc2 = MagicMock()
        assoc2.disease_name = "Disease B"
        assoc2.mim_number = "600001"
        assoc2.inheritance_mode = "AR"

        ctx = MagicMock()
        ctx.disease_associations = [assoc1, assoc2]

        obj = MagicMock()
        obj.gene_context = ctx

        result = _format_disease_associations(obj)
        assert ";" in result
        assert "Disease A" in result
        assert "Disease B [MIM:600001] (AR)" in result


class TestGetGeneContextField:
    """_get_gene_context_field handles missing context gracefully."""

    def test_returns_none_when_no_gene_context_attr(self) -> None:
        obj = MagicMock(spec=[])
        assert _get_gene_context_field(obj, "clingen_validity") is None

    def test_returns_none_when_gene_context_is_none(self) -> None:
        obj = MagicMock()
        obj.gene_context = None
        assert _get_gene_context_field(obj, "clingen_validity") is None

    def test_returns_field_value(self) -> None:
        ctx = MagicMock()
        ctx.clingen_validity = "Definitive"
        obj = MagicMock()
        obj.gene_context = ctx
        assert _get_gene_context_field(obj, "clingen_validity") == "Definitive"


class TestGetConstraintField:
    """_get_constraint_field handles nested constraint objects."""

    def test_returns_none_when_no_gene_context(self) -> None:
        obj = MagicMock(spec=[])
        assert _get_constraint_field(obj, "pli") is None

    def test_returns_none_when_constraint_is_none(self) -> None:
        ctx = MagicMock()
        ctx.constraint = None
        obj = MagicMock()
        obj.gene_context = ctx
        assert _get_constraint_field(obj, "pli") is None

    def test_returns_constraint_value(self) -> None:
        constraint = MagicMock()
        constraint.pli = 0.99
        ctx = MagicMock()
        ctx.constraint = constraint
        obj = MagicMock()
        obj.gene_context = ctx
        assert _get_constraint_field(obj, "pli") == 0.99


class TestVariantToRow:
    """_variant_to_row extracts all fields in canonical order."""

    def test_all_none_optional_fields(self) -> None:
        cv = _make_classified(
            allele_frequency=None,
            composite_rank=None,
            revel_score=None,
            clinvar=None,
            evidence_tags=frozenset(),
        )
        row = _variant_to_row(cv)
        assert len(row) == 21
        # AF, revel, composite, clinvar, evidence_tags all empty
        assert row[6] == ""
        assert row[7] == ""
        assert row[8] == ""
        assert row[9] == ""
        assert row[11] == ""

    def test_evidence_tags_sorted_by_value(self) -> None:
        tags = frozenset({EvidenceTag.PP3, EvidenceTag.PVS1, EvidenceTag.PM2})
        cv = _make_classified(evidence_tags=tags)
        row = _variant_to_row(cv)
        # Tags should be semicolon-joined and sorted alphabetically
        assert row[11] == "PM2;PP3;PVS1"


class TestWriteCsvIterator:
    """write_csv handles iterator inputs (not just sequences)."""

    def test_accepts_iterator(self, tmp_path: Path) -> None:
        output = tmp_path / "output.csv"

        def gen_variants():
            yield _make_classified(pos=1)
            yield _make_classified(pos=2)

        result = write_csv(gen_variants(), output)
        assert result == output
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 rows

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        output = tmp_path / "deep" / "nested" / "report.csv"
        write_csv([], output)
        assert output.exists()


class TestReportGeneratorUnsupportedFormat:
    """ReportGenerator raises on unknown format strings."""

    def test_raises_on_unknown_format(self, tmp_path: Path) -> None:
        config = ReportConfig.__new__(ReportConfig)
        object.__setattr__(config, "output_format", "xlsx")

        gen = ReportGenerator(config)
        output = tmp_path / "report.xlsx"

        with pytest.raises(IOError, match="Unsupported output format"):
            gen.generate([], output)


class TestReportGeneratorVcfFormat:
    """ReportGenerator VCF path requires source_vcf_path."""

    def test_raises_without_source_vcf_path(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="vcf")
        gen = ReportGenerator(config)
        output = tmp_path / "output.vcf"

        with pytest.raises(IOError, match="source_vcf_path"):
            gen.generate([], output)

    def test_vcf_format_delegates_to_write_vcf(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="vcf")
        gen = ReportGenerator(config)
        output = tmp_path / "output.vcf"
        source = tmp_path / "input.vcf"
        source.write_text("")

        with patch("vartriage.reporting.vcf_writer.write_vcf") as mock_write:
            result = gen.generate([], output, source_vcf_path=source)

        mock_write.assert_called_once()
        assert result == output


class TestReportGeneratorClinicalFormat:
    """ReportGenerator clinical format delegation."""

    def test_raises_without_clinical_config(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="clinical-pdf")
        gen = ReportGenerator(config, clinical_config=None)
        output = tmp_path / "report.pdf"

        with pytest.raises(IOError, match="ClinicalReportConfig"):
            gen.generate([], output)

    def test_delegates_to_clinical_generator(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="clinical-html")
        clinical_config = ClinicalReportConfig(
            patient_id="P001",
            panel_name="Panel A",
            output_format="clinical-html",
        )
        gen = ReportGenerator(
            config,
            clinical_config=clinical_config,
            reference_checksums={"ref.gtf": "abc123"},
        )
        output = tmp_path / "report.html"

        with patch(
            "vartriage.reporting.clinical.generator.ClinicalReportGenerator"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_instance.generate.return_value = output
            mock_cls.return_value = mock_instance

            result = gen.generate([], output)

        mock_cls.assert_called_once()
        call_kwargs = mock_cls.call_args.kwargs
        assert call_kwargs["config"] == clinical_config
        assert call_kwargs["reference_checksums"] == {"ref.gtf": "abc123"}
        assert result == output


class TestReportGeneratorPdfFallback:
    """PDF generation falls back when reportlab is absent."""

    def test_uses_fallback_when_reportlab_missing(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="pdf")
        gen = ReportGenerator(config)
        output = tmp_path / "report.pdf"
        variants = [_make_classified()]

        with patch(
            "vartriage.reporting.generator.ReportGenerator._write_pdf"
        ) as mock_pdf:
            mock_pdf.return_value = output
            # Trigger through the main generate method
            tmp_result = gen.generate(variants, output)

        assert tmp_result == output


class TestReportGeneratorErrorWrapping:
    """Non-IOError exceptions get wrapped in IOError."""

    def test_wraps_unexpected_exceptions(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        output = tmp_path / "report.json"

        with (
            patch(
                "vartriage.reporting.generator.write_json",
                side_effect=RuntimeError("unexpected"),
            ),
            pytest.raises(IOError, match="Failed to generate JSON"),
        ):
            gen.generate([_make_classified()], output)

        # Temp file should be cleaned up
        tmp_files = list(tmp_path.glob(".report_*"))
        assert tmp_files == []


class TestReportGeneratorTempCleanup:
    """Temp file cleanup on various failure modes."""

    def test_csv_failure_cleans_temp(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)
        output = tmp_path / "report.csv"

        with (
            patch(
                "vartriage.reporting.generator.write_csv",
                side_effect=OSError("disk full"),
            ),
            pytest.raises(IOError),
        ):
            gen.generate([_make_classified()], output)

        assert not output.exists()
        tmp_files = list(tmp_path.glob(".report_*"))
        assert tmp_files == []
