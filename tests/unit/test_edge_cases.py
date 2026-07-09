"""Unit tests for edge cases and error paths.

Covers boundary values, empty input handling, error message content,
and file I/O error simulation across all pipeline stages.
"""

from __future__ import annotations

import os
import tempfile
import warnings
from pathlib import Path
from typing import Iterator
from unittest.mock import patch

import pytest

from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.exceptions import ParseError
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.config import (
    AnnotationConfig,
    PrioritizationConfig,
    QualityFilterConfig,
    ReportConfig,
)
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
from vartriage.prioritization.frequency_filter import FrequencyFilter
from vartriage.prioritization.scoring import (
    normalize_cadd_scores,
    score_variants,
    validate_revel_scores,
)
from vartriage.classification.acmg import ACMGClassifier
from vartriage.reporting.generator import ReportGenerator


def _make_variant(
    chrom: str = "chr1",
    pos: int = 100,
    qual: float | None = 30.0,
    filter_status: str = "PASS",
) -> Variant:
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=qual,
        filter_status=filter_status,
    )


def _make_annotated(
    chrom: str = "chr1",
    pos: int = 100,
    allele_frequency: float | None = 0.005,
    frequency_unknown: bool = False,
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
) -> AnnotatedVariant:
    variant = _make_variant(chrom=chrom, pos=pos)
    return AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=allele_frequency,
        frequency_unknown=frequency_unknown,
    )


def _make_scored(
    chrom: str = "chr1",
    pos: int = 100,
    composite_rank: float | None = 0.5,
    revel_score: float | None = None,
    allele_frequency: float | None = 0.001,
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    clinvar_assertion: ClinVarAssertion | None = None,
) -> ScoredVariant:
    annotated = AnnotatedVariant(
        variant=_make_variant(chrom=chrom, pos=pos),
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
    )
    return ScoredVariant(
        annotated=annotated,
        composite_rank=composite_rank,
        revel_score=revel_score,
    )


def _make_classified(
    chrom: str = "chr1",
    pos: int = 100,
) -> ClassifiedVariant:
    scored = _make_scored(chrom=chrom, pos=pos)
    return ClassifiedVariant(
        scored=scored,
        evidence_tags=frozenset(),
        classification=ACMGClassification.VUS,
    )


# --------------------------------------------------------------------------
# Boundary Value Tests
# --------------------------------------------------------------------------


class TestQualBoundaryValues:
    """QUAL threshold boundary testing."""

    def test_qual_exactly_at_default_threshold_20(self) -> None:
        """QUAL == 20.0 passes with default config (threshold is 20)."""
        qf = QualityFilter()
        variant = _make_variant(qual=20.0, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == [variant]

    def test_qual_just_below_default_threshold(self) -> None:
        """QUAL = 19.999... is excluded at default threshold of 20."""
        qf = QualityFilter()
        variant = _make_variant(qual=19.999, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == []

    def test_qual_exactly_zero_passes_with_zero_threshold(self) -> None:
        """QUAL == 0 passes when threshold is 0."""
        qf = QualityFilter(QualityFilterConfig(min_qual=0.0))
        variant = _make_variant(qual=0.0, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == [variant]

    def test_qual_at_max_threshold_boundary(self) -> None:
        """QUAL == 1_000_000 passes at max threshold."""
        qf = QualityFilter(QualityFilterConfig(min_qual=1_000_000))
        variant = _make_variant(qual=1_000_000, filter_status="PASS")
        result = list(qf.apply(iter([variant])))
        assert result == [variant]


class TestAFBoundaryValues:
    """Allele frequency threshold boundary testing."""

    def test_af_exactly_at_default_threshold_0_01(self) -> None:
        """AF == 0.01 is retained (threshold excludes strictly greater)."""
        ff = FrequencyFilter()
        variant = _make_annotated(allele_frequency=0.01)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_af_barely_above_default_threshold(self) -> None:
        """AF == 0.0100001 is excluded."""
        ff = FrequencyFilter()
        variant = _make_annotated(allele_frequency=0.0100001)
        result = list(ff.apply(iter([variant])))
        assert result == []

    def test_af_exactly_zero_retained(self) -> None:
        """AF == 0.0 is always retained."""
        ff = FrequencyFilter()
        variant = _make_annotated(allele_frequency=0.0)
        result = list(ff.apply(iter([variant])))
        assert result == [variant]

    def test_af_threshold_at_zero_excludes_everything_above(self) -> None:
        """With max_af=0.0, only variants at exactly 0.0 pass."""
        ff = FrequencyFilter(PrioritizationConfig(max_allele_frequency=0.0))
        at_zero = _make_annotated(allele_frequency=0.0)
        above_zero = _make_annotated(allele_frequency=0.0001)
        result = list(ff.apply(iter([at_zero, above_zero])))
        assert result == [at_zero]


class TestRevelBoundaryValues:
    """REVEL score boundary testing (valid range 0.0 to 1.0)."""

    def test_revel_exactly_0_is_valid(self) -> None:
        """REVEL == 0.0 is accepted."""
        result = validate_revel_scores([0.0])
        assert result == [0.0]

    def test_revel_exactly_1_is_valid(self) -> None:
        """REVEL == 1.0 is accepted."""
        result = validate_revel_scores([1.0])
        assert result == [1.0]

    def test_revel_slightly_above_1_rejected(self) -> None:
        """REVEL == 1.0001 triggers a validation warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_revel_scores([1.0001])
            assert result == [None]
            assert len(w) == 1

    def test_revel_slightly_below_0_rejected(self) -> None:
        """REVEL == -0.0001 triggers a validation warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_revel_scores([-0.0001])
            assert result == [None]
            assert len(w) == 1


class TestCaddBoundaryValues:
    """CADD Phred score boundary testing."""

    def test_cadd_exactly_0_normalizes_to_0(self) -> None:
        """CADD == 0.0 normalizes to 0.0."""
        result = normalize_cadd_scores([0.0])
        assert result == [0.0]

    def test_cadd_exactly_99_normalizes_to_1(self) -> None:
        """CADD == 99.0 normalizes to 1.0."""
        result = normalize_cadd_scores([99.0])
        assert result == [pytest.approx(1.0)]

    def test_cadd_above_99_capped_at_1(self) -> None:
        """CADD == 120.0 caps at 1.0."""
        result = normalize_cadd_scores([120.0])
        assert result == [1.0]

    def test_cadd_slightly_negative_rejected(self) -> None:
        """CADD == -0.001 triggers warning and returns None."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = normalize_cadd_scores([-0.001])
            assert result == [None]
            assert len(w) == 1


# --------------------------------------------------------------------------
# Empty Input Handling
# --------------------------------------------------------------------------


class TestEmptyInputHandling:
    """All pipeline stages produce empty output without error on empty input."""

    def test_quality_filter_empty_input(self) -> None:
        qf = QualityFilter()
        result = list(qf.apply(iter([])))
        assert result == []

    def test_frequency_filter_empty_input(self) -> None:
        ff = FrequencyFilter()
        result = list(ff.apply(iter([])))
        assert result == []

    def test_scoring_empty_input(self) -> None:
        result = score_variants([], [], [])
        assert result == []

    def test_acmg_classifier_empty_input(self) -> None:
        classifier = ACMGClassifier()
        result = list(classifier.classify(iter([])))
        assert result == []

    def test_annotation_engine_empty_input(self, tmp_path: Path) -> None:
        """AnnotationEngine.annotate yields nothing on empty input."""
        from unittest.mock import MagicMock, patch as mock_patch
        from vartriage.annotation.engine import AnnotationEngine

        with mock_patch.object(
            AnnotationEngine, "__init__", lambda self, config: None
        ):
            engine = AnnotationEngine.__new__(AnnotationEngine)
            engine._config = AnnotationConfig(
                gene_annotation_path=tmp_path / "fake.gtf",
                gnomad_path=tmp_path / "fake.vcf",
                batch_size=1000,
            )
            engine._warnings = []
            engine._consequence_annotator = MagicMock()
            engine._frequency_db = MagicMock()
            engine._clinvar_db = None

            result = list(engine.annotate(iter([])))
            assert result == []

    def test_report_generator_empty_json(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        output = tmp_path / "empty.json"
        result = gen.generate([], output)
        assert result == output
        assert output.exists()
        import json
        data = json.loads(output.read_text())
        assert data == []

    def test_report_generator_empty_csv(self, tmp_path: Path) -> None:
        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)
        output = tmp_path / "empty.csv"
        result = gen.generate([], output)
        assert result == output
        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 1  # header only

    def test_report_generator_empty_pdf(self, tmp_path: Path) -> None:
        """Empty variant list produces a valid PDF (structure only)."""
        config = ReportConfig(output_format="pdf")
        gen = ReportGenerator(config)
        output = tmp_path / "empty.pdf"
        result = gen.generate([], output)
        assert result == output
        assert output.exists()
        assert output.stat().st_size > 0


# --------------------------------------------------------------------------
# Error Message Content Tests
# --------------------------------------------------------------------------


class TestParseErrorContent:
    """ParseError includes line_number and detail."""

    def test_parse_error_has_line_number_attribute(self) -> None:
        err = ParseError(line_number=42, detail="Something went wrong")
        assert err.line_number == 42

    def test_parse_error_has_detail_attribute(self) -> None:
        err = ParseError(
            line_number=7, detail="Missing mandatory QUAL column"
        )
        assert err.detail == "Missing mandatory QUAL column"

    def test_parse_error_has_field_attribute_when_set(self) -> None:
        err = ParseError(
            line_number=10,
            field="POS",
            detail="Non-integer value 'abc'",
        )
        assert err.field == "POS"

    def test_parse_error_field_none_when_unset(self) -> None:
        err = ParseError(line_number=1, detail="Missing header")
        assert err.field is None

    def test_parse_error_str_includes_line_number(self) -> None:
        err = ParseError(line_number=5, detail="Bad data")
        assert "5" in str(err)

    def test_parse_error_str_includes_detail(self) -> None:
        err = ParseError(line_number=5, detail="Bad data")
        assert "Bad data" in str(err)

    def test_parse_error_str_includes_field_when_present(self) -> None:
        err = ParseError(
            line_number=12, field="QUAL", detail="Non-numeric"
        )
        assert "QUAL" in str(err)

    def test_missing_fileformat_raises_with_appropriate_detail(
        self, tmp_path: Path
    ) -> None:
        """VCFParser raises ParseError mentioning fileformat on missing declaration."""
        bad_content = (
            "##INFO=<ID=DP,Number=1,Type=Integer,"
            'Description="Depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n"
        )
        vcf_path = tmp_path / "bad.vcf"
        vcf_path.write_text(bad_content)
        with pytest.raises(ParseError) as exc_info:
            VCFParser(vcf_path)
        assert exc_info.value.line_number >= 1
        assert "fileformat" in exc_info.value.detail.lower()


class TestValueErrorContent:
    """ValueError from configs specifies the valid range."""

    def test_quality_filter_config_error_specifies_range(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            QualityFilterConfig(min_qual=-1.0)
        msg = str(exc_info.value)
        assert "0" in msg
        assert "1000000" in msg

    def test_quality_filter_config_error_includes_invalid_value(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            QualityFilterConfig(min_qual=2_000_000)
        msg = str(exc_info.value)
        assert "2000000" in msg

    def test_prioritization_config_af_error_specifies_range(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            PrioritizationConfig(max_allele_frequency=-0.5)
        msg = str(exc_info.value)
        assert "0.0" in msg
        assert "1.0" in msg

    def test_prioritization_config_af_error_includes_value(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            PrioritizationConfig(max_allele_frequency=2.0)
        msg = str(exc_info.value)
        assert "2.0" in msg

    def test_annotation_config_batch_size_error_specifies_range(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            AnnotationConfig(
                gene_annotation_path=Path("/fake"),
                gnomad_path=Path("/fake"),
                batch_size=500,
            )
        msg = str(exc_info.value)
        assert "1000" in msg
        assert "100000" in msg

    def test_prioritization_config_batch_size_error(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            PrioritizationConfig(batch_size=200_000)
        msg = str(exc_info.value)
        assert "1000" in msg
        assert "100000" in msg


class TestFileNotFoundErrorContent:
    """FileNotFoundError includes the file path."""

    def test_vcf_parser_file_not_found_includes_path(
        self, tmp_path: Path
    ) -> None:
        fake_path = tmp_path / "nonexistent.vcf"
        with pytest.raises(FileNotFoundError) as exc_info:
            VCFParser(fake_path)
        assert str(fake_path) in str(exc_info.value)

    def test_vcf_parser_missing_tbi_includes_index_path(
        self, tmp_path: Path
    ) -> None:
        gz_path = tmp_path / "test.vcf.gz"
        gz_path.write_bytes(b"fake")
        expected_tbi = str(gz_path) + ".tbi"
        with pytest.raises(FileNotFoundError) as exc_info:
            VCFParser(gz_path)
        assert expected_tbi in str(exc_info.value)


# --------------------------------------------------------------------------
# ReportGenerator I/O Error Handling
# --------------------------------------------------------------------------


class TestReportGeneratorIOErrors:
    """ReportGenerator raises IOError on write failure without partial output."""

    def test_raises_ioerror_on_unwritable_directory(
        self, tmp_path: Path
    ) -> None:
        """When parent directory is not writable, should raise IOError."""
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        variants = [_make_classified()]

        readonly_dir = tmp_path / "readonly"
        readonly_dir.mkdir()
        output = readonly_dir / "report.json"

        os.chmod(readonly_dir, 0o444)
        try:
            with pytest.raises((IOError, OSError, PermissionError)):
                gen.generate(variants, output)
        finally:
            os.chmod(readonly_dir, 0o755)
        # After restoring perms, verify no output was written
        assert not output.exists()

    def test_no_partial_output_on_write_error(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """On write failure, no partial file remains at target path."""
        from vartriage.reporting import generator

        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        target = tmp_path / "output.json"

        def _failing_write(*args, **kwargs):
            raise OSError("Disk full simulation")

        monkeypatch.setattr(generator, "write_json", _failing_write)

        with pytest.raises(IOError, match="Disk full simulation"):
            gen.generate([_make_classified()], target)

        assert not target.exists()

    def test_csv_write_error_no_partial(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """CSV write failures also leave no partial file."""
        from vartriage.reporting import generator

        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)
        target = tmp_path / "output.csv"

        def _failing_write(*args, **kwargs):
            raise OSError("Permission denied")

        monkeypatch.setattr(generator, "write_csv", _failing_write)

        with pytest.raises(IOError, match="Permission denied"):
            gen.generate([_make_classified()], target)

        assert not target.exists()

    def test_ioerror_message_specifies_failure_reason(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The IOError message includes the specific failure reason."""
        from vartriage.reporting import generator

        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)
        target = tmp_path / "output.json"

        def _failing_write(*args, **kwargs):
            raise RuntimeError("Encoding failure on non-UTF8 char")

        monkeypatch.setattr(generator, "write_json", _failing_write)

        with pytest.raises(IOError) as exc_info:
            gen.generate([_make_classified()], target)

        assert "Encoding failure" in str(exc_info.value)
