"""Tests for PyPI release readiness: streaming, ScoreLoader, CLI, warnings."""

from __future__ import annotations

import logging
import tempfile
import warnings
from pathlib import Path
from typing import Optional
from unittest.mock import patch

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from vartriage.models.config import ReportConfig
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
from vartriage.reporting.generator import ReportGenerator
from vartriage.prioritization.score_loader import CoordinateKey, ScoreLoader
from vartriage.exceptions import VarTriageWarning
from vartriage.prioritization.scoring import ScoreValidationWarning
from vartriage._internal.warning_accumulator import MissingDataSummaryWarning

from tests.generators.variants import scored_variant, evidence_tag_set


# ---------------------------------------------------------------------------
# Shared strategies
# ---------------------------------------------------------------------------

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


@st.composite
def st_classified_variant(draw: st.DrawFn) -> ClassifiedVariant:
    """Build a random ClassifiedVariant."""
    sv = draw(scored_variant())
    tags = draw(evidence_tag_set())
    classification = draw(st.sampled_from(list(ACMGClassification)))
    missing_sources = draw(
        st.frozensets(
            st.sampled_from(["gnomAD", "REVEL", "ClinVar", "CADD"]),
            min_size=0,
            max_size=3,
        )
    )
    return ClassifiedVariant(
        scored=sv,
        evidence_tags=tags,
        classification=classification,
        missing_data_sources=missing_sources,
    )


@st.composite
def st_coordinate_key(draw: st.DrawFn) -> CoordinateKey:
    """Build a random (chrom, pos, ref, alt) tuple."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    pos = draw(st.integers(min_value=1, max_value=250_000_000))
    ref = draw(
        st.text(
            alphabet="ACGT", min_size=1, max_size=5
        )
    )
    alt = draw(
        st.text(
            alphabet="ACGT", min_size=1, max_size=5
        )
    )
    return (chrom, pos, ref, alt)


# ===========================================================================
# Task 9.1: Property test for iterator/sequence output equivalence
# ===========================================================================


class TestIteratorSequenceEquivalence:
    """Iterator vs sequence output should be byte-identical."""

    @given(
        variants=st.lists(
            st_classified_variant(), min_size=0, max_size=30
        )
    )
    @settings(max_examples=100, deadline=None)
    def test_json_output_identical_for_iter_and_sequence(
        self, variants: list[ClassifiedVariant]
    ) -> None:
        """JSON bytes match regardless of iter() vs list input."""
        config = ReportConfig(output_format="json")
        gen = ReportGenerator(config)

        with tempfile.TemporaryDirectory() as tmpdir:
            path_seq = Path(tmpdir) / "seq.json"
            path_iter = Path(tmpdir) / "iter.json"

            gen.generate(variants, path_seq)
            gen.generate(iter(variants), path_iter)

            bytes_seq = path_seq.read_bytes()
            bytes_iter = path_iter.read_bytes()

        assert bytes_seq == bytes_iter

    @given(
        variants=st.lists(
            st_classified_variant(), min_size=0, max_size=30
        )
    )
    @settings(max_examples=100, deadline=None)
    def test_csv_output_identical_for_iter_and_sequence(
        self, variants: list[ClassifiedVariant]
    ) -> None:
        """CSV bytes match regardless of iter() vs list input."""
        config = ReportConfig(output_format="csv")
        gen = ReportGenerator(config)

        with tempfile.TemporaryDirectory() as tmpdir:
            path_seq = Path(tmpdir) / "seq.csv"
            path_iter = Path(tmpdir) / "iter.csv"

            gen.generate(variants, path_seq)
            gen.generate(iter(variants), path_iter)

            bytes_seq = path_seq.read_bytes()
            bytes_iter = path_iter.read_bytes()

        assert bytes_seq == bytes_iter


# ===========================================================================
# Task 9.4: Property test for ScoreLoader round-trip correctness
# ===========================================================================


class TestScoreLoaderRoundTrip:
    """Round-trip and lookup correctness for ScoreLoader."""

    @given(
        data=st.dictionaries(
            keys=st_coordinate_key(),
            values=st.floats(
                min_value=0.0, max_value=99.0,
                allow_nan=False, allow_infinity=False,
            ),
            min_size=0,
            max_size=30,
        )
    )
    @settings(max_examples=100, deadline=None)
    def test_load_cadd_round_trip(
        self, data: dict[CoordinateKey, float]
    ) -> None:
        """Write a TSV then load_cadd. Each key maps to its original score."""
        loader = ScoreLoader()

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "scores.tsv"
            self._write_tsv(tsv_path, data)

            loaded = loader.load_cadd(tsv_path)

        for key, expected_score in data.items():
            assert key in loaded, f"Key {key} missing from loaded data"
            assert loaded[key] == pytest.approx(expected_score, rel=1e-9), (
                f"Score mismatch for {key}: "
                f"expected {expected_score}, got {loaded[key]}"
            )

    @given(
        data=st.dictionaries(
            keys=st_coordinate_key(),
            values=st.floats(
                min_value=0.0, max_value=99.0,
                allow_nan=False, allow_infinity=False,
            ),
            min_size=1,
            max_size=20,
        ),
        missing_keys=st.lists(
            st_coordinate_key(), min_size=1, max_size=5
        ),
    )
    @settings(max_examples=100, deadline=None)
    def test_lookup_batch_returns_none_for_missing(
        self,
        data: dict[CoordinateKey, float],
        missing_keys: list[CoordinateKey],
    ) -> None:
        """Keys not in the loaded dict come back as None."""
        loader = ScoreLoader()

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "scores.tsv"
            self._write_tsv(tsv_path, data)
            loaded = loader.load_cadd(tsv_path)

        truly_missing = [k for k in missing_keys if k not in data]
        if not truly_missing:
            return

        results = loader.lookup_batch(truly_missing, loaded)
        for i, result in enumerate(results):
            assert result is None, (
                f"Expected None for missing key {truly_missing[i]}, "
                f"got {result}"
            )

    @given(
        data=st.dictionaries(
            keys=st_coordinate_key(),
            values=st.floats(
                min_value=0.0, max_value=1.0,
                allow_nan=False, allow_infinity=False,
            ),
            min_size=1,
            max_size=20,
        )
    )
    @settings(max_examples=100, deadline=None)
    def test_load_revel_round_trip(
        self, data: dict[CoordinateKey, float]
    ) -> None:
        """load_revel behaves identically to load_cadd for round-trips."""
        loader = ScoreLoader()

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "revel.tsv"
            self._write_tsv(tsv_path, data)
            loaded = loader.load_revel(tsv_path)

        keys = list(data.keys())
        results = loader.lookup_batch(keys, loaded)

        for i, key in enumerate(keys):
            assert results[i] == pytest.approx(data[key], rel=1e-9)

    @staticmethod
    def _write_tsv(path: Path, data: dict[CoordinateKey, float]) -> None:
        """Helper to write coordinate-score data as a TSV file."""
        with open(path, "w", encoding="utf-8") as f:
            f.write("#chrom\tpos\tref\talt\tscore\n")
            for (chrom, pos, ref, alt), score in data.items():
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{score}\n")


# ===========================================================================
# Task 9.6: Unit tests for CLI argument parsing and error handling
# ===========================================================================


class TestCLI:
    """Unit tests for CLI argument parsing and error handling."""

    def test_vcf_and_output_produces_exit_0(
        self, tmp_path: Path
    ) -> None:
        """Valid --vcf + --output exits 0 with a mocked pipeline."""
        from vartriage.cli import main

        vcf_file = tmp_path / "input.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")
        output_file = tmp_path / "report.json"

        with patch("vartriage.cli._run_pipeline") as mock_pipeline:
            mock_pipeline.return_value = output_file
            with pytest.raises(SystemExit) as exc_info:
                main(["--vcf", str(vcf_file), "--output", str(output_file)])

        assert exc_info.value.code == 0

    def test_missing_required_args_produces_nonzero_exit(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """No args at all → non-zero exit."""
        from vartriage.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main([])

        assert exc_info.value.code != 0

    def test_version_flag_prints_version_and_exits_0(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--version prints version string and exits 0."""
        from vartriage.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])

        assert exc_info.value.code == 0
        captured = capsys.readouterr()
        assert "vartriage" in captured.out

    def test_invalid_vcf_path_produces_error_on_stderr(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Non-existent VCF path → exit 1 with error on stderr."""
        from vartriage.cli import main

        nonexistent = tmp_path / "no_such_file.vcf"
        output_file = tmp_path / "report.json"

        with pytest.raises(SystemExit) as exc_info:
            main([
                "--vcf", str(nonexistent),
                "--output", str(output_file),
            ])

        assert exc_info.value.code == 1
        captured = capsys.readouterr()
        assert "Error" in captured.err
        assert "not found" in captured.err


# ===========================================================================
# Task 9.7: Unit tests for warning hierarchy
# ===========================================================================


class TestWarningHierarchy:
    """Warning class inheritance works as expected."""

    def test_score_validation_warning_inherits_vartriage_warning(self) -> None:
        """ScoreValidationWarning is a subclass of VarTriageWarning."""
        assert issubclass(ScoreValidationWarning, VarTriageWarning)

    def test_missing_data_summary_warning_inherits_vartriage_warning(
        self,
    ) -> None:
        """MissingDataSummaryWarning is a subclass of VarTriageWarning."""
        assert issubclass(MissingDataSummaryWarning, VarTriageWarning)

    def test_vartriage_warning_inherits_user_warning(self) -> None:
        """VarTriageWarning is a subclass of UserWarning."""
        assert issubclass(VarTriageWarning, UserWarning)

    def test_filter_vartriage_warning_suppresses_subclass_warnings(
        self,
    ) -> None:
        """filterwarnings with VarTriageWarning suppresses both subclasses."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            warnings.filterwarnings("ignore", category=VarTriageWarning)

            warnings.warn("score issue", ScoreValidationWarning)
            warnings.warn(
                MissingDataSummaryWarning(
                    total_count=5,
                    sources=frozenset({"gnomAD"}),
                    not_found_count=3,
                    connection_failure_count=2,
                )
            )

        vartriage_warnings = [
            w for w in caught
            if issubclass(w.category, VarTriageWarning)
        ]
        assert len(vartriage_warnings) == 0

    def test_import_vartriage_warning_from_package(self) -> None:
        """from vartriage import VarTriageWarning succeeds."""
        from vartriage import VarTriageWarning as VTW

        assert VTW is VarTriageWarning


# ===========================================================================
# Task 9.8: Unit tests for ScoreLoader
# ===========================================================================


class TestScoreLoaderUnit:
    """Unit tests for ScoreLoader edge cases and error handling."""

    def test_valid_cadd_tsv_loading_and_lookup(
        self, tmp_path: Path
    ) -> None:
        """Well-formed CADD TSV loads correctly."""
        tsv_content = (
            "#chrom\tpos\tref\talt\tscore\n"
            "chr1\t100\tA\tT\t25.3\n"
            "chr2\t200\tG\tC\t30.0\n"
            "chrX\t500\tAA\tTT\t15.7\n"
        )
        tsv_path = tmp_path / "cadd.tsv"
        tsv_path.write_text(tsv_content, encoding="utf-8")

        loader = ScoreLoader()
        scores = loader.load_cadd(tsv_path)

        assert scores[("chr1", 100, "A", "T")] == pytest.approx(25.3)
        assert scores[("chr2", 200, "G", "C")] == pytest.approx(30.0)
        assert scores[("chrX", 500, "AA", "TT")] == pytest.approx(15.7)

    def test_valid_revel_tsv_loading_and_lookup(
        self, tmp_path: Path
    ) -> None:
        """Well-formed REVEL TSV loads correctly."""
        tsv_content = (
            "#chrom\tpos\tref\talt\tscore\n"
            "chr1\t100\tA\tG\t0.85\n"
            "chr3\t999\tC\tT\t0.12\n"
        )
        tsv_path = tmp_path / "revel.tsv"
        tsv_path.write_text(tsv_content, encoding="utf-8")

        loader = ScoreLoader()
        scores = loader.load_revel(tsv_path)

        assert scores[("chr1", 100, "A", "G")] == pytest.approx(0.85)
        assert scores[("chr3", 999, "C", "T")] == pytest.approx(0.12)

    def test_missing_coordinate_returns_none(
        self, tmp_path: Path
    ) -> None:
        """Absent coordinates come back as None in batch lookups."""
        tsv_content = (
            "#chrom\tpos\tref\talt\tscore\n"
            "chr1\t100\tA\tT\t25.3\n"
        )
        tsv_path = tmp_path / "cadd.tsv"
        tsv_path.write_text(tsv_content, encoding="utf-8")

        loader = ScoreLoader()
        scores = loader.load_cadd(tsv_path)

        results = loader.lookup_batch(
            [("chr1", 100, "A", "T"), ("chr99", 1, "G", "C")],
            scores,
        )
        assert results[0] == pytest.approx(25.3)
        assert results[1] is None

    def test_nonexistent_file_raises_value_error(
        self, tmp_path: Path
    ) -> None:
        """Non-existent path raises ValueError."""
        loader = ScoreLoader()
        fake_path = tmp_path / "nonexistent.tsv"

        with pytest.raises(ValueError, match="not found"):
            loader.load_cadd(fake_path)

    def test_empty_file_returns_empty_dict(
        self, tmp_path: Path
    ) -> None:
        """Empty (or header-only) file gives an empty dict."""
        tsv_path = tmp_path / "empty.tsv"
        tsv_path.write_text("", encoding="utf-8")

        loader = ScoreLoader()
        scores = loader.load_cadd(tsv_path)

        assert scores == {}

    def test_malformed_lines_skipped_with_warning(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Bad lines are skipped and logged; good lines still load."""
        tsv_content = (
            "#chrom\tpos\tref\talt\tscore\n"
            "chr1\t100\tA\tT\t25.3\n"
            "chr2\tNOT_A_NUMBER\tG\tC\t30.0\n"
            "chr3\t300\tA\n"
            "chr4\t400\tG\tC\tnot_numeric\n"
            "chr5\t500\tA\tT\t10.0\n"
        )
        tsv_path = tmp_path / "mixed.tsv"
        tsv_path.write_text(tsv_content, encoding="utf-8")

        loader = ScoreLoader()
        with caplog.at_level(logging.WARNING):
            scores = loader.load_cadd(tsv_path)

        assert ("chr1", 100, "A", "T") in scores
        assert ("chr5", 500, "A", "T") in scores
        assert len(scores) == 2

        assert any("malformed" in r.message.lower() for r in caplog.records)
