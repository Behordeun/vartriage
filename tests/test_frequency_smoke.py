"""Smoke test for frequency lookup implementations."""

import tempfile
from pathlib import Path

from vartriage.annotation.frequency import DictFrequencyDatabase
from vartriage.io.exceptions import ReferenceFileError
from vartriage.models.warnings import MissingDataWarning


def create_test_gnomad_file() -> Path:
    """Create a temporary gnomAD reference TSV file."""
    content = (
        "chrom\tpos\tref\talt\taf\n"
        "chr1\t100\tA\tT\t0.001\n"
        "chr1\t200\tG\tC\t0.05\n"
        "chr2\t500\tAT\tA\t0.0001\n"
        "chrX\t1000\tC\tG\t0.25\n"
    )
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    tmp.write(content)
    tmp.close()
    return Path(tmp.name)


def test_dict_frequency_database():
    """Test DictFrequencyDatabase with a valid reference file."""
    ref_path = create_test_gnomad_file()

    db = DictFrequencyDatabase()
    db.load(ref_path)

    variants = [
        ("chr1", 100, "A", "T"),
        ("chr1", 200, "G", "C"),
        ("chr1", 300, "A", "G"),
        ("chr2", 500, "AT", "A"),
        ("chr3", 999, "T", "C"),
    ]

    results = db.lookup_batch(variants)

    assert results[0] == 0.001
    assert results[1] == 0.05
    assert results[2] is None
    assert results[3] == 0.0001
    assert results[4] is None

    assert len(db.warnings) == 2
    assert all(isinstance(w, MissingDataWarning) for w in db.warnings)
    assert db.warnings[0].chrom == "chr1"
    assert db.warnings[0].pos == 300
    assert db.warnings[0].source == "gnomAD"
    assert db.warnings[1].chrom == "chr3"
    assert db.warnings[1].pos == 999


def test_dict_frequency_missing_file():
    """Test that ReferenceFileError is raised for missing files."""
    db = DictFrequencyDatabase()
    try:
        db.load(Path("/nonexistent/gnomad.tsv"))
        assert False, "Should have raised ReferenceFileError"
    except ReferenceFileError as exc:
        assert "file not found" in str(exc)


def test_dict_frequency_bad_format():
    """Test ReferenceFileError for files with wrong columns."""
    content = "col_a\tcol_b\n1\t2\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    tmp.write(content)
    tmp.close()

    db = DictFrequencyDatabase()
    try:
        db.load(Path(tmp.name))
        assert False, "Should have raised ReferenceFileError"
    except ReferenceFileError as exc:
        assert "missing required columns" in str(exc)


def test_dict_frequency_empty_batch():
    """Test that empty batch returns empty list."""
    ref_path = create_test_gnomad_file()
    db = DictFrequencyDatabase()
    db.load(ref_path)

    results = db.lookup_batch([])
    assert results == []


def test_polars_frequency_database():
    """Test PolarsFrequencyDatabase if polars is available."""
    try:
        from vartriage.annotation.frequency_polars import (
            POLARS_AVAILABLE, PolarsFrequencyDatabase)
    except ImportError:
        return

    if not POLARS_AVAILABLE:
        return

    ref_path = create_test_gnomad_file()
    db = PolarsFrequencyDatabase()
    db.load(ref_path)

    variants = [
        ("chr1", 100, "A", "T"),
        ("chr1", 200, "G", "C"),
        ("chr1", 300, "A", "G"),
        ("chr2", 500, "AT", "A"),
        ("chr3", 999, "T", "C"),
    ]

    results = db.lookup_batch(variants)

    assert results[0] == 0.001
    assert results[1] == 0.05
    assert results[2] is None
    assert results[3] == 0.0001
    assert results[4] is None

    assert len(db.warnings) == 2
    assert db.warnings[0].chrom == "chr1"
    assert db.warnings[0].pos == 300
    assert db.warnings[1].chrom == "chr3"
    assert db.warnings[1].pos == 999


def test_polars_frequency_missing_file():
    """Test PolarsFrequencyDatabase with missing file."""
    try:
        from vartriage.annotation.frequency_polars import (
            POLARS_AVAILABLE, PolarsFrequencyDatabase)
    except ImportError:
        return

    if not POLARS_AVAILABLE:
        return

    db = PolarsFrequencyDatabase()
    try:
        db.load(Path("/nonexistent/gnomad.tsv"))
        assert False, "Should have raised ReferenceFileError"
    except ReferenceFileError as exc:
        assert "file not found" in str(exc)


def test_polars_frequency_empty_batch():
    """Test PolarsFrequencyDatabase with empty batch."""
    try:
        from vartriage.annotation.frequency_polars import (
            POLARS_AVAILABLE, PolarsFrequencyDatabase)
    except ImportError:
        return

    if not POLARS_AVAILABLE:
        return

    ref_path = create_test_gnomad_file()
    db = PolarsFrequencyDatabase()
    db.load(ref_path)

    results = db.lookup_batch([])
    assert results == []
