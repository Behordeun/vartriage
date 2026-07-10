"""Unit and property tests for cache infrastructure, TabixFrequencyDatabase,
and AnnotationEngine backend selection.

Covers tasks 7.1–7.6 of the reference-loading-performance spec.
"""

from __future__ import annotations

import logging
import os
import pickle
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from vartriage._internal.cache import (
    CacheEnvelope,
    cache_path_for,
    try_load_cache,
    try_write_cache,
)


# ---------------------------------------------------------------------------
# Task 7.1: Property test for cache path computation (Property 1)
# Validates: Requirements 1.2, 5.2
# ---------------------------------------------------------------------------


# Strategy: generate realistic filesystem path components
path_segments = st.text(
    alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.",
    min_size=1,
    max_size=40,
)

file_extensions = st.sampled_from([
    ".gtf", ".gtf.gz", ".tsv", ".tsv.gz", ".vcf.bgz", ".vcf.gz",
    ".bed", ".txt", ".cache", "",
])


@st.composite
def valid_paths(draw: st.DrawFn) -> Path:
    """Generate valid filesystem paths for testing cache_path_for."""
    depth = draw(st.integers(min_value=1, max_value=4))
    segments = [draw(path_segments) for _ in range(depth)]
    filename = draw(path_segments) + draw(file_extensions)
    return Path("/", *segments, filename)


@given(p=valid_paths())
@settings(max_examples=200)
def test_cache_path_appends_suffix(p: Path) -> None:
    """cache_path_for(p) == Path(str(p) + '.vartriage.cache') for any path."""
    result = cache_path_for(p)
    expected = Path(str(p) + ".vartriage.cache")
    assert result == expected


@given(p=valid_paths())
@settings(max_examples=200)
def test_cache_path_same_parent_directory(p: Path) -> None:
    """The cache file always sits in the same parent directory as p."""
    result = cache_path_for(p)
    assert result.parent == p.parent


# ---------------------------------------------------------------------------
# Task 7.2: Property test for cache envelope round-trip (Property 2)
# Validates: Requirements 1.3, 5.3, 7.1
# ---------------------------------------------------------------------------

# Strategy for picklable data: dicts, lists, tuples, strings, floats
picklable_data = st.recursive(
    st.one_of(
        st.none(),
        st.booleans(),
        st.integers(min_value=-1_000_000, max_value=1_000_000),
        st.floats(allow_nan=False, allow_infinity=False),
        st.text(min_size=0, max_size=50),
    ),
    lambda children: st.one_of(
        st.lists(children, max_size=10),
        st.tuples(children, children),
        st.dictionaries(
            st.text(min_size=1, max_size=10),
            children,
            max_size=10,
        ),
    ),
    max_leaves=30,
)

version_strings = st.from_regex(r"[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}", fullmatch=True)
python_version_strings = st.from_regex(r"[0-9]{1,2}\.[0-9]{1,2}", fullmatch=True)


@given(
    data=picklable_data,
    mtime=st.floats(
        min_value=0.0,
        max_value=2_000_000_000.0,
        allow_nan=False,
        allow_infinity=False,
    ),
    vt_version=version_strings,
    py_version=python_version_strings,
)
@settings(max_examples=200)
def test_cache_envelope_round_trip(
    data: object,
    mtime: float,
    vt_version: str,
    py_version: str,
) -> None:
    """Serialize a CacheEnvelope via pickle, deserialize, verify all fields."""
    envelope = CacheEnvelope(
        vartriage_version=vt_version,
        python_version=py_version,
        source_mtime=mtime,
        data=data,
    )

    serialized = pickle.dumps(envelope, protocol=pickle.HIGHEST_PROTOCOL)
    restored = pickle.loads(serialized)  # noqa: S301

    assert isinstance(restored, CacheEnvelope)
    assert restored.vartriage_version == vt_version
    assert restored.python_version == py_version
    assert restored.source_mtime == mtime
    assert restored.data == data


# ---------------------------------------------------------------------------
# Task 7.3: Property test for cache version envelope matching runtime
# (Property 7)
# Validates: Requirements 7.1, 7.2
# ---------------------------------------------------------------------------


@given(
    data=picklable_data,
)
@settings(max_examples=200)
def test_cache_write_stamps_current_versions(data: object) -> None:
    """try_write_cache produces an envelope with the current vartriage_version
    and python_version.
    """
    from vartriage import __version__ as current_vt_version

    expected_py_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    with tempfile.TemporaryDirectory() as tmpdir:
        source = Path(tmpdir) / "source.tsv"
        source.write_text("placeholder", encoding="utf-8")

        try_write_cache(source, data)

        cache_file = cache_path_for(source)
        assert cache_file.exists(), "Cache file should have been written"

        with open(cache_file, "rb") as f:
            envelope: CacheEnvelope = pickle.load(f)  # noqa: S301

        assert envelope.vartriage_version == current_vt_version
        assert envelope.python_version == expected_py_version


# ---------------------------------------------------------------------------
# Task 7.4: Unit tests for cache failure scenarios
# Requirements: 1.4, 2.4, 5.4, 6.4, 7.4
# ---------------------------------------------------------------------------


class TestCacheWriteFailure:
    """Write failure (os.replace raises OSError) logs warning, no exception."""

    def test_write_failure_logs_warning_no_exception(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        source = tmp_path / "data.tsv"
        source.write_text("content", encoding="utf-8")

        with patch("os.replace", side_effect=OSError("disk full")):
            with caplog.at_level(logging.WARNING):
                try_write_cache(source, {"key": "value"})

        assert "Failed to write cache" in caplog.text
        assert not cache_path_for(source).exists()


class TestCorruptedCacheFile:
    """Corrupted cache file triggers warning, deletion, returns None."""

    def test_corrupted_cache_warns_deletes_returns_none(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        source = tmp_path / "data.tsv"
        source.write_text("content", encoding="utf-8")

        cache_file = cache_path_for(source)
        cache_file.write_bytes(b"not valid pickle data at all \x00\x01\x02")

        with caplog.at_level(logging.WARNING):
            result = try_load_cache(source)

        assert result is None
        assert "Failed to deserialize" in caplog.text or "Cannot read" in caplog.text
        assert not cache_file.exists()


class TestVersionMismatchInvalidatesCache:
    """Version mismatch in envelope invalidates cache, returns None."""

    def test_vartriage_version_mismatch_returns_none(
        self, tmp_path: Path
    ) -> None:
        source = tmp_path / "data.tsv"
        source.write_text("content", encoding="utf-8")

        # Write a cache with a fake old version
        envelope = CacheEnvelope(
            vartriage_version="0.0.0-fake",
            python_version=f"{sys.version_info.major}.{sys.version_info.minor}",
            source_mtime=source.stat().st_mtime,
            data={"cached": True},
        )
        cache_file = cache_path_for(source)
        with open(cache_file, "wb") as f:
            pickle.dump(envelope, f, protocol=pickle.HIGHEST_PROTOCOL)

        result = try_load_cache(source)
        assert result is None
        assert not cache_file.exists()

    def test_python_version_mismatch_returns_none(
        self, tmp_path: Path
    ) -> None:
        from vartriage import __version__ as vt_version

        source = tmp_path / "data.tsv"
        source.write_text("content", encoding="utf-8")

        envelope = CacheEnvelope(
            vartriage_version=vt_version,
            python_version="2.7",  # obviously wrong
            source_mtime=source.stat().st_mtime,
            data={"cached": True},
        )
        cache_file = cache_path_for(source)
        with open(cache_file, "wb") as f:
            pickle.dump(envelope, f, protocol=pickle.HIGHEST_PROTOCOL)

        result = try_load_cache(source)
        assert result is None
        assert not cache_file.exists()


class TestPermissionErrorOnRead:
    """Permission error on cache read logs warning, returns None."""

    def test_permission_error_warns_returns_none(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        source = tmp_path / "data.tsv"
        source.write_text("content", encoding="utf-8")

        # Write a valid cache, then make it unreadable
        try_write_cache(source, {"key": "value"})
        cache_file = cache_path_for(source)
        assert cache_file.exists()

        os.chmod(cache_file, 0o000)
        try:
            with caplog.at_level(logging.WARNING):
                result = try_load_cache(source)
            assert result is None
            assert "Cannot read cache" in caplog.text
        finally:
            # File may have been deleted by _delete_cache; restore perms only if still present
            if cache_file.exists():
                os.chmod(cache_file, 0o644)


# ---------------------------------------------------------------------------
# Task 7.5: Unit tests for TabixFrequencyDatabase
# Requirements: 3.2, 3.3, 3.4, 3.5
# ---------------------------------------------------------------------------


class TestTabixMissingIndex:
    """Missing .tbi index raises ReferenceFileError with descriptive message."""

    def test_missing_tbi_raises_reference_file_error(
        self, tmp_path: Path
    ) -> None:
        from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase
        from vartriage.io.exceptions import ReferenceFileError

        vcf_path = tmp_path / "gnomad.vcf.bgz"
        vcf_path.write_bytes(b"fake vcf content")
        # No .tbi file created

        db = TabixFrequencyDatabase()
        with pytest.raises(ReferenceFileError) as exc_info:
            db.load(vcf_path)

        error_msg = str(exc_info.value)
        assert "tabix index file not found" in error_msg
        assert ".tbi" in error_msg


class TestTabixLookupAbsentVariant:
    """Lookup for absent variant returns None."""

    def test_absent_variant_returns_none(self) -> None:
        from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

        db = TabixFrequencyDatabase()
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter([])
        db._tabix = mock_tabix

        results = db.lookup_batch([("chr1", 12345, "A", "G")])
        assert results == [None]


class TestTabixMultiallelicRecord:
    """Multiallelic VCF record returns correct AF for specific alt allele."""

    def test_correct_af_for_second_alt(self) -> None:
        from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

        db = TabixFrequencyDatabase()
        mock_tabix = MagicMock()

        # Multiallelic record: REF=A, ALT=G,T  AF=0.001,0.05
        record = "chr1\t100\t.\tA\tG,T\t.\t.\tAF=0.001,0.05"
        mock_tabix.fetch.return_value = iter([record])
        db._tabix = mock_tabix

        # Query for the second alt allele (T)
        results = db.lookup_batch([("chr1", 100, "A", "T")])
        assert results[0] == pytest.approx(0.05)

    def test_correct_af_for_first_alt(self) -> None:
        from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

        db = TabixFrequencyDatabase()
        mock_tabix = MagicMock()

        record = "chr1\t200\t.\tC\tA,G,T\t.\t.\tAF=0.1,0.02,0.003"
        mock_tabix.fetch.return_value = iter([record])
        db._tabix = mock_tabix

        results = db.lookup_batch([("chr1", 200, "C", "A")])
        assert results[0] == pytest.approx(0.1)


class TestTabixMalformedAFField:
    """Malformed AF field logs warning and returns None."""

    def test_malformed_af_logs_warning_returns_none(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

        db = TabixFrequencyDatabase()
        mock_tabix = MagicMock()

        # AF field has a non-numeric value
        record = "chr1\t100\t.\tA\tG\t.\t.\tAF=NOT_A_NUMBER"
        mock_tabix.fetch.return_value = iter([record])
        db._tabix = mock_tabix

        with caplog.at_level(logging.WARNING, logger="vartriage.annotation.frequency_tabix"):
            results = db.lookup_batch([("chr1", 100, "A", "G")])

        assert results == [None]
        assert "Malformed AF value" in caplog.text


# ---------------------------------------------------------------------------
# Task 7.6: Unit tests for AnnotationEngine backend selection
# Requirements: 4.1, 4.2, 4.3
# ---------------------------------------------------------------------------


class TestAnnotationEngineBackendSelection:
    """AnnotationEngine._build_frequency_db selects correct backend by extension."""

    def _make_engine(self, tmp_path: Path) -> object:
        """Create an AnnotationEngine instance with mocked internals."""
        from vartriage.annotation.engine import AnnotationEngine
        from vartriage.models.config import AnnotationConfig

        with patch.object(AnnotationEngine, "__init__", lambda self, config: None):
            engine = AnnotationEngine.__new__(AnnotationEngine)
            return engine

    def test_vcf_bgz_selects_tabix(self, tmp_path: Path) -> None:
        from vartriage.annotation.engine import AnnotationEngine

        engine = self._make_engine(tmp_path)
        gnomad_path = tmp_path / "gnomad.vcf.bgz"
        gnomad_path.write_bytes(b"fake")

        with patch(
            "vartriage.annotation.frequency_tabix.TabixFrequencyDatabase"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance

            result = engine._build_frequency_db(gnomad_path)

        mock_cls.assert_called_once()
        mock_instance.load.assert_called_once_with(gnomad_path)
        assert result is mock_instance

    def test_vcf_gz_selects_tabix(self, tmp_path: Path) -> None:
        from vartriage.annotation.engine import AnnotationEngine

        engine = self._make_engine(tmp_path)
        gnomad_path = tmp_path / "gnomad.vcf.gz"
        gnomad_path.write_bytes(b"fake")

        with patch(
            "vartriage.annotation.frequency_tabix.TabixFrequencyDatabase"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance

            result = engine._build_frequency_db(gnomad_path)

        mock_cls.assert_called_once()
        assert result is mock_instance

    def test_tsv_does_not_select_tabix(self, tmp_path: Path) -> None:
        from vartriage.annotation.engine import AnnotationEngine

        engine = self._make_engine(tmp_path)
        gnomad_path = tmp_path / "gnomad_frequencies.tsv"
        gnomad_path.write_text("chr1\t100\tA\tG\t0.01\n", encoding="utf-8")

        with patch(
            "vartriage.annotation.frequency_tabix.TabixFrequencyDatabase"
        ) as mock_tabix_cls:
            # Mock the dict/polars backend that will actually be used
            with patch(
                "vartriage.annotation.frequency.DictFrequencyDatabase"
            ) as mock_dict_cls:
                mock_dict_instance = MagicMock()
                mock_dict_cls.return_value = mock_dict_instance

                result = engine._build_frequency_db(gnomad_path)

        mock_tabix_cls.assert_not_called()
        assert result is mock_dict_instance

    def test_backend_name_logged_at_info_level(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        from vartriage.annotation.engine import AnnotationEngine

        engine = self._make_engine(tmp_path)
        gnomad_path = tmp_path / "gnomad.vcf.bgz"
        gnomad_path.write_bytes(b"fake")

        with patch(
            "vartriage.annotation.frequency_tabix.TabixFrequencyDatabase"
        ) as mock_cls:
            mock_instance = MagicMock()
            mock_cls.return_value = mock_instance

            with caplog.at_level(logging.INFO, logger="vartriage.annotation.engine"):
                engine._build_frequency_db(gnomad_path)

        assert "tabix" in caplog.text.lower() or "backend" in caplog.text.lower()
