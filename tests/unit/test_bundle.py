"""Tests for the score bundle downloader module."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from vartriage.bundle._checksums import (ChecksumMismatchError, compute_sha256,
                                         verify_checksum)
from vartriage.bundle._disk import (available_space_bytes, check_disk_space,
                                    format_bytes)
from vartriage.bundle.config import BundleConfig
from vartriage.bundle.downloader import BundleDownloader, DownloadError
from vartriage.bundle.manifest import BundleManifest
from vartriage.bundle.registry import BundleEntry, BundleRegistry
from vartriage.bundle.storage import BundleStorage
from vartriage.bundle.transformer import (TRANSFORM_REGISTRY,
                                          PassthroughTransformer,
                                          get_transformer)


class TestBundleRegistry:
    def test_load_bundled_registry(self) -> None:
        reg = BundleRegistry.load()
        assert reg.version == "1.0.0"
        assert len(reg.bundles) >= 5

    def test_get_existing_bundle(self) -> None:
        reg = BundleRegistry.load()
        entry = reg.get("clinvar", "grch38")
        assert entry is not None
        assert entry.display_name == "ClinVar"
        assert "grch38" in entry.source_urls

    def test_get_nonexistent_bundle(self) -> None:
        reg = BundleRegistry.load()
        assert reg.get("nonexistent", "grch38") is None

    def test_get_wrong_build(self) -> None:
        reg = BundleRegistry.load()
        # gnomad-exomes-chr22 only supports grch38
        assert reg.get("gnomad-exomes-chr22", "grch37") is None

    def test_available_for_build(self) -> None:
        reg = BundleRegistry.load()
        grch38 = reg.available_for_build("grch38")
        assert len(grch38) >= 5

        grch37 = reg.available_for_build("grch37")
        # ClinVar and SpliceAI support grch37
        assert len(grch37) >= 2

    def test_bundle_names(self) -> None:
        reg = BundleRegistry.load()
        names = reg.bundle_names()
        assert "clinvar" in names
        assert names == sorted(names)

    def test_parse_invalid_json(self) -> None:
        with pytest.raises(ValueError, match="Invalid registry JSON"):
            BundleRegistry._parse("not valid json")

    def test_parse_missing_bundles(self) -> None:
        reg = BundleRegistry._parse('{"version": "1.0"}')
        assert len(reg.bundles) == 0


class TestBundleManifest:
    def test_save_and_load(self, tmp_path: Path) -> None:
        manifest = BundleManifest(
            bundle_name="clinvar",
            version="2026-07-01",
            genome_build="grch38",
            download_timestamp="2026-07-13T10:00:00Z",
            source_url="https://example.com/clinvar.vcf.gz",
            raw_checksum="sha256:abc123",
            transformed_filename="clinvar.tsv",
            raw_filename="clinvar.vcf.gz",
        )
        manifest_path = tmp_path / "manifest.json"
        manifest.save(manifest_path)

        loaded = BundleManifest.load(manifest_path)
        assert loaded.bundle_name == "clinvar"
        assert loaded.version == "2026-07-01"
        assert loaded.genome_build == "grch38"
        assert loaded.raw_checksum == "sha256:abc123"

    def test_load_nonexistent_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            BundleManifest.load(tmp_path / "missing.json")

    def test_load_invalid_json(self, tmp_path: Path) -> None:
        bad_file = tmp_path / "bad.json"
        bad_file.write_text("not json")
        with pytest.raises(ValueError, match="Invalid manifest JSON"):
            BundleManifest.load(bad_file)


class TestBundleStorage:
    def test_bundle_dir_path(self, tmp_path: Path) -> None:
        storage = BundleStorage(tmp_path)
        path = storage.bundle_dir("grch38", "clinvar")
        assert path == tmp_path / "grch38" / "clinvar"

    def test_not_installed_by_default(self, tmp_path: Path) -> None:
        storage = BundleStorage(tmp_path)
        assert storage.is_installed("grch38", "clinvar") is False
        assert storage.installed_version("grch38", "clinvar") is None

    def test_resolve_path_returns_none_when_not_installed(self, tmp_path: Path) -> None:
        storage = BundleStorage(tmp_path)
        assert storage.resolve_path("grch38", "clinvar") is None

    def test_resolve_path_after_install(self, tmp_path: Path) -> None:
        storage = BundleStorage(tmp_path)
        storage.ensure_dirs("grch38", "clinvar")

        # Create a fake transformed file
        bundle_dir = storage.bundle_dir("grch38", "clinvar")
        (bundle_dir / "clinvar.tsv").write_text("chrom\tpos\n")

        # Write manifest
        manifest = BundleManifest(
            bundle_name="clinvar",
            version="2026-07-01",
            genome_build="grch38",
            transformed_filename="clinvar.tsv",
        )
        manifest.save(storage.manifest_path("grch38", "clinvar"))

        resolved = storage.resolve_path("grch38", "clinvar")
        assert resolved is not None
        assert resolved.name == "clinvar.tsv"

    def test_list_installed_empty(self, tmp_path: Path) -> None:
        storage = BundleStorage(tmp_path)
        assert storage.list_installed("grch38") == []

    def test_disk_usage(self, tmp_path: Path) -> None:
        storage = BundleStorage(tmp_path)
        storage.ensure_dirs("grch38", "clinvar")
        (storage.bundle_dir("grch38", "clinvar") / "test.tsv").write_text("data")

        usage = storage.disk_usage("grch38")
        assert "clinvar" in usage
        assert usage["clinvar"] > 0


class TestChecksums:
    def test_compute_sha256(self, tmp_path: Path) -> None:
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello world")
        checksum = compute_sha256(test_file)
        assert checksum.startswith("sha256:")
        assert len(checksum) == len("sha256:") + 64

    def test_verify_checksum_matches(self, tmp_path: Path) -> None:
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello world")
        checksum = compute_sha256(test_file)
        assert verify_checksum(test_file, checksum) is True

    def test_verify_checksum_mismatch(self, tmp_path: Path) -> None:
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello world")
        assert verify_checksum(test_file, "sha256:wrong") is False

    def test_verify_empty_expected_skips(self, tmp_path: Path) -> None:
        test_file = tmp_path / "test.txt"
        test_file.write_text("anything")
        assert verify_checksum(test_file, "") is True

    def test_checksum_mismatch_error(self) -> None:
        err = ChecksumMismatchError(
            Path("/tmp/test"), "sha256:expected", "sha256:actual"
        )
        assert "expected" in str(err)
        assert "actual" in str(err)


class TestDiskUtilities:
    def test_available_space_positive(self) -> None:
        space = available_space_bytes(Path.home())
        assert space > 0

    def test_format_bytes_values(self) -> None:
        assert format_bytes(0) == "0 B"
        assert format_bytes(512) == "512 B"
        assert "KB" in format_bytes(2048)
        assert "MB" in format_bytes(5 * 1024 * 1024)
        assert "GB" in format_bytes(3 * 1024 * 1024 * 1024)

    def test_format_bytes_negative(self) -> None:
        assert format_bytes(-1) == "0 B"


class TestBundleConfig:
    def test_default_config(self) -> None:
        config = BundleConfig.load()
        assert config.default_build == "grch38"
        assert config.download_concurrency == 2
        assert config.auto_verify is True

    def test_env_override(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("VARTRIAGE_DEFAULT_BUILD", "grch37")
        config = BundleConfig.load()
        assert config.default_build == "grch37"

    def test_storage_path_env(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        monkeypatch.setenv("VARTRIAGE_BUNDLE_STORAGE", str(tmp_path))
        config = BundleConfig.load()
        assert config.storage_path == tmp_path


class TestTransformerRegistry:
    def test_registry_has_all_types(self) -> None:
        assert "vcf_to_tsv" in TRANSFORM_REGISTRY
        assert "csv_to_tsv" in TRANSFORM_REGISTRY
        assert "spliceai_vcf" in TRANSFORM_REGISTRY
        assert "none" in TRANSFORM_REGISTRY

    def test_get_clinvar_transformer(self) -> None:
        t = get_transformer("vcf_to_tsv", "clinvar")
        assert type(t).__name__ == "ClinvarVcfTransformer"

    def test_get_unknown_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown transform type"):
            get_transformer("nonexistent")

    def test_passthrough_decompresses_gz(self, tmp_path: Path) -> None:
        import gzip

        source = tmp_path / "test.gtf.gz"
        content = b"chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tgene_id\n" * 10
        with gzip.open(source, "wb") as f:
            f.write(content)

        dest = tmp_path / "test.gtf"
        transformer = PassthroughTransformer()
        result = transformer.transform(source, dest, "grch38")

        assert dest.exists()
        assert result.rows_written >= 9
        assert dest.read_bytes() == content


class TestBundleDownloader:
    """Tests for BundleDownloader using monkeypatched urlopen."""

    def test_successful_download(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Mock urlopen returning data, verify file written and .partial removed."""
        import io
        from unittest.mock import MagicMock

        fake_response = MagicMock()
        fake_response.status = 200
        fake_response.read = MagicMock(side_effect=[b"hello world data", b""])

        def mock_urlopen(request: object, timeout: object = None) -> MagicMock:
            return fake_response

        monkeypatch.setattr("vartriage.bundle.downloader.urlopen", mock_urlopen)

        dest = tmp_path / "test_file.txt"
        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=0, show_progress=False
        )
        result = downloader.download(url="http://example.com/test.txt", dest=dest)

        assert dest.exists()
        assert dest.read_bytes() == b"hello world data"
        assert not Path(str(dest) + ".partial").exists()
        assert result.bytes_downloaded == 16
        assert result.checksum_verified is True

    def test_retry_on_transient_error(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Mock urlopen raising HTTPError(500) then succeeding."""
        from unittest.mock import MagicMock
        from urllib.error import HTTPError

        call_count = 0

        def mock_urlopen(request: object, timeout: object = None) -> MagicMock:
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise HTTPError(
                    "http://example.com/f.txt",
                    500,
                    "Server Error",
                    {},
                    None,  # type: ignore[arg-type]
                )
            resp = MagicMock()
            resp.status = 200
            resp.read = MagicMock(side_effect=[b"data", b""])
            return resp

        monkeypatch.setattr("vartriage.bundle.downloader.urlopen", mock_urlopen)
        # Patch time.sleep to avoid actual waiting
        monkeypatch.setattr("vartriage.bundle.downloader.time.sleep", lambda _: None)

        dest = tmp_path / "retry_file.txt"
        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=2, show_progress=False
        )
        result = downloader.download(url="http://example.com/f.txt", dest=dest)

        assert dest.exists()
        assert result.bytes_downloaded == 4
        assert call_count == 2

    def test_checksum_mismatch_raises(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Provide wrong expected checksum, verify DownloadError raised."""
        from unittest.mock import MagicMock

        fake_response = MagicMock()
        fake_response.status = 200
        fake_response.read = MagicMock(side_effect=[b"some data", b""])

        def mock_urlopen(request: object, timeout: object = None) -> MagicMock:
            return fake_response

        monkeypatch.setattr("vartriage.bundle.downloader.urlopen", mock_urlopen)

        dest = tmp_path / "checksum_file.txt"
        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=0, show_progress=False
        )

        with pytest.raises(DownloadError, match="Checksum mismatch"):
            downloader.download(
                url="http://example.com/f.txt",
                dest=dest,
                expected_checksum="sha256:0000000000000000000000000000000000000000000000000000000000000000",
            )


class TestDiskSpaceInsufficient:
    """Test check_disk_space raises on insufficient space."""

    def test_check_disk_space_insufficient_raises(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Monkeypatch available_space_bytes to return less than required."""
        monkeypatch.setattr(
            "vartriage.bundle._disk.available_space_bytes",
            lambda path: 100,
        )
        with pytest.raises(OSError, match="Insufficient disk space"):
            check_disk_space(tmp_path, 1_000_000)


class TestParallelDownload:
    """Tests for download_many parallel download function."""

    def test_download_many_empty_requests(self) -> None:
        """Empty request list returns empty result immediately."""
        from vartriage.bundle.downloader import download_many

        result = download_many(requests=[], concurrency=2)
        assert result.all_succeeded is True
        assert result.success_count == 0
        assert result.error_count == 0
        assert result.total_bytes == 0

    def test_download_many_two_files(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Download two files in parallel with mocked urlopen."""
        from unittest.mock import MagicMock

        from vartriage.bundle.downloader import DownloadRequest, download_many

        call_count = 0

        def mock_urlopen(request: object, timeout: object = None) -> MagicMock:
            nonlocal call_count
            call_count += 1
            resp = MagicMock()
            resp.status = 200
            resp.read = MagicMock(side_effect=[b"file content here", b""])
            return resp

        monkeypatch.setattr("vartriage.bundle.downloader.urlopen", mock_urlopen)

        requests = [
            DownloadRequest(
                url="http://example.com/a.txt",
                dest=tmp_path / "a.txt",
                label="file-a",
            ),
            DownloadRequest(
                url="http://example.com/b.txt",
                dest=tmp_path / "b.txt",
                label="file-b",
            ),
        ]

        result = download_many(requests, concurrency=2, show_progress=False)

        assert result.all_succeeded is True
        assert result.success_count == 2
        assert result.error_count == 0
        assert "file-a" in result.results
        assert "file-b" in result.results
        assert (tmp_path / "a.txt").exists()
        assert (tmp_path / "b.txt").exists()
        assert result.total_bytes == 34  # 17 bytes * 2 files

    def test_download_many_partial_failure(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """One download fails, the other succeeds — partial results returned."""
        from unittest.mock import MagicMock
        from urllib.error import HTTPError

        from vartriage.bundle.downloader import DownloadRequest, download_many

        call_urls: list[str] = []

        def mock_urlopen(request: object, timeout: object = None) -> MagicMock:
            url = request.full_url if hasattr(request, "full_url") else str(request)
            call_urls.append(url)
            if "fail" in url:
                raise HTTPError(url, 403, "Forbidden", {}, None)  # type: ignore[arg-type]
            resp = MagicMock()
            resp.status = 200
            resp.read = MagicMock(side_effect=[b"ok", b""])
            return resp

        monkeypatch.setattr("vartriage.bundle.downloader.urlopen", mock_urlopen)
        # Speed up retries
        monkeypatch.setattr("vartriage.bundle.downloader.time.sleep", lambda _: None)

        requests = [
            DownloadRequest(
                url="http://example.com/good.txt",
                dest=tmp_path / "good.txt",
                label="good",
            ),
            DownloadRequest(
                url="http://example.com/fail.txt",
                dest=tmp_path / "fail.txt",
                label="bad",
            ),
        ]

        result = download_many(
            requests, concurrency=2, max_retries=0, show_progress=False
        )

        assert result.all_succeeded is False
        assert result.success_count == 1
        assert result.error_count == 1
        assert "good" in result.results
        assert "bad" in result.errors
        assert (tmp_path / "good.txt").exists()
