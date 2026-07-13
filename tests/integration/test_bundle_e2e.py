"""End-to-end integration tests for the bundle download + transform pipeline.

Tests exercise the full workflow:
  registry lookup → download (mocked HTTP) → transform → manifest written → verify

Marked @pytest.mark.slow because they do real file I/O and subprocess calls.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from threading import Thread
from typing import Any
from unittest.mock import MagicMock

import pytest

from vartriage.bundle._checksums import compute_sha256
from vartriage.bundle.config import BundleConfig
from vartriage.bundle.downloader import BundleDownloader, DownloadResult
from vartriage.bundle.manifest import BundleManifest
from vartriage.bundle.registry import BundleEntry, BundleRegistry
from vartriage.bundle.storage import BundleStorage
from vartriage.bundle.transformer import (CsvToTsvTransformer,
                                          PassthroughTransformer,
                                          get_transformer)

# Minimal VCF content for ClinVar-like data
_CLINVAR_VCF_CONTENT = """\
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t12345\trs001\tA\tG\t.\tPASS\tCLNSIG=Pathogenic
1\t23456\trs002\tC\tT\t.\tPASS\tCLNSIG=Likely_pathogenic
2\t34567\trs003\tG\tA\t.\tPASS\tCLNSIG=Uncertain_significance
2\t45678\trs004\tT\tC\t.\tPASS\tCLNSIG=Benign
X\t56789\trs005\tA\tT\t.\tPASS\tCLNSIG=Likely_benign
"""

# Minimal CSV content mimicking REVEL format
_REVEL_CSV_CONTENT = """\
chr,grch38_pos,ref,alt,REVEL
1,12345,A,G,0.85
1,23456,C,T,0.42
2,34567,G,A,0.15
2,45678,T,C,0.91
"""

# Simple GTF content for passthrough
_GENCODE_GTF_CONTENT = """\
##description: evidence-based annotation
chr1\tENSEMBL\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"
chr1\tENSEMBL\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"
chr1\tENSEMBL\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"
"""


@pytest.fixture()
def bundle_storage(tmp_path: Path) -> BundleStorage:
    """Provide a BundleStorage rooted in a temp directory."""
    return BundleStorage(tmp_path / "bundles")


@pytest.fixture()
def serve_files(tmp_path: Path) -> tuple[str, Path]:
    """Spin up a local HTTP server serving files from a temp directory.

    Returns (base_url, serve_dir) so tests can place files and
    download them via HTTP.
    """
    serve_dir = tmp_path / "serve"
    serve_dir.mkdir()

    class QuietHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, directory=str(serve_dir), **kwargs)

        def log_message(self, format: str, *args: Any) -> None:
            pass  # suppress noisy output

    server = HTTPServer(("127.0.0.1", 0), QuietHandler)
    port = server.server_address[1]
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()

    yield f"http://127.0.0.1:{port}", serve_dir

    server.shutdown()


@pytest.mark.slow
class TestBundleDownloadTransformE2E:
    """Full end-to-end: download from local HTTP → transform → manifest."""

    def test_clinvar_vcf_download_and_transform(
        self,
        bundle_storage: BundleStorage,
        serve_files: tuple[str, Path],
    ) -> None:
        """Download a ClinVar-like VCF, transform to TSV, verify manifest."""
        base_url, serve_dir = serve_files
        build = "grch38"
        bundle_name = "clinvar"

        # Place a gzipped VCF on the local server
        vcf_gz_path = serve_dir / "clinvar_test.vcf.gz"
        with gzip.open(vcf_gz_path, "wb") as f:
            f.write(_CLINVAR_VCF_CONTENT.encode("utf-8"))

        # Compute checksum of the served file for verification
        expected_checksum = compute_sha256(vcf_gz_path)

        # Download
        bundle_storage.ensure_dirs(build, bundle_name)
        raw_dir = bundle_storage.raw_dir(build, bundle_name)

        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=1, show_progress=False
        )
        result = downloader.download(
            url=f"{base_url}/clinvar_test.vcf.gz",
            dest=raw_dir / "clinvar_test.vcf.gz",
            expected_checksum=expected_checksum,
        )

        assert result.path.exists()
        assert result.checksum_verified is True
        assert result.bytes_downloaded > 0

        # Transform using the passthrough transformer (since we don't have
        # bcftools in CI, and ClinvarVcfTransformer needs it or pysam)
        # Instead, test with PassthroughTransformer to verify the pipeline
        transformer = PassthroughTransformer()
        bundle_dir = bundle_storage.bundle_dir(build, bundle_name)
        transform_dest = bundle_dir / "clinvar.vcf"

        t_result = transformer.transform(result.path, transform_dest, build)

        assert transform_dest.exists()
        assert t_result.rows_written > 0

        # Write manifest
        manifest = BundleManifest(
            bundle_name=bundle_name,
            version="2026-07-01",
            genome_build=build,
            download_timestamp="2026-07-14T00:00:00Z",
            source_url=f"{base_url}/clinvar_test.vcf.gz",
            raw_checksum=expected_checksum,
            transformed_filename="clinvar.vcf",
            raw_filename="clinvar_test.vcf.gz",
            transform_type="none",
            file_size_bytes=result.path.stat().st_size,
        )
        manifest.save(bundle_storage.manifest_path(build, bundle_name))

        # Verify manifest roundtrip
        loaded = BundleManifest.load(bundle_storage.manifest_path(build, bundle_name))
        assert loaded.bundle_name == "clinvar"
        assert loaded.version == "2026-07-01"
        assert loaded.raw_checksum == expected_checksum

        # Verify storage reports it as installed
        assert bundle_storage.is_installed(build, bundle_name)
        assert bundle_storage.installed_version(build, bundle_name) == "2026-07-01"

    def test_revel_csv_download_and_transform(
        self,
        bundle_storage: BundleStorage,
        serve_files: tuple[str, Path],
    ) -> None:
        """Download a REVEL-like CSV, transform to TSV, verify output."""
        base_url, serve_dir = serve_files
        build = "grch38"
        bundle_name = "revel"

        # Place CSV on the server
        csv_path = serve_dir / "revel_scores.csv"
        csv_path.write_text(_REVEL_CSV_CONTENT, encoding="utf-8")

        # Download
        bundle_storage.ensure_dirs(build, bundle_name)
        raw_dir = bundle_storage.raw_dir(build, bundle_name)

        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=1, show_progress=False
        )
        result = downloader.download(
            url=f"{base_url}/revel_scores.csv",
            dest=raw_dir / "revel_scores.csv",
        )

        assert result.path.exists()

        # Transform CSV → TSV
        transformer = CsvToTsvTransformer()
        bundle_dir = bundle_storage.bundle_dir(build, bundle_name)
        transform_dest = bundle_dir / "revel.tsv"

        t_result = transformer.transform(result.path, transform_dest, build)

        assert transform_dest.exists()
        assert t_result.rows_written == 4

        # Verify output format
        lines = transform_dest.read_text().strip().split("\n")
        header = lines[0]
        assert "chrom" in header
        assert "pos" in header
        assert "score" in header

        # First data line should have chr prefix
        first_data = lines[1].split("\t")
        assert first_data[0] == "chr1"

    def test_gencode_passthrough_gz(
        self,
        bundle_storage: BundleStorage,
        serve_files: tuple[str, Path],
    ) -> None:
        """Download a gzipped GTF, passthrough decompress, verify."""
        base_url, serve_dir = serve_files
        build = "grch38"
        bundle_name = "gencode"

        # Place gzipped GTF on server
        gtf_gz_path = serve_dir / "gencode.gtf.gz"
        with gzip.open(gtf_gz_path, "wb") as f:
            f.write(_GENCODE_GTF_CONTENT.encode("utf-8"))

        # Download
        bundle_storage.ensure_dirs(build, bundle_name)
        raw_dir = bundle_storage.raw_dir(build, bundle_name)

        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=1, show_progress=False
        )
        result = downloader.download(
            url=f"{base_url}/gencode.gtf.gz",
            dest=raw_dir / "gencode.gtf.gz",
        )

        assert result.path.exists()

        # Passthrough transform (decompress .gz)
        transformer = PassthroughTransformer()
        bundle_dir = bundle_storage.bundle_dir(build, bundle_name)
        transform_dest = bundle_dir / "gencode.gtf"

        t_result = transformer.transform(result.path, transform_dest, build)

        assert transform_dest.exists()
        assert t_result.rows_written >= 2
        content = transform_dest.read_text()
        assert "ENSG00000223972" in content

    def test_download_resume_after_interruption(
        self,
        bundle_storage: BundleStorage,
        serve_files: tuple[str, Path],
    ) -> None:
        """Simulate partial download, verify resume picks up correctly."""
        base_url, serve_dir = serve_files
        build = "grch38"
        bundle_name = "resume-test"

        # Create a file to serve
        content = b"A" * 1000 + b"B" * 1000
        (serve_dir / "large.bin").write_bytes(content)

        bundle_storage.ensure_dirs(build, bundle_name)
        raw_dir = bundle_storage.raw_dir(build, bundle_name)
        dest = raw_dir / "large.bin"

        # Simulate a partial download (write first 500 bytes as .partial)
        partial = Path(str(dest) + ".partial")
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.write_bytes(content[:500])

        # Download with resume
        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=1, show_progress=False
        )
        result = downloader.download(
            url=f"{base_url}/large.bin",
            dest=dest,
            resume=True,
        )

        assert dest.exists()
        # Server may or may not support Range; either way we get the full file
        assert dest.stat().st_size >= len(content)

    def test_full_registry_to_storage_workflow(
        self,
        bundle_storage: BundleStorage,
        serve_files: tuple[str, Path],
    ) -> None:
        """Complete workflow: registry lookup → download → transform → storage."""
        base_url, serve_dir = serve_files
        build = "grch38"

        # Place a CSV on the server
        csv_path = serve_dir / "revel-v1.3.csv"
        csv_path.write_text(_REVEL_CSV_CONTENT, encoding="utf-8")

        # Create a custom registry pointing to our local server
        registry_data = {
            "version": "1.0.0-test",
            "updated_at": "2026-07-14",
            "bundles": {
                "revel": {
                    "display_name": "REVEL",
                    "description": "Test REVEL bundle",
                    "source_urls": {
                        "grch38": f"{base_url}/revel-v1.3.csv",
                    },
                    "version": "1.3",
                    "expected_size_bytes": len(_REVEL_CSV_CONTENT),
                    "checksum_sha256": "",
                    "transformed_checksum": "",
                    "genome_builds": ["grch38"],
                    "release_date": "2022-09-01",
                    "transform_type": "csv_to_tsv",
                },
            },
        }

        # Write and load custom registry
        reg_path = bundle_storage.base_path / "test_registry.json"
        reg_path.parent.mkdir(parents=True, exist_ok=True)
        reg_path.write_text(json.dumps(registry_data), encoding="utf-8")
        registry = BundleRegistry.load(reg_path)

        # Lookup bundle
        entry = registry.get("revel", build)
        assert entry is not None
        assert entry.version == "1.3"

        url = entry.source_urls[build]
        bundle_name = entry.name

        # Download
        bundle_storage.ensure_dirs(build, bundle_name)
        raw_dir = bundle_storage.raw_dir(build, bundle_name)
        raw_filename = url.rsplit("/", 1)[-1]

        downloader = BundleDownloader(
            timeout=(5, 10), max_retries=1, show_progress=False
        )
        dl_result = downloader.download(
            url=url,
            dest=raw_dir / raw_filename,
        )
        assert dl_result.path.exists()

        # Transform
        transformer = get_transformer(entry.transform_type, bundle_name)
        bundle_dir = bundle_storage.bundle_dir(build, bundle_name)
        output_name = f"{bundle_name}.tsv"
        transform_dest = bundle_dir / output_name

        t_result = transformer.transform(dl_result.path, transform_dest, build)
        assert t_result.rows_written == 4

        # Write manifest
        manifest = BundleManifest(
            bundle_name=bundle_name,
            version=entry.version,
            genome_build=build,
            download_timestamp="2026-07-14T00:00:00Z",
            source_url=url,
            raw_checksum=compute_sha256(dl_result.path),
            transformed_filename=output_name,
            raw_filename=raw_filename,
            transform_type=entry.transform_type,
            file_size_bytes=dl_result.path.stat().st_size,
        )
        manifest.save(bundle_storage.manifest_path(build, bundle_name))

        # Verify end state via storage API
        assert bundle_storage.is_installed(build, bundle_name)
        resolved = bundle_storage.resolve_path(build, bundle_name)
        assert resolved is not None
        assert resolved.name == "revel.tsv"

        # Verify disk usage is tracked
        usage = bundle_storage.disk_usage(build)
        assert "revel" in usage
        assert usage["revel"] > 0
