"""Per-bundle manifest: tracks installed version, checksums, and metadata.

Each installed bundle has a manifest.json in its directory recording
when it was downloaded, what version, and the integrity checksums.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class BundleManifest:
    """On-disk manifest for an installed bundle.

    Parameters
    ----------
    bundle_name : str
        Bundle identifier matching the registry name.
    version : str
        Installed version string.
    genome_build : str
        Genome build this bundle was downloaded for.
    download_timestamp : str
        ISO 8601 timestamp of when the download completed.
    source_url : str
        URL the raw file was downloaded from.
    raw_checksum : str
        SHA-256 hex digest of the raw downloaded file.
    transformed_checksum : str
        SHA-256 hex digest of the transformed output file.
    transformed_filename : str
        Filename of the transformed output (e.g., "clinvar.tsv").
    raw_filename : str
        Filename of the raw download in the raw/ subdirectory.
    transform_type : str
        Transformation applied ("vcf_to_tsv", "csv_to_tsv", "none", etc.).
    file_size_bytes : int
        Size of the transformed file in bytes.
    """

    bundle_name: str
    version: str
    genome_build: str
    download_timestamp: str = ""
    source_url: str = ""
    raw_checksum: str = ""
    transformed_checksum: str = ""
    transformed_filename: str = ""
    raw_filename: str = ""
    transform_type: str = "none"
    file_size_bytes: int = 0

    def save(self, path: Path) -> None:
        """Write manifest to a JSON file.

        Parameters
        ----------
        path : Path
            File path to write the manifest JSON.
        """
        data = {
            "bundle_name": self.bundle_name,
            "version": self.version,
            "genome_build": self.genome_build,
            "download_timestamp": self.download_timestamp,
            "source_url": self.source_url,
            "raw_checksum": self.raw_checksum,
            "transformed_checksum": self.transformed_checksum,
            "transformed_filename": self.transformed_filename,
            "raw_filename": self.raw_filename,
            "transform_type": self.transform_type,
            "file_size_bytes": self.file_size_bytes,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "BundleManifest":
        """Load manifest from a JSON file.

        Parameters
        ----------
        path : Path
            Path to the manifest JSON file.

        Returns
        -------
        BundleManifest
            Parsed manifest instance.

        Raises
        ------
        FileNotFoundError
            If the manifest file doesn't exist.
        ValueError
            If the JSON is invalid.
        """
        if not path.exists():
            raise FileNotFoundError(f"Manifest not found: {path}")

        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid manifest JSON at {path}: {exc}") from exc

        return cls(
            bundle_name=data.get("bundle_name", ""),
            version=data.get("version", ""),
            genome_build=data.get("genome_build", ""),
            download_timestamp=data.get("download_timestamp", ""),
            source_url=data.get("source_url", ""),
            raw_checksum=data.get("raw_checksum", ""),
            transformed_checksum=data.get("transformed_checksum", ""),
            transformed_filename=data.get("transformed_filename", ""),
            raw_filename=data.get("raw_filename", ""),
            transform_type=data.get("transform_type", "none"),
            file_size_bytes=data.get("file_size_bytes", 0),
        )
