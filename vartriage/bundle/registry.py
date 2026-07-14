"""Bundle registry: metadata for available score bundles.

The registry is a JSON manifest shipped with the package and
updatable from GitHub releases. It defines each bundle's source
URL, expected checksums, genome build support, and transform type.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class BundleEntry:
    """Metadata for a single downloadable score bundle.

    Parameters
    ----------
    name : str
        Machine identifier (e.g., "clinvar", "gnomad-exomes").
    display_name : str
        Human-readable name.
    description : str
        Brief description of what this bundle provides.
    source_urls : dict[str, str]
        Mapping of genome build to download URL.
    version : str
        Bundle version string (e.g., "2025-06-01", "v4.1.1").
    expected_size_bytes : int
        Expected download size in bytes (for disk space checks).
    checksum_sha256 : str
        SHA-256 hex digest of the raw downloaded file.
    transformed_checksum : str
        SHA-256 hex digest of the transformed TSV output.
    genome_builds : list[str]
        Supported genome builds (e.g., ["grch37", "grch38"]).
    release_date : str
        ISO 8601 date of the source data release.
    transform_type : str
        Transformation to apply after download. One of:
        "vcf_to_tsv", "csv_to_tsv", "spliceai_vcf", "none".
    """

    name: str
    display_name: str
    description: str
    source_urls: dict[str, str]
    version: str
    expected_size_bytes: int
    checksum_sha256: str
    transformed_checksum: str
    genome_builds: list[str]
    release_date: str
    transform_type: str


@dataclass
class BundleRegistry:
    """Registry of all available score bundles.

    Loads from a JSON file shipped with the package or fetched
    from a remote URL for updates.
    """

    version: str
    updated_at: str
    bundles: dict[str, BundleEntry]

    @classmethod
    def load(cls, path: Optional[Path] = None) -> "BundleRegistry":
        """Load the registry from a JSON file.

        Parameters
        ----------
        path : Path, optional
            Path to registry JSON file. If None, loads the
            bundled registry shipped with the package.

        Returns
        -------
        BundleRegistry
            Parsed registry instance.

        Raises
        ------
        FileNotFoundError
            If the registry file does not exist.
        ValueError
            If the JSON structure is invalid.
        """
        if path is None:
            path = cls._bundled_registry_path()

        try:
            raw = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise FileNotFoundError(
                f"Registry file not found at {path}: {exc}"
            ) from exc

        return cls._parse(raw)

    @classmethod
    def _parse(cls, raw_json: str) -> "BundleRegistry":
        """Parse raw JSON into a BundleRegistry."""
        try:
            data = json.loads(raw_json)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid registry JSON: {exc}") from exc

        if not isinstance(data, dict):
            raise ValueError("Registry must be a JSON object")

        version = data.get("version", "unknown")
        updated_at = data.get("updated_at", "")
        bundles_raw = data.get("bundles", {})

        if not isinstance(bundles_raw, dict):
            raise ValueError("'bundles' must be a JSON object")

        bundles: dict[str, BundleEntry] = {}
        for name, entry_data in bundles_raw.items():
            try:
                bundles[name] = BundleEntry(
                    name=name,
                    display_name=entry_data.get("display_name", name),
                    description=entry_data.get("description", ""),
                    source_urls=entry_data.get("source_urls", {}),
                    version=entry_data.get("version", "unknown"),
                    expected_size_bytes=entry_data.get("expected_size_bytes", 0),
                    checksum_sha256=entry_data.get("checksum_sha256", ""),
                    transformed_checksum=entry_data.get("transformed_checksum", ""),
                    genome_builds=entry_data.get("genome_builds", []),
                    release_date=entry_data.get("release_date", ""),
                    transform_type=entry_data.get("transform_type", "none"),
                )
            except (TypeError, KeyError) as exc:
                raise ValueError(f"Invalid entry for bundle '{name}': {exc}") from exc

        return cls(version=version, updated_at=updated_at, bundles=bundles)

    def get(self, name: str, build: str = "grch38") -> Optional[BundleEntry]:
        """Look up a bundle by name, filtered by genome build.

        Returns None if the bundle doesn't exist or doesn't
        support the requested build.
        """
        entry = self.bundles.get(name)
        if entry is None:
            return None
        if build.lower() not in [b.lower() for b in entry.genome_builds]:
            return None
        return entry

    def available_for_build(self, build: str = "grch38") -> list[BundleEntry]:
        """Return all bundles that support a given genome build."""
        build_lower = build.lower()
        return [
            entry
            for entry in self.bundles.values()
            if build_lower in [b.lower() for b in entry.genome_builds]
        ]

    def bundle_names(self) -> list[str]:
        """Return sorted list of all bundle names."""
        return sorted(self.bundles.keys())

    @staticmethod
    def _bundled_registry_path() -> Path:
        """Resolve the path to the bundled registry.json."""
        pkg_dir = Path(__file__).parent
        return pkg_dir / "registry.json"
