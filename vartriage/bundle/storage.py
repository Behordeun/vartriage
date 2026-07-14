"""Bundle storage: on-disk layout and path resolution.

Manages the directory structure where downloaded and transformed
reference files are stored. Default location: ~/.vartriage/bundles/
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from vartriage.bundle.manifest import BundleManifest


class BundleStorage:
    """Manages on-disk storage layout for score bundles.

    Directory structure::

        {base_path}/
            {build}/
                {bundle_name}/
                    raw/           # Original downloaded files
                    manifest.json  # Version, checksums, timestamps
                    {output_file}  # Transformed TSV ready for pipeline
    """

    def __init__(self, base_path: Optional[Path] = None) -> None:
        """Initialize storage manager.

        Parameters
        ----------
        base_path : Path, optional
            Root directory for bundle storage. Defaults to
            ~/.vartriage/bundles or VARTRIAGE_BUNDLE_STORAGE env var.
        """
        if base_path is not None:
            self._base = base_path
        else:
            env_path = os.environ.get("VARTRIAGE_BUNDLE_STORAGE")
            if env_path:
                self._base = Path(env_path)
            else:
                self._base = Path.home() / ".vartriage" / "bundles"

    @property
    def base_path(self) -> Path:
        """Root storage directory."""
        return self._base

    def bundle_dir(self, build: str, bundle_name: str) -> Path:
        """Path to a specific bundle's directory."""
        return self._base / build.lower() / bundle_name

    def raw_dir(self, build: str, bundle_name: str) -> Path:
        """Path to the raw/ subdirectory for original downloads."""
        return self.bundle_dir(build, bundle_name) / "raw"

    def manifest_path(self, build: str, bundle_name: str) -> Path:
        """Path to the bundle's manifest.json."""
        return self.bundle_dir(build, bundle_name) / "manifest.json"

    def is_installed(self, build: str, bundle_name: str) -> bool:
        """Check if a bundle is installed (has a valid manifest)."""
        manifest = self.manifest_path(build, bundle_name)
        return manifest.exists()

    def installed_version(self, build: str, bundle_name: str) -> Optional[str]:
        """Return the installed version, or None if not installed."""
        try:
            manifest = BundleManifest.load(self.manifest_path(build, bundle_name))
            return manifest.version
        except (FileNotFoundError, ValueError):
            return None

    def resolve_path(self, build: str, bundle_name: str) -> Optional[Path]:
        """Resolve the path to the transformed file if installed.

        Returns None if the bundle is not installed or the
        transformed file doesn't exist on disk.
        """
        try:
            manifest = BundleManifest.load(self.manifest_path(build, bundle_name))
        except (FileNotFoundError, ValueError):
            return None

        if not manifest.transformed_filename:
            return None

        transformed = (
            self.bundle_dir(build, bundle_name) / manifest.transformed_filename
        )
        if transformed.exists():
            return transformed
        return None

    def ensure_dirs(self, build: str, bundle_name: str) -> None:
        """Create the bundle directory structure if it doesn't exist."""
        self.raw_dir(build, bundle_name).mkdir(parents=True, exist_ok=True)

    def disk_usage(self, build: str = "grch38") -> dict[str, int]:
        """Calculate disk usage per bundle for a given build.

        Returns a dict mapping bundle name to total bytes on disk.
        """
        build_dir = self._base / build.lower()
        usage: dict[str, int] = {}

        if not build_dir.exists():
            return usage

        for bundle_dir in build_dir.iterdir():
            if bundle_dir.is_dir():
                total = sum(
                    f.stat().st_size for f in bundle_dir.rglob("*") if f.is_file()
                )
                usage[bundle_dir.name] = total

        return usage

    def list_installed(self, build: str = "grch38") -> list[BundleManifest]:
        """List all installed bundles for a given build."""
        build_dir = self._base / build.lower()
        manifests: list[BundleManifest] = []

        if not build_dir.exists():
            return manifests

        for bundle_dir in sorted(build_dir.iterdir()):
            manifest_file = bundle_dir / "manifest.json"
            if manifest_file.exists():
                try:
                    manifests.append(BundleManifest.load(manifest_file))
                except (ValueError, FileNotFoundError):
                    continue

        return manifests
