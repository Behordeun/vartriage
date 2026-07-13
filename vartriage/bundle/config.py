"""Bundle configuration: TOML parsing with env var overrides.

Reads user preferences from ~/.vartriage/config.toml and
overrides from environment variables.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Python 3.11+ has tomllib in stdlib; for 3.10 we use tomli
import sys

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        tomllib = None  # type: ignore[assignment,unused-ignore]


@dataclass
class BundleConfig:
    """User configuration for the bundle subsystem.

    Parameters
    ----------
    default_build : str
        Default genome build for downloads (default: "grch38").
    download_concurrency : int
        Number of parallel downloads (default: 2).
    storage_path : Path
        Root storage directory for bundles.
    auto_verify : bool
        Verify checksums automatically after download (default: True).
    http_proxy : str
        HTTP proxy URL, empty string if none.
    https_proxy : str
        HTTPS proxy URL, empty string if none.
    """

    default_build: str = "grch38"
    download_concurrency: int = 2
    storage_path: Path = field(
        default_factory=lambda: Path.home() / ".vartriage" / "bundles"
    )
    auto_verify: bool = True
    http_proxy: str = ""
    https_proxy: str = ""

    @classmethod
    def load(cls, path: Optional[Path] = None) -> "BundleConfig":
        """Load config from TOML file with env var overrides.

        Parameters
        ----------
        path : Path, optional
            Path to config.toml. Defaults to ~/.vartriage/config.toml.

        Returns
        -------
        BundleConfig
            Parsed configuration with env var overrides applied.
        """
        if path is None:
            path = Path.home() / ".vartriage" / "config.toml"

        config = cls()

        if path.exists() and tomllib is not None:
            config._apply_toml(path)

        config._apply_env_overrides()
        return config

    def _apply_toml(self, path: Path) -> None:
        """Parse TOML file and apply values to config."""
        try:
            with open(path, "rb") as f:
                data = tomllib.load(f)
        except (OSError, ValueError):
            return

        bundle_section = data.get("bundle", {})
        if "default_build" in bundle_section:
            self.default_build = str(bundle_section["default_build"])
        if "download_concurrency" in bundle_section:
            self.download_concurrency = int(bundle_section["download_concurrency"])
        if "storage_path" in bundle_section:
            self.storage_path = Path(
                os.path.expanduser(str(bundle_section["storage_path"]))
            )
        if "auto_verify" in bundle_section:
            self.auto_verify = bool(bundle_section["auto_verify"])

        proxy_section = bundle_section.get("proxy", {})
        if "http_proxy" in proxy_section:
            self.http_proxy = str(proxy_section["http_proxy"])
        if "https_proxy" in proxy_section:
            self.https_proxy = str(proxy_section["https_proxy"])

    def _apply_env_overrides(self) -> None:
        """Apply environment variable overrides."""
        env_storage = os.environ.get("VARTRIAGE_BUNDLE_STORAGE")
        if env_storage:
            self.storage_path = Path(env_storage)

        env_build = os.environ.get("VARTRIAGE_DEFAULT_BUILD")
        if env_build:
            self.default_build = env_build
