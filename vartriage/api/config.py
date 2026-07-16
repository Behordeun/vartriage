"""Configuration for the API annotation backend.

Loads settings from (in priority order):
1. Explicit constructor arguments
2. Environment variables (NCBI_API_KEY, HTTP_PROXY, HTTPS_PROXY)
3. TOML config file (~/.vartriage/config.toml [api] section)
4. Hardcoded defaults

All rate limits and timeouts are validated at construction time.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional


def _load_toml_api_section(config_path: Path) -> dict[str, object]:
    """Load the [api] section from a TOML config file.

    Returns an empty dict if the file doesn't exist or has no [api] section.
    """
    if not config_path.exists():
        return {}

    if sys.version_info >= (3, 11):
        import tomllib
    else:
        try:
            import tomli as tomllib
        except ImportError:
            return {}

    try:
        with open(config_path, "rb") as f:
            data = tomllib.load(f)
        return dict(data.get("api", {}))
    except (OSError, ValueError, KeyError):
        return {}


# Mapping from nested TOML section keys to (source_key -> target_key) pairs.
# None value means merge the entire sub-dict directly.
_TOML_SECTION_KEYS: dict[str, dict[str, str] | None] = {
    "rate_limits": None,
    "timeouts": {"connect_seconds": "connect_timeout", "read_seconds": "read_timeout"},
    "proxy": {"url": "proxy_url"},
}


def _merge_toml_values(toml_values: dict[str, object]) -> dict[str, object]:
    """Merge nested TOML sections into a flat config dict."""
    merged: dict[str, object] = {}
    for key, val in toml_values.items():
        mapping = _TOML_SECTION_KEYS.get(key)
        if mapping is None and key not in _TOML_SECTION_KEYS:
            merged[key] = val
        elif not isinstance(val, dict):
            merged[key] = val
        elif mapping is None:
            merged.update(val)
        else:
            for src, dest in mapping.items():
                if src in val:
                    merged[dest] = val[src]
    return merged


def _apply_env_vars(merged: dict[str, object]) -> None:
    """Apply environment variable overrides to the merged config dict."""
    env_api_key = os.environ.get("NCBI_API_KEY")
    if env_api_key:
        merged["ncbi_api_key"] = env_api_key

    env_proxy = os.environ.get("HTTPS_PROXY") or os.environ.get("HTTP_PROXY")
    if env_proxy and "proxy_url" not in merged:
        merged["proxy_url"] = env_proxy


@dataclass(frozen=True)
class APIConfig:
    """Configuration for the API-based annotation backend.

    Parameters
    ----------
    mode
        Pipeline annotation mode. "local" uses file-based backends only.
        "api" queries remote services. "hybrid" fills gaps via API.
    genome_build
        Target genome assembly for coordinate-aware APIs.
    ncbi_api_key
        NCBI API key for higher ClinVar rate limits (10 req/sec vs 3).
        Also read from NCBI_API_KEY environment variable.
    cache_path
        SQLite cache database location.
    cache_ttl_days
        Default cache entry lifetime. -1 pins entries indefinitely.
    vep_batch_size
        Variants per VEP POST request (Ensembl max is 200).
    max_retries
        Retry attempts for transient HTTP failures.
    connect_timeout
        TCP connection timeout in seconds.
    read_timeout
        HTTP response read timeout in seconds.
    vep_rate_limit
        Ensembl VEP requests per second.
    clinvar_rate_limit
        NCBI ClinVar requests per second (with API key).
    cadd_rate_limit
        CADD API requests per second.
    spliceai_rate_limit
        SpliceAI Lookup requests per second.
    vep_daily_limit
        Maximum VEP requests per UTC day. None for unlimited.
    proxy_url
        Explicit HTTP/HTTPS proxy URL. Also reads HTTP_PROXY/HTTPS_PROXY env.
    preferred_frequency_source
        Which gnomAD frequency to extract from VEP colocated_variants.

    Raises
    ------
    ValueError
        If any parameter is outside its valid range.
    """

    mode: Literal["local", "api", "hybrid"] = "local"
    genome_build: Literal["grch37", "grch38"] = "grch38"
    ncbi_api_key: Optional[str] = None
    cache_path: Path = field(
        default_factory=lambda: Path.home() / ".vartriage" / "api_cache.db"
    )
    cache_ttl_days: int = 7
    vep_batch_size: int = 200
    max_retries: int = 3
    connect_timeout: float = 10.0
    read_timeout: float = 30.0
    vep_rate_limit: float = 15.0
    clinvar_rate_limit: float = 10.0
    cadd_rate_limit: float = 2.0
    spliceai_rate_limit: float = 0.08
    vep_daily_limit: Optional[int] = 55_000
    proxy_url: Optional[str] = None
    preferred_frequency_source: Literal["gnomad_exome", "gnomad_genome"] = (
        "gnomad_exome"
    )

    def __post_init__(self) -> None:
        if self.vep_batch_size < 1 or self.vep_batch_size > 200:
            raise ValueError(f"vep_batch_size must be 1-200, got {self.vep_batch_size}")
        if self.max_retries < 0 or self.max_retries > 10:
            raise ValueError(f"max_retries must be 0-10, got {self.max_retries}")
        if self.connect_timeout <= 0:
            raise ValueError("connect_timeout must be positive")
        if self.read_timeout <= 0:
            raise ValueError("read_timeout must be positive")
        if self.vep_rate_limit <= 0:
            raise ValueError("vep_rate_limit must be positive")
        if self.clinvar_rate_limit <= 0:
            raise ValueError("clinvar_rate_limit must be positive")
        if self.cadd_rate_limit <= 0:
            raise ValueError("cadd_rate_limit must be positive")
        if self.spliceai_rate_limit <= 0:
            raise ValueError("spliceai_rate_limit must be positive")
        if self.cache_ttl_days < -1 or self.cache_ttl_days == 0:
            raise ValueError(
                "cache_ttl_days must be positive or -1 (pinned), got "
                f"{self.cache_ttl_days}"
            )

    @classmethod
    def load(
        cls,
        config_path: Path | None = None,
        **overrides: object,
    ) -> "APIConfig":
        """Build APIConfig from TOML file + env vars + explicit overrides.

        Priority (highest wins): overrides > env vars > TOML > defaults.

        Parameters
        ----------
        config_path
            Path to config.toml. Defaults to ~/.vartriage/config.toml.
        **overrides
            Explicit field values that take highest priority.
        """
        toml_path = config_path or Path.home() / ".vartriage" / "config.toml"
        toml_values = _load_toml_api_section(toml_path)

        merged = _merge_toml_values(toml_values)
        _apply_env_vars(merged)

        # Explicit overrides take final priority
        merged.update({k: v for k, v in overrides.items() if v is not None})

        # Convert cache_path string to Path if needed
        if "cache_path" in merged and isinstance(merged["cache_path"], str):
            merged["cache_path"] = Path(merged["cache_path"]).expanduser()

        # Filter to only valid APIConfig fields
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in merged.items() if k in valid_fields}

        return cls(**filtered)  # type: ignore[arg-type]
