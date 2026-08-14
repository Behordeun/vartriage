"""Configuration for the remote tabix score backend."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

_DEFAULT_CACHE_PATH = Path("~/.vartriage/remote_cache.db")


@dataclass(frozen=True)
class RemoteTabixConfig:
    """Configuration for remote tabix score lookups.

    Parameters
    ----------
    cadd_remote_url : str | None
        URL or named preset for remote CADD score queries. When None,
        remote CADD is disabled. Accepts a full URL or a preset name
        (e.g., "cadd-v1.7-grch38").
    gnomad_remote_url : str | None
        URL or named preset for remote gnomAD frequency queries. When
        None, remote gnomAD is disabled. Supports {chrom} placeholder
        for per-chromosome files.
    cache_ttl_days : int
        Score cache TTL in days. -1 pins entries indefinitely (clinical
        reproducibility mode). Must be >= -1.
    cache_path : Path
        Path to the SQLite cache database file.
    connect_timeout : float
        TCP connect timeout in seconds for the initial tabix connection.
        Must be > 0.
    read_timeout : float
        Per-query read timeout in seconds. Must be > 0.
    max_retries : int
        Maximum retry attempts for transient failures (5xx, timeout).
        Must be >= 0.
    batch_window_bp : int
        When variants cluster within this window (bp), they are grouped
        into a single range query to reduce HTTP round-trips.

    Raises
    ------
    ValueError
        If any parameter is outside its valid range.
    """

    cadd_remote_url: str | None = None
    gnomad_remote_url: str | None = None
    cache_ttl_days: int = 30
    cache_path: Path = field(default=_DEFAULT_CACHE_PATH)
    connect_timeout: float = 10.0
    read_timeout: float = 30.0
    max_retries: int = 3
    batch_window_bp: int = 10_000

    def __post_init__(self) -> None:
        if self.cache_ttl_days < -1:
            raise ValueError(f"cache_ttl_days must be >= -1, got {self.cache_ttl_days}")
        if self.connect_timeout <= 0:
            raise ValueError(f"connect_timeout must be > 0, got {self.connect_timeout}")
        if self.read_timeout <= 0:
            raise ValueError(f"read_timeout must be > 0, got {self.read_timeout}")
        if self.max_retries < 0:
            raise ValueError(f"max_retries must be >= 0, got {self.max_retries}")
        if self.batch_window_bp < 1:
            raise ValueError(
                f"batch_window_bp must be >= 1, got {self.batch_window_bp}"
            )

    @property
    def is_cadd_active(self) -> bool:
        """True when remote CADD scoring is configured."""
        return self.cadd_remote_url is not None

    @property
    def is_gnomad_active(self) -> bool:
        """True when remote gnomAD frequency lookup is configured."""
        return self.gnomad_remote_url is not None

    @property
    def has_any_remote(self) -> bool:
        """True when at least one remote backend is configured."""
        return self.is_cadd_active or self.is_gnomad_active
