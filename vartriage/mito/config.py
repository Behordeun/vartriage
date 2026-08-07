"""Configuration for the mitochondrial variant analysis pipeline.

Provides MitoConfig as a frozen dataclass with startup validation,
consistent with the pattern used by other pipeline configs in
vartriage.models.config.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MitoConfig:
    """Configuration for the mitochondrial analysis sub-pipeline.

    Parameters
    ----------
    enabled : bool
        Whether mitochondrial analysis is active. When False, chrM
        variants pass through without mtDNA-specific processing.
        Defaults to True (auto-enabled when chrM variants are present).
    min_heteroplasmy : float
        Minimum heteroplasmy percentage for reporting. Variants below
        this threshold are filtered as sub-threshold noise. Must be
        in [0.0, 100.0]. Default is 1.0%.
    gene_map_path : Path or None
        Custom path to mt_gene_map.tsv. When None, uses the
        package-bundled default.
    mitomap_path : Path or None
        Custom path to mitomap_pathogenic.tsv. When None, uses the
        package-bundled default.
    helixmtdb_path : Path or None
        Custom path to helixmtdb_frequency.tsv. When None, uses the
        package-bundled default.

    Raises
    ------
    ValueError
        If min_heteroplasmy is outside [0.0, 100.0].
    """

    enabled: bool = True
    min_heteroplasmy: float = 1.0
    gene_map_path: Path | None = None
    mitomap_path: Path | None = None
    helixmtdb_path: Path | None = None

    def __post_init__(self) -> None:
        if not (0.0 <= self.min_heteroplasmy <= 100.0):
            raise ValueError(
                f"min_heteroplasmy must be between 0.0 and 100.0, "
                f"got {self.min_heteroplasmy}"
            )
