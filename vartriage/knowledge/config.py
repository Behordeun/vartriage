"""Configuration for the gene knowledge base."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


def _default_data_dir() -> Path:
    """Resolve the bundled knowledge data directory."""
    return Path(__file__).resolve().parent.parent / "data" / "knowledge"


@dataclass(frozen=True)
class KnowledgeBaseConfig:
    """Configuration for gene-disease linkage knowledge base.

    Parameters
    ----------
    data_dir : Path | None
        Directory containing pre-processed TSV knowledge files.
        When None, uses the bundled package data at
        ``vartriage/data/knowledge/``.
    hpo_terms : frozenset[str]
        Patient HPO term IDs for phenotype-driven prioritization.
        Empty set disables phenotype boosting.

    Raises
    ------
    ValueError
        If any HPO term doesn't match the HP:NNNNNNN format.
    """

    data_dir: Optional[Path] = field(default=None)
    hpo_terms: frozenset[str] = field(default_factory=frozenset)

    def __post_init__(self) -> None:
        for term in self.hpo_terms:
            if not term.startswith("HP:") or len(term) != 10:
                raise ValueError(
                    f"HPO terms must match HP:NNNNNNN format, got '{term}'"
                )

    @property
    def resolved_data_dir(self) -> Path:
        """Return the effective data directory (explicit or bundled default)."""
        if self.data_dir is not None:
            return self.data_dir
        return _default_data_dir()
