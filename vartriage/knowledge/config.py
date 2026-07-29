"""Configuration for the gene knowledge base."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

VALID_INHERITANCE_MODES = frozenset({"AD", "AR", "XL", "XLD", "XLR", "MT"})


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
    inheritance_mode : str | None
        Filter variants to genes matching this inheritance pattern
        (AD, AR, XL, XLD, XLR, MT). None disables the filter.
    flag_actionable : bool
        When True, gene_context.is_actionable is populated from
        ClinGen actionability curations. When False, actionability
        annotation still occurs but no filtering is applied.

    Raises
    ------
    ValueError
        If any HPO term doesn't match the HP:NNNNNNN format.
    ValueError
        If inheritance_mode is not a recognized mode.
    """

    data_dir: Optional[Path] = field(default=None)
    hpo_terms: frozenset[str] = field(default_factory=frozenset)
    inheritance_mode: Optional[str] = field(default=None)
    flag_actionable: bool = False

    def __post_init__(self) -> None:
        for term in self.hpo_terms:
            if not term.startswith("HP:") or len(term) != 10:
                raise ValueError(
                    f"HPO terms must match HP:NNNNNNN format, got '{term}'"
                )
        if (
            self.inheritance_mode is not None
            and self.inheritance_mode not in VALID_INHERITANCE_MODES
        ):
            raise ValueError(
                f"inheritance_mode must be one of {sorted(VALID_INHERITANCE_MODES)}, "
                f"got '{self.inheritance_mode}'"
            )

    @property
    def resolved_data_dir(self) -> Path:
        """Return the effective data directory (explicit or bundled default)."""
        if self.data_dir is not None:
            return self.data_dir
        return _default_data_dir()
