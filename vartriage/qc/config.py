"""QC configuration dataclass with startup validation."""

from __future__ import annotations

import sys
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class QCConfig:
    """Configuration for VCF quality control checks.

    Parameters
    ----------
    assay_type : Literal["wgs", "wes", "panel"]
        Sequencing assay type governing expected metric ranges.
        Default is "wes".
    strict : bool
        When True, any FAIL-level metric halts the pipeline with exit
        code 3 before annotation begins. Default is False.
    skip : bool
        When True, QC is bypassed entirely (for pre-validated files).
        Default is False.
    expected_ti_tv : tuple[float, float] | None
        Override warn-level Ti/Tv range (min, max). When None, uses the
        assay-type default.
    expected_het_hom : tuple[float, float] | None
        Override warn-level het/hom range (min, max). When None, uses the
        assay-type default.
    sample_id : str | None
        Sample name for het/hom extraction from multi-sample VCFs. When
        None and the VCF has exactly one sample, that sample is used
        automatically.

    Raises
    ------
    ValueError
        If expected_ti_tv or expected_het_hom has min >= max.
    ValueError
        If assay_type is not one of the supported values.
    """

    assay_type: Literal["wgs", "wes", "panel"] = "wes"
    strict: bool = False
    skip: bool = False
    expected_ti_tv: tuple[float, float] | None = None
    expected_het_hom: tuple[float, float] | None = None
    sample_id: str | None = None

    _VALID_ASSAY_TYPES: frozenset[str] = field(
        default=frozenset({"wgs", "wes", "panel"}),
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        if self.assay_type not in self._VALID_ASSAY_TYPES:
            raise ValueError(
                f"assay_type must be one of {sorted(self._VALID_ASSAY_TYPES)}, "
                f"got '{self.assay_type}'"
            )
        if self.expected_ti_tv is not None:
            lo, hi = self.expected_ti_tv
            if lo >= hi:
                raise ValueError(
                    f"expected_ti_tv min must be less than max, got ({lo}, {hi})"
                )
        if self.expected_het_hom is not None:
            lo, hi = self.expected_het_hom
            if lo >= hi:
                raise ValueError(
                    f"expected_het_hom min must be less than max, got ({lo}, {hi})"
                )

    def merged_with_toml(self, config_path: Path | None = None) -> QCConfig:
        """Return a copy with TOML [qc] values filling unset fields.

        CLI-provided values (already on this instance) take precedence over
        TOML values. Only the expected_titv and expected_het_hom warn ranges
        are read from TOML, and only when not already set on the CLI.

        Parameters
        ----------
        config_path : Path | None
            Path to a TOML config file. Defaults to
            ~/.vartriage/config.toml when None.

        Returns
        -------
        QCConfig
            A new config with TOML defaults applied where CLI values
            were absent.
        """
        toml_values = _load_toml_qc_section(config_path)
        if not toml_values:
            return self

        updates: dict[str, object] = {}

        if self.expected_ti_tv is None and "expected_titv" in toml_values:
            updates["expected_ti_tv"] = _parse_toml_range(toml_values["expected_titv"])
        if self.expected_het_hom is None and "expected_het_hom" in toml_values:
            updates["expected_het_hom"] = _parse_toml_range(
                toml_values["expected_het_hom"]
            )

        if not updates:
            return self

        return replace(self, **updates)  # type: ignore[arg-type]


def _load_toml_qc_section(config_path: Path | None) -> dict[str, object]:
    """Load the [qc] section from a TOML config file.

    Returns an empty dict when the file is missing or has no [qc] section.
    """
    path = config_path or Path.home() / ".vartriage" / "config.toml"
    if not path.exists():
        return {}

    if sys.version_info >= (3, 11):
        import tomllib
    else:
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ImportError:
            return {}

    try:
        with open(path, "rb") as f:
            data = tomllib.load(f)
        return dict(data.get("qc", {}))
    except (OSError, ValueError, KeyError):
        return {}


def _parse_toml_range(value: object) -> tuple[float, float]:
    """Parse a TOML range value (list or tuple of two numbers)."""
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(
            f"TOML range must be a two-element array [min, max], got {value!r}"
        )
    return (float(value[0]), float(value[1]))
