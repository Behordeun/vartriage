"""Assay-specific QC threshold definitions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class AssayType(Enum):
    """Supported sequencing assay types."""

    WGS = "wgs"
    WES = "wes"
    PANEL = "panel"


@dataclass(frozen=True)
class AssayThresholds:
    """Expected metric ranges for a given assay type.

    Warn-level ranges represent borderline values. Fail-level ranges
    represent critically abnormal values indicating a likely sample or
    pipeline problem.
    """

    ti_tv_warn: tuple[float, float] = (1.8, 2.5)
    ti_tv_fail: tuple[float, float] = (1.5, 3.0)
    het_hom_warn: tuple[float, float] = (1.0, 3.0)
    het_hom_fail: tuple[float, float] = (0.5, 5.0)
    variant_count_min: int = 0
    variant_count_max: int = 999_999_999
    ins_del_warn: tuple[float, float] = (0.5, 1.5)
    ins_del_fail: tuple[float, float] = (0.3, 2.0)


WGS_THRESHOLDS = AssayThresholds(
    ti_tv_warn=(1.8, 2.5),
    ti_tv_fail=(1.5, 3.0),
    het_hom_warn=(1.0, 3.0),
    het_hom_fail=(0.5, 5.0),
    variant_count_min=3_500_000,
    variant_count_max=6_000_000,
    ins_del_warn=(0.5, 1.5),
    ins_del_fail=(0.3, 2.0),
)

WES_THRESHOLDS = AssayThresholds(
    ti_tv_warn=(1.8, 2.5),
    ti_tv_fail=(1.5, 3.0),
    het_hom_warn=(1.0, 3.0),
    het_hom_fail=(0.5, 5.0),
    variant_count_min=20_000,
    variant_count_max=150_000,
    ins_del_warn=(0.5, 1.5),
    ins_del_fail=(0.3, 2.0),
)

PANEL_THRESHOLDS = AssayThresholds(
    ti_tv_warn=(1.5, 3.0),
    ti_tv_fail=(1.2, 3.5),
    het_hom_warn=(0.8, 4.0),
    het_hom_fail=(0.3, 6.0),
    variant_count_min=0,
    variant_count_max=999_999_999,
    ins_del_warn=(0.3, 2.0),
    ins_del_fail=(0.2, 3.0),
)


def get_thresholds(assay_type: str | AssayType) -> AssayThresholds:
    """Resolve thresholds for a given assay type string or enum.

    Parameters
    ----------
    assay_type : str | AssayType
        One of "wgs", "wes", "panel" (case-insensitive string) or an
        AssayType enum member.

    Returns
    -------
    AssayThresholds
        Frozen threshold configuration.

    Raises
    ------
    ValueError
        If assay_type is not recognized.
    """
    if isinstance(assay_type, AssayType):
        key = assay_type
    else:
        try:
            key = AssayType(assay_type.lower())
        except ValueError as exc:
            raise ValueError(
                f"Unknown assay type '{assay_type}'. Supported: wgs, wes, panel"
            ) from exc

    mapping = {
        AssayType.WGS: WGS_THRESHOLDS,
        AssayType.WES: WES_THRESHOLDS,
        AssayType.PANEL: PANEL_THRESHOLDS,
    }
    return mapping[key]
