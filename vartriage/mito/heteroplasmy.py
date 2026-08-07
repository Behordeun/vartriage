"""Heteroplasmy extraction and level classification for mtDNA variants.

Mitochondrial DNA exists in hundreds to thousands of copies per cell.
Variants can be present in a fraction of copies (heteroplasmy) rather
than the binary het/hom states of nuclear diploid variants. This module
extracts the alternate allele fraction from VCF FORMAT fields and
classifies it into clinically meaningful categories.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Literal

logger = logging.getLogger(__name__)

HeteroplasmyCategory = Literal[
    "homoplasmic",
    "high",
    "moderate",
    "low",
    "sub_threshold",
]

# Threshold boundaries (applied to percentage 0-100)
_HOMOPLASMIC_THRESHOLD: float = 95.0
_HIGH_THRESHOLD: float = 60.0
_MODERATE_THRESHOLD: float = 20.0
_LOW_THRESHOLD: float = 1.0


@dataclass(frozen=True, slots=True)
class HeteroplasmyLevel:
    """Heteroplasmy measurement for a mitochondrial variant.

    Parameters
    ----------
    fraction
        Alternate allele fraction (0.0 to 1.0).
    percentage
        Same value expressed as a percentage (0.0 to 100.0).
    category
        Clinical classification bucket based on the percentage.
    depth
        Total read depth at this position (sum of all allele depths).
    """

    fraction: float
    percentage: float
    category: HeteroplasmyCategory
    depth: int


def classify_level(percentage: float) -> HeteroplasmyCategory:
    """Classify a heteroplasmy percentage into a clinical category.

    Thresholds per the spec:
    - homoplasmic: >= 95%
    - high: 60% to <95%
    - moderate: 20% to <60%
    - low: 1% to <20%
    - sub_threshold: < 1% (likely artifact or sequencing noise)

    Parameters
    ----------
    percentage
        Heteroplasmy as a percentage (0.0-100.0).

    Returns
    -------
    HeteroplasmyCategory
        The classification bucket.
    """
    if percentage >= _HOMOPLASMIC_THRESHOLD:
        return "homoplasmic"
    if percentage >= _HIGH_THRESHOLD:
        return "high"
    if percentage >= _MODERATE_THRESHOLD:
        return "moderate"
    if percentage >= _LOW_THRESHOLD:
        return "low"
    return "sub_threshold"


def extract_heteroplasmy(info: dict[str, Any]) -> HeteroplasmyLevel | None:
    """Extract heteroplasmy level from VCF variant info dict.

    Attempts two extraction strategies in order:
    1. AD field (allele depths): computes ALT_depth / total_depth
    2. AF field (allele fraction): uses the value directly

    The info dict is expected to carry FORMAT-level data under the keys
    used by the VCFParser's sample extraction (AD as a tuple/list of ints,
    AF as a float or list of floats).

    Parameters
    ----------
    info
        Variant info dict, potentially containing AD or AF from
        FORMAT fields (populated by VCFParser when extract_samples=True).

    Returns
    -------
    HeteroplasmyLevel or None
        Extracted heteroplasmy, or None if neither AD nor AF is available
        or the data is malformed.
    """
    # Strategy 1: AD field (Mutect2 mitochondrial mode output)
    ad = info.get("AD")
    if ad is not None:
        result = _from_allele_depths(ad)
        if result is not None:
            return result

    # Strategy 2: AF field (caller-provided allele fraction)
    af = info.get("AF")
    if af is not None:
        result = _from_allele_fraction(af)
        if result is not None:
            return result

    return None


def _from_allele_depths(ad: Any) -> HeteroplasmyLevel | None:
    """Compute heteroplasmy from allele depth values.

    Expects AD as a sequence of at least 2 integers: [REF_depth, ALT_depth].
    For multi-allelic sites (already split upstream), only index 1 is used.
    """
    if not hasattr(ad, "__len__") or len(ad) < 2:
        return None

    try:
        ref_depth = int(ad[0])
        alt_depth = int(ad[1])
    except (TypeError, ValueError):
        return None

    if ref_depth < 0 or alt_depth < 0:
        return None

    total = ref_depth + alt_depth
    if total <= 0:
        return None

    fraction = alt_depth / total
    percentage = fraction * 100.0
    category = classify_level(percentage)

    return HeteroplasmyLevel(
        fraction=fraction,
        percentage=percentage,
        category=category,
        depth=total,
    )


def _from_allele_fraction(af: Any) -> HeteroplasmyLevel | None:
    """Build heteroplasmy from a pre-computed allele fraction value.

    AF may be a single float or a list (take first element for the ALT allele).
    Depth is unknown in this path, so we report 0.
    """
    try:
        value = float(af[0]) if hasattr(af, "__len__") else float(af)
    except (TypeError, ValueError, IndexError):
        return None

    if value < 0.0 or value > 1.0:
        return None

    percentage = value * 100.0
    category = classify_level(percentage)

    return HeteroplasmyLevel(
        fraction=value,
        percentage=percentage,
        category=category,
        depth=0,
    )
