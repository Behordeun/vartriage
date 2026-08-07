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


def extract_heteroplasmy(
    info: dict[str, Any], sample_name: str | None = None
) -> HeteroplasmyLevel | None:
    """Extract heteroplasmy level from VCF variant info dict.

    Searches for AD/AF in three locations (first match wins):
    1. Direct keys in info dict (INFO-level AD/AF from some callers)
    2. Per-sample FORMAT fields in _pysam_samples[sample_name]
    3. Per-sample FORMAT fields in _pysam_samples (first sample if
       sample_name is None)

    Parameters
    ----------
    info
        Variant info dict, potentially containing AD or AF directly,
        or _pysam_samples with per-sample FORMAT data from VCFParser.
    sample_name
        Sample to extract FORMAT data from. When None, uses the first
        available sample in _pysam_samples.

    Returns
    -------
    HeteroplasmyLevel or None
        Extracted heteroplasmy, or None if neither AD nor AF is available
        or the data is malformed.
    """
    # Strategy 1: direct INFO-level AD/AF (some callers put these in INFO)
    ad = info.get("AD")
    if ad is not None:
        result = _from_allele_depths(ad)
        if result is not None:
            return result

    af = info.get("AF")
    if af is not None:
        result = _from_allele_fraction(af)
        if result is not None:
            return result

    # Strategy 2: look in _pysam_samples FORMAT fields
    sample_data = info.get("_pysam_samples")
    if sample_data is not None:
        sample_entry = _get_sample_entry(sample_data, sample_name)
        if sample_entry is not None:
            ad = sample_entry.get("AD")
            if ad is not None:
                result = _from_allele_depths(ad)
                if result is not None:
                    return result
            af = sample_entry.get("AF")
            if af is not None:
                result = _from_allele_fraction(af)
                if result is not None:
                    return result

    return None


def _get_sample_entry(
    sample_data: dict[str, Any], sample_name: str | None
) -> dict[str, Any] | None:
    """Get the sample entry dict from _pysam_samples."""
    if not sample_data:
        return None
    if sample_name is not None:
        return sample_data.get(sample_name)
    # Fall back to first sample
    for entry in sample_data.values():
        if isinstance(entry, dict):
            return entry
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
