"""ACMG/AMP 2015 evidence combining rules.

Maps a set of evidence tags (with strength tiers) to a final classification.
Counts evidence at each level and checks against threshold patterns.

Combining Rules
---------------
PATHOGENIC when any of:
- ≥1 Very Strong + ≥1 Strong
- ≥2 Strong + ≥1 Supporting
- ≥1 Very Strong + ≥2 Supporting

LIKELY_PATHOGENIC when any of:
- ≥1 Very Strong + ≥1 Moderate
- ≥1 Strong + 1–2 Moderate
- ≥1 Strong + ≥2 Supporting

Everything else (including empty evidence) → VUS.
LIKELY_BENIGN and BENIGN are not produced yet. No benign tags are assigned.
"""

from __future__ import annotations

from vartriage.models.variant import (
    ACMGClassification,
    EvidenceStrength,
    EvidenceTag,
    EVIDENCE_STRENGTH_MAP,
)


def combine_evidence(tags: frozenset[EvidenceTag]) -> ACMGClassification:
    """Combine evidence tags into a final ACMG classification.

    Parameters
    ----------
    tags : frozenset[EvidenceTag]
        Evidence tags assigned to a variant.

    Returns
    -------
    ACMGClassification
        PATHOGENIC, LIKELY_PATHOGENIC, or VUS.
    """
    if not tags:
        return ACMGClassification.VUS

    counts = _count_by_strength(tags)

    if _meets_pathogenic(counts):
        return ACMGClassification.PATHOGENIC

    if _meets_likely_pathogenic(counts):
        return ACMGClassification.LIKELY_PATHOGENIC

    return ACMGClassification.VUS


def _count_by_strength(
    tags: frozenset[EvidenceTag],
) -> dict[EvidenceStrength, int]:
    """Tally tags by strength tier."""
    counts: dict[EvidenceStrength, int] = {
        EvidenceStrength.VERY_STRONG: 0,
        EvidenceStrength.STRONG: 0,
        EvidenceStrength.MODERATE: 0,
        EvidenceStrength.SUPPORTING: 0,
    }
    for tag in tags:
        strength = EVIDENCE_STRENGTH_MAP[tag]
        counts[strength] += 1
    return counts


def _meets_pathogenic(counts: dict[EvidenceStrength, int]) -> bool:
    """True if counts satisfy any Pathogenic rule."""
    vs = counts[EvidenceStrength.VERY_STRONG]
    s = counts[EvidenceStrength.STRONG]
    sup = counts[EvidenceStrength.SUPPORTING]

    if vs >= 1 and s >= 1:
        return True
    if s >= 2 and sup >= 1:
        return True
    if vs >= 1 and sup >= 2:
        return True

    return False


def _meets_likely_pathogenic(counts: dict[EvidenceStrength, int]) -> bool:
    """True if counts satisfy any Likely Pathogenic rule."""
    vs = counts[EvidenceStrength.VERY_STRONG]
    s = counts[EvidenceStrength.STRONG]
    m = counts[EvidenceStrength.MODERATE]
    sup = counts[EvidenceStrength.SUPPORTING]

    if vs >= 1 and m >= 1:
        return True
    if s >= 1 and 1 <= m <= 2:
        return True
    if s >= 1 and sup >= 2:
        return True

    return False
