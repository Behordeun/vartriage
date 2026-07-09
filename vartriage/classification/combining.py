"""ACMG/AMP 2015 evidence combining rules for variant classification.

This module implements the combining logic that maps a set of evidence tags
(each carrying a strength tier) to a final ACMG classification. The rules
follow the 2015 ACMG/AMP guidelines for interpreting sequence variants.

The combining function counts evidence at each strength level and checks
whether any of the defined rule thresholds are met, returning the most
severe classification that applies.
"""

from __future__ import annotations

from vartriage.models.variant import (
    ACMGClassification,
    EvidenceStrength,
    EvidenceTag,
    EVIDENCE_STRENGTH_MAP,
)


def combine_evidence(tags: frozenset[EvidenceTag]) -> ACMGClassification:
    """Combine ACMG evidence tags into a final classification.

    Applies ACMG/AMP 2015 combining rules by counting evidence at each
    strength tier and checking whether pathogenic or likely pathogenic
    thresholds are met.

    Parameters
    ----------
    tags : frozenset[EvidenceTag]
        Set of evidence tags assigned to a variant. Each tag maps to a
        strength tier via ``EVIDENCE_STRENGTH_MAP``.

    Returns
    -------
    ACMGClassification
        The final classification: PATHOGENIC, LIKELY_PATHOGENIC, or VUS.
        Likely_Benign and Benign are not produced in v1 since benign
        evidence tags are not assigned by the current rule set.
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
    """Count how many tags fall into each strength tier.

    Parameters
    ----------
    tags : frozenset[EvidenceTag]
        Evidence tags to categorize.

    Returns
    -------
    dict[EvidenceStrength, int]
        Counts per strength tier, defaulting to zero for absent tiers.
    """
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
    """Check whether evidence meets any Pathogenic combining rule.

    Pathogenic rules (ACMG/AMP 2015):
    - >=1 Very Strong AND >=1 Strong
    - >=2 Strong AND >=1 Supporting
    - >=1 Very Strong AND >=2 Supporting

    Parameters
    ----------
    counts : dict[EvidenceStrength, int]
        Evidence counts by strength tier.

    Returns
    -------
    bool
        True if any pathogenic rule is satisfied.
    """
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
    """Check whether evidence meets any Likely Pathogenic combining rule.

    Likely Pathogenic rules (ACMG/AMP 2015):
    - 1 Very Strong AND 1 Moderate
    - 1 Strong AND 1-2 Moderate
    - 1 Strong AND >=2 Supporting

    Parameters
    ----------
    counts : dict[EvidenceStrength, int]
        Evidence counts by strength tier.

    Returns
    -------
    bool
        True if any likely pathogenic rule is satisfied.
    """
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
