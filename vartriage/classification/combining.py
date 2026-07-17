"""ACMG/AMP 2015 evidence combining rules.

Maps a set of evidence tags (with strength tiers) to a final classification.
Separates pathogenic and benign evidence, checks for conflicts, then
applies threshold patterns from ACMG Table 5.

Pathogenic Combining Rules:
- PATHOGENIC: 1 VS + 1 S, or 2 S + 1 Sup, or 1 VS + 2 Sup
- LIKELY_PATHOGENIC: 1 VS + 1 M, or 1 S + 1-2 M, or 1 S + 2 Sup

Benign Combining Rules:
- BENIGN: 1 BA (standalone), or 2 BS
- LIKELY_BENIGN: 1 BS + 1 BP

Note: 2 BP alone is insufficient for classification change per ACMG
Table 5 (returns VUS). This is intentional.

Conflicting Evidence:
- Pathogenic + Benign evidence present = VUS with conflict flag
"""

from __future__ import annotations

from vartriage.models.variant import (EVIDENCE_STRENGTH_MAP,
                                      ACMGClassification, EvidenceStrength,
                                      EvidenceTag)

# Benign tags are those with STANDALONE, STRONG (benign), or SUPPORTING (benign) strength
_BENIGN_TAGS: frozenset[EvidenceTag] = frozenset({
    EvidenceTag.BA1, EvidenceTag.BS1, EvidenceTag.BS2,
    EvidenceTag.BP4, EvidenceTag.BP7,
})


def combine_evidence(
    tags: frozenset[EvidenceTag],
) -> ACMGClassification:
    """Combine evidence tags into a final ACMG classification.

    Separates pathogenic and benign evidence. When both are present,
    the result is VUS (conflicting evidence). Otherwise applies the
    appropriate combining rules.

    Parameters
    ----------
    tags : frozenset[EvidenceTag]
        Evidence tags assigned to a variant.

    Returns
    -------
    ACMGClassification
        Final classification: PATHOGENIC, LIKELY_PATHOGENIC, VUS,
        LIKELY_BENIGN, or BENIGN.
    """
    if not tags:
        return ACMGClassification.VUS

    pathogenic_tags = tags - _BENIGN_TAGS
    benign_tags = tags & _BENIGN_TAGS

    has_pathogenic = len(pathogenic_tags) > 0
    has_benign = len(benign_tags) > 0

    # Conflicting evidence: both pathogenic and benign present
    if has_pathogenic and has_benign:
        return ACMGClassification.VUS

    # Pure benign evidence
    if has_benign and not has_pathogenic:
        return _classify_benign(benign_tags)

    # Pure pathogenic evidence
    counts = _count_pathogenic_strengths(pathogenic_tags)

    if _meets_pathogenic(counts):
        return ACMGClassification.PATHOGENIC

    if _meets_likely_pathogenic(counts):
        return ACMGClassification.LIKELY_PATHOGENIC

    return ACMGClassification.VUS


def has_conflicting_evidence(tags: frozenset[EvidenceTag]) -> bool:
    """Check if both pathogenic and benign tags are present."""
    pathogenic_tags = tags - _BENIGN_TAGS
    benign_tags = tags & _BENIGN_TAGS
    return len(pathogenic_tags) > 0 and len(benign_tags) > 0


def _classify_benign(benign_tags: frozenset[EvidenceTag]) -> ACMGClassification:
    """Apply benign combining rules."""
    # BA1 standalone = Benign
    if EvidenceTag.BA1 in benign_tags:
        return ACMGClassification.BENIGN

    bs_count = sum(
        1 for t in benign_tags
        if EVIDENCE_STRENGTH_MAP[t] == EvidenceStrength.STRONG
    )
    bp_count = sum(
        1 for t in benign_tags
        if EVIDENCE_STRENGTH_MAP[t] == EvidenceStrength.SUPPORTING
    )

    # 2 BS = Benign
    if bs_count >= 2:
        return ACMGClassification.BENIGN

    # 1 BS + 1 BP = Likely Benign
    if bs_count >= 1 and bp_count >= 1:
        return ACMGClassification.LIKELY_BENIGN

    # 2 BP alone = not sufficient for classification change
    return ACMGClassification.VUS


def _count_pathogenic_strengths(
    tags: frozenset[EvidenceTag],
) -> dict[EvidenceStrength, int]:
    """Tally pathogenic tags by strength tier."""
    counts: dict[EvidenceStrength, int] = {
        EvidenceStrength.VERY_STRONG: 0,
        EvidenceStrength.STRONG: 0,
        EvidenceStrength.MODERATE: 0,
        EvidenceStrength.SUPPORTING: 0,
    }
    for tag in tags:
        strength = EVIDENCE_STRENGTH_MAP[tag]
        if strength in counts:
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
