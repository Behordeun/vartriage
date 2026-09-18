"""ClinGen SVI Bayesian point system for ACMG/AMP evidence combining.

The point-based framework of Tavtigian et al. (2018, 2020), adopted by the
ClinGen Sequence Variant Interpretation working group, assigns a signed
integer to each evidence criterion by its strength and sums them. A single
threshold ladder then maps the total to a classification. This replaces the
combinatorial rule table with arithmetic that weighs opposing evidence
directly rather than vetoing on the mere co-occurrence of both signs.

Point values (pathogenic positive, benign negative):

    Supporting    1
    Moderate      2
    Strong        4
    Very Strong   8

Tier thresholds on the summed total:

    Pathogenic         >= 10
    Likely Pathogenic    6 .. 9
    VUS                 -5 .. 5
    Likely Benign       -6 .. -1
    Benign             <= -7

BA1 (stand-alone benign) short-circuits to Benign regardless of the total,
matching ACMG 2015 Table 5.

This module is a pure, side-effect-free scoring function. It is shipped
alongside the existing combining rules; wiring the classifier onto it is a
separate change.
"""

from __future__ import annotations

from vartriage.models.variant import (
    EVIDENCE_STRENGTH_MAP,
    ACMGClassification,
    EvidenceStrength,
    EvidenceTag,
)

_PATHOGENIC_POINTS: dict[EvidenceStrength, int] = {
    EvidenceStrength.VERY_STRONG: 8,
    EvidenceStrength.STRONG: 4,
    EvidenceStrength.MODERATE: 2,
    EvidenceStrength.SUPPORTING: 1,
}

_BENIGN_TAGS: frozenset[EvidenceTag] = frozenset(
    {
        EvidenceTag.BA1,
        EvidenceTag.BS1,
        EvidenceTag.BS2,
        EvidenceTag.BP4,
        EvidenceTag.BP4_MODERATE,
        EvidenceTag.BP7,
    }
)

_PATHOGENIC_THRESHOLD = 10
_LIKELY_PATHOGENIC_THRESHOLD = 6
_LIKELY_BENIGN_THRESHOLD = -1
_BENIGN_THRESHOLD = -7


def score_points(tags: frozenset[EvidenceTag]) -> int:
    """Sum the signed points for a set of evidence tags.

    Pathogenic tags contribute positive points by strength; benign tags
    contribute the same magnitude negated. BA1 carries no point value (it
    is a stand-alone override handled by :func:`classify_by_points`) and
    contributes zero here.
    """
    total = 0
    for tag in tags:
        if tag is EvidenceTag.BA1:
            continue
        magnitude = _PATHOGENIC_POINTS[EVIDENCE_STRENGTH_MAP[tag]]
        if tag in _BENIGN_TAGS:
            total -= magnitude
        else:
            total += magnitude
    return total


def classify_by_points(tags: frozenset[EvidenceTag]) -> ACMGClassification:
    """Classify a variant by summed evidence points.

    BA1 short-circuits to Benign. Otherwise the signed total is mapped
    through the SVI threshold ladder.
    """
    if EvidenceTag.BA1 in tags:
        return ACMGClassification.BENIGN

    total = score_points(tags)

    if total >= _PATHOGENIC_THRESHOLD:
        return ACMGClassification.PATHOGENIC
    if total >= _LIKELY_PATHOGENIC_THRESHOLD:
        return ACMGClassification.LIKELY_PATHOGENIC
    if total <= _BENIGN_THRESHOLD:
        return ACMGClassification.BENIGN
    if total <= _LIKELY_BENIGN_THRESHOLD:
        return ACMGClassification.LIKELY_BENIGN
    return ACMGClassification.VUS
