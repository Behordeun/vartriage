"""ACMG/AMP evidence combining via the ClinGen SVI point system.

A set of evidence tags is mapped to a final classification by the
point-based framework of Tavtigian et al. (2018, 2020): each criterion
contributes a signed point value by strength, the points are summed, and
one threshold ladder decides the tier. This weighs opposing evidence by
arithmetic rather than resolving it on the bare co-occurrence of a
pathogenic and a benign criterion.

The point assignment and thresholds live in
``vartriage.classification.points``. BA1 remains a stand-alone Benign
override per ACMG 2015 Table 5.

``has_conflicting_evidence`` reports whether both a pathogenic and a benign
criterion are present. It is a derived annotation for reporting; it does
not by itself force a classification, because the point sum already
accounts for evidence on both sides.
"""

from __future__ import annotations

from vartriage.classification.points import classify_by_points
from vartriage.models.variant import (
    ACMGClassification,
    EvidenceTag,
)

# Tags carrying benign-direction evidence, used only by the reporting flag.
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


def combine_evidence(
    tags: frozenset[EvidenceTag],
) -> ACMGClassification:
    """Combine evidence tags into a final ACMG classification.

    Delegates to the SVI point engine: the signed points for the tag set
    are summed and mapped through the threshold ladder, with BA1 as a
    stand-alone Benign override.

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
    return classify_by_points(tags)


def has_conflicting_evidence(tags: frozenset[EvidenceTag]) -> bool:
    """Report whether both pathogenic and benign tags are present.

    A derived flag for reporting. The point sum, not this flag, decides
    the classification, so conflicting evidence resolves to the net tier
    rather than being forced to VUS.
    """
    pathogenic_tags = tags - _BENIGN_TAGS
    benign_tags = tags & _BENIGN_TAGS
    return len(pathogenic_tags) > 0 and len(benign_tags) > 0
