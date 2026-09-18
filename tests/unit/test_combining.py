"""Combining rules under the ClinGen SVI point system.

combine_evidence sums signed points (Supporting 1, Moderate 2, Strong 4,
Very Strong 8; benign negated) and maps the total: Pathogenic >= 10,
Likely Pathogenic 6..9, Likely Benign -6..-1, Benign <= -7, VUS between,
with BA1 as a stand-alone Benign override.
"""

from __future__ import annotations

from vartriage.classification.combining import (
    combine_evidence,
    has_conflicting_evidence,
)
from vartriage.models.variant import ACMGClassification, EvidenceTag


class TestCombineEvidenceVUS:
    def test_empty_tag_set_returns_vus(self) -> None:
        assert combine_evidence(frozenset()) == ACMGClassification.VUS

    def test_single_supporting_tag_returns_vus(self) -> None:
        # 1 point
        assert combine_evidence(frozenset({EvidenceTag.PP3})) == (
            ACMGClassification.VUS
        )

    def test_single_moderate_tag_returns_vus(self) -> None:
        # 2 points
        assert combine_evidence(frozenset({EvidenceTag.PM2})) == (
            ACMGClassification.VUS
        )

    def test_two_supporting_returns_vus(self) -> None:
        # 2 points
        assert (
            combine_evidence(frozenset({EvidenceTag.PP3, EvidenceTag.PP5}))
            == ACMGClassification.VUS
        )

    def test_moderate_plus_supporting_returns_vus(self) -> None:
        # 3 points
        assert (
            combine_evidence(frozenset({EvidenceTag.PM2, EvidenceTag.PP3}))
            == ACMGClassification.VUS
        )

    def test_two_moderate_returns_vus(self) -> None:
        # 4 points: the central correction, never Likely Pathogenic.
        assert (
            combine_evidence(frozenset({EvidenceTag.PM2, EvidenceTag.PM4}))
            == ACMGClassification.VUS
        )


class TestCombineEvidenceLikelyPathogenic:
    def test_very_strong_alone_is_likely_pathogenic(self) -> None:
        # 8 points falls in the 6..9 Likely Pathogenic band.
        assert combine_evidence(frozenset({EvidenceTag.PVS1})) == (
            ACMGClassification.LIKELY_PATHOGENIC
        )

    def test_very_strong_plus_supporting_is_likely_pathogenic(self) -> None:
        # 9 points
        assert (
            combine_evidence(frozenset({EvidenceTag.PVS1, EvidenceTag.PP3}))
            == ACMGClassification.LIKELY_PATHOGENIC
        )

    def test_strong_plus_moderate_is_likely_pathogenic(self) -> None:
        # 4 + 2 = 6 points
        assert (
            combine_evidence(frozenset({EvidenceTag.PS1, EvidenceTag.PM2}))
            == ACMGClassification.LIKELY_PATHOGENIC
        )


class TestCombineEvidencePathogenic:
    def test_very_strong_plus_moderate_is_pathogenic(self) -> None:
        # 8 + 2 = 10 points
        assert (
            combine_evidence(frozenset({EvidenceTag.PVS1, EvidenceTag.PM2}))
            == ACMGClassification.PATHOGENIC
        )

    def test_very_strong_plus_two_supporting_is_pathogenic(self) -> None:
        # 8 + 1 + 1 = 10 points
        assert (
            combine_evidence(
                frozenset({EvidenceTag.PVS1, EvidenceTag.PP3, EvidenceTag.PP5})
            )
            == ACMGClassification.PATHOGENIC
        )

    def test_all_four_tags_returns_pathogenic(self) -> None:
        # 8 + 2 + 1 + 1 = 12 points
        assert (
            combine_evidence(
                frozenset(
                    {
                        EvidenceTag.PVS1,
                        EvidenceTag.PM2,
                        EvidenceTag.PP3,
                        EvidenceTag.PP5,
                    }
                )
            )
            == ACMGClassification.PATHOGENIC
        )

    def test_two_very_strong_is_pathogenic(self) -> None:
        # 4 + 4 + 8 = 16 (PS1 + PS1-equivalent not expressible; use PVS1+PS1+PS1
        # surrogate). Two strong plus very strong all exceed threshold.
        assert (
            combine_evidence(frozenset({EvidenceTag.PVS1, EvidenceTag.PS1}))
            == ACMGClassification.PATHOGENIC
        )


class TestCombineEvidenceBenign:
    def test_ba1_is_standalone_benign(self) -> None:
        assert combine_evidence(frozenset({EvidenceTag.BA1})) == (
            ACMGClassification.BENIGN
        )

    def test_two_strong_benign_is_benign(self) -> None:
        # -4 + -4 = -8
        assert (
            combine_evidence(frozenset({EvidenceTag.BS1, EvidenceTag.BS2}))
            == ACMGClassification.BENIGN
        )

    def test_strong_benign_alone_is_likely_benign(self) -> None:
        # -4
        assert combine_evidence(frozenset({EvidenceTag.BS1})) == (
            ACMGClassification.LIKELY_BENIGN
        )


class TestConflictResolvedBySummation:
    def test_opposing_evidence_nets_to_vus(self) -> None:
        # PVS1 (8) + BS1 (-4) = 4 -> VUS, by arithmetic not by veto.
        tags = frozenset({EvidenceTag.PVS1, EvidenceTag.BS1})
        assert combine_evidence(tags) == ACMGClassification.VUS
        assert has_conflicting_evidence(tags) is True

    def test_dominant_pathogenic_survives_a_supporting_benign(self) -> None:
        # PVS1 (8) + PM2 (2) + BP4 (-1) = 9 -> Likely Pathogenic, not vetoed.
        tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PM2, EvidenceTag.BP4})
        assert combine_evidence(tags) == ACMGClassification.LIKELY_PATHOGENIC
        assert has_conflicting_evidence(tags) is True

    def test_conflict_flag_false_without_benign(self) -> None:
        assert has_conflicting_evidence(frozenset({EvidenceTag.PVS1})) is False
