"""ClinGen SVI Bayesian point engine: point assignment and tier thresholds.

Point values and thresholds follow Tavtigian et al. 2018/2020 and the
ClinGen SVI point-based framework: Supporting 1, Moderate 2, Strong 4,
Very Strong 8; benign criteria symmetric negative; BA1 standalone Benign.
Tier thresholds on summed points: Pathogenic >= 10, Likely Pathogenic
6..9, VUS -5..5, Likely Benign -6..-1, Benign <= -7.
"""

from __future__ import annotations

from vartriage.classification.points import (
    classify_by_points,
    score_points,
)
from vartriage.models.variant import ACMGClassification, EvidenceTag


class TestPointAssignment:
    def test_supporting_pathogenic_is_one_point(self) -> None:
        assert score_points(frozenset({EvidenceTag.PP3})) == 1

    def test_moderate_pathogenic_is_two_points(self) -> None:
        assert score_points(frozenset({EvidenceTag.PM2})) == 2

    def test_strong_pathogenic_is_four_points(self) -> None:
        assert score_points(frozenset({EvidenceTag.PS1})) == 4

    def test_very_strong_pathogenic_is_eight_points(self) -> None:
        assert score_points(frozenset({EvidenceTag.PVS1})) == 8

    def test_strong_benign_is_minus_four_points(self) -> None:
        assert score_points(frozenset({EvidenceTag.BS1})) == -4

    def test_supporting_benign_is_minus_one_point(self) -> None:
        assert score_points(frozenset({EvidenceTag.BP4})) == -1

    def test_moderate_benign_is_minus_two_points(self) -> None:
        assert score_points(frozenset({EvidenceTag.BP4_MODERATE})) == -2

    def test_points_sum_across_tags(self) -> None:
        # 1 Strong (4) + 1 Moderate (2) + 1 Supporting (1) = 7
        tags = frozenset({EvidenceTag.PS1, EvidenceTag.PM2, EvidenceTag.PP3})
        assert score_points(tags) == 7

    def test_opposing_evidence_nets_out(self) -> None:
        # PVS1 (8) + BS1 (-4) = 4
        tags = frozenset({EvidenceTag.PVS1, EvidenceTag.BS1})
        assert score_points(tags) == 4


class TestTierThresholds:
    def test_empty_is_vus(self) -> None:
        assert classify_by_points(frozenset()) == ACMGClassification.VUS

    def test_pathogenic_at_ten_points(self) -> None:
        # PVS1 (8) + PM2 (2) = 10 -> Pathogenic
        tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PM2})
        assert classify_by_points(tags) == ACMGClassification.PATHOGENIC

    def test_likely_pathogenic_at_six_points(self) -> None:
        # PS1 (4) + PM2 (2) = 6 -> Likely Pathogenic
        tags = frozenset({EvidenceTag.PS1, EvidenceTag.PM2})
        assert classify_by_points(tags) == ACMGClassification.LIKELY_PATHOGENIC

    def test_likely_pathogenic_upper_boundary_nine(self) -> None:
        # PVS1 (8) + PP3 (1) = 9 -> Likely Pathogenic (not Pathogenic)
        tags = frozenset({EvidenceTag.PVS1, EvidenceTag.PP3})
        assert classify_by_points(tags) == ACMGClassification.LIKELY_PATHOGENIC

    def test_two_moderate_is_vus_not_likely_pathogenic(self) -> None:
        # The central correction: 2 Moderate = 4 points = VUS, never LP.
        tags = frozenset({EvidenceTag.PM2, EvidenceTag.PM4})
        assert classify_by_points(tags) == ACMGClassification.VUS

    def test_three_moderate_reaches_likely_pathogenic(self) -> None:
        # 3 Moderate = 6 points -> Likely Pathogenic (the real point route)
        tags = frozenset({EvidenceTag.PM1, EvidenceTag.PM2, EvidenceTag.PM4})
        assert classify_by_points(tags) == ACMGClassification.LIKELY_PATHOGENIC

    def test_single_moderate_is_vus(self) -> None:
        assert classify_by_points(frozenset({EvidenceTag.PM2})) == (
            ACMGClassification.VUS
        )

    def test_benign_at_minus_seven(self) -> None:
        # BS1 (-4) + BS2 (-4) = -8 -> Benign
        tags = frozenset({EvidenceTag.BS1, EvidenceTag.BS2})
        assert classify_by_points(tags) == ACMGClassification.BENIGN

    def test_likely_benign_at_minus_six(self) -> None:
        # BS1 (-4) + BP4_MODERATE (-2) = -6 -> Likely Benign
        tags = frozenset({EvidenceTag.BS1, EvidenceTag.BP4_MODERATE})
        assert classify_by_points(tags) == ACMGClassification.LIKELY_BENIGN

    def test_likely_benign_at_minus_one(self) -> None:
        assert classify_by_points(frozenset({EvidenceTag.BP4})) == (
            ACMGClassification.LIKELY_BENIGN
        )

    def test_single_supporting_pathogenic_is_vus(self) -> None:
        assert classify_by_points(frozenset({EvidenceTag.PP3})) == (
            ACMGClassification.VUS
        )


class TestStandaloneAndConflict:
    def test_ba1_is_benign_regardless_of_points(self) -> None:
        # BA1 standalone override, even alongside strong pathogenic points.
        tags = frozenset({EvidenceTag.BA1, EvidenceTag.PVS1, EvidenceTag.PM2})
        assert classify_by_points(tags) == ACMGClassification.BENIGN

    def test_conflicting_evidence_resolves_by_summation(self) -> None:
        # PVS1 (8) + BS1 (-4) = 4 -> VUS, decided by the arithmetic, not a
        # boolean "both signs present" veto.
        tags = frozenset({EvidenceTag.PVS1, EvidenceTag.BS1})
        assert classify_by_points(tags) == ACMGClassification.VUS
