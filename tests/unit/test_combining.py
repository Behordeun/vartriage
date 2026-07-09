"""Unit tests for ACMG/AMP 2015 evidence combining rules."""

from __future__ import annotations

import pytest

from vartriage.classification.combining import combine_evidence
from vartriage.models.variant import (
    ACMGClassification,
    EvidenceTag,
)


class TestCombineEvidenceVUS:
    """Cases that should produce VUS classification."""

    def test_empty_tag_set_returns_vus(self) -> None:
        result = combine_evidence(frozenset())
        assert result == ACMGClassification.VUS

    def test_single_supporting_tag_returns_vus(self) -> None:
        result = combine_evidence(frozenset({EvidenceTag.PP3}))
        assert result == ACMGClassification.VUS

    def test_single_moderate_tag_returns_vus(self) -> None:
        result = combine_evidence(frozenset({EvidenceTag.PM2}))
        assert result == ACMGClassification.VUS

    def test_two_supporting_without_strong_or_very_strong_returns_vus(self) -> None:
        result = combine_evidence(frozenset({EvidenceTag.PP3, EvidenceTag.PP5}))
        assert result == ACMGClassification.VUS

    def test_moderate_plus_supporting_returns_vus(self) -> None:
        result = combine_evidence(frozenset({EvidenceTag.PM2, EvidenceTag.PP3}))
        assert result == ACMGClassification.VUS


class TestCombineEvidenceLikelyPathogenic:
    """Cases that should produce Likely_Pathogenic classification."""

    def test_very_strong_plus_moderate(self) -> None:
        """PVS1 (Very Strong) + PM2 (Moderate) -> Likely Pathogenic."""
        result = combine_evidence(frozenset({EvidenceTag.PVS1, EvidenceTag.PM2}))
        assert result == ACMGClassification.LIKELY_PATHOGENIC


class TestCombineEvidencePathogenic:
    """Cases that should produce Pathogenic classification."""

    def test_very_strong_plus_two_supporting(self) -> None:
        """PVS1 (Very Strong) + PP3 + PP5 (2 Supporting) -> Pathogenic."""
        result = combine_evidence(
            frozenset({EvidenceTag.PVS1, EvidenceTag.PP3, EvidenceTag.PP5})
        )
        assert result == ACMGClassification.PATHOGENIC

    def test_all_four_tags_returns_pathogenic(self) -> None:
        """PVS1 + PM2 + PP3 + PP5 -> Pathogenic (Very Strong + 2 Supporting)."""
        result = combine_evidence(
            frozenset(
                {EvidenceTag.PVS1, EvidenceTag.PM2, EvidenceTag.PP3, EvidenceTag.PP5}
            )
        )
        assert result == ACMGClassification.PATHOGENIC


class TestCombineEvidenceSpecificScenarios:
    """Domain-specific scenarios from the library's tag set."""

    def test_pvs1_alone_is_vus(self) -> None:
        """A single Very Strong tag without additional evidence is VUS."""
        result = combine_evidence(frozenset({EvidenceTag.PVS1}))
        assert result == ACMGClassification.VUS

    def test_pvs1_plus_single_supporting_is_vus(self) -> None:
        """Very Strong + 1 Supporting doesn't meet any rule threshold."""
        result = combine_evidence(frozenset({EvidenceTag.PVS1, EvidenceTag.PP3}))
        assert result == ACMGClassification.VUS

    def test_pvs1_plus_moderate_plus_supporting_is_pathogenic(self) -> None:
        """PVS1 + PM2 + PP3: Very Strong meets LP via Moderate, but also
        Very Strong + Moderate + Supporting -> the LP rule fires, but we also
        check if Pathogenic applies. PVS1(VS) + PP3(Sup) = only 1 Supporting,
        not enough for Pathogenic. So result is Likely Pathogenic.

        Wait: VS=1, M=1, Sup=1. Pathogenic needs VS>=1 + S>=1 (no Strong),
        or S>=2 + Sup>=1 (no), or VS>=1 + Sup>=2 (only 1 Supporting).
        None satisfied -> check LP: VS>=1 + M>=1 -> yes -> Likely Pathogenic.
        """
        result = combine_evidence(
            frozenset({EvidenceTag.PVS1, EvidenceTag.PM2, EvidenceTag.PP3})
        )
        assert result == ACMGClassification.LIKELY_PATHOGENIC
