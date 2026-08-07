"""Unit tests for the HelixMTdb frequency database module."""

from __future__ import annotations

import pytest

from vartriage.mito.frequency import HelixMTdbDatabase


@pytest.fixture
def helix_db() -> HelixMTdbDatabase:
    """Load the bundled HelixMTdb database."""
    return HelixMTdbDatabase()


class TestHelixMTdbLoading:
    """Verify the HelixMTdb database loads correctly."""

    def test_loads_entries(self, helix_db: HelixMTdbDatabase):
        assert helix_db.size > 0

    def test_loads_at_least_100_entries(self, helix_db: HelixMTdbDatabase):
        assert helix_db.size >= 100


class TestHelixMTdbLookup:
    """Verify population frequency lookups."""

    def test_common_variant_has_high_af(self, helix_db: HelixMTdbDatabase):
        # Position 73 A>G is a very common variant (AF ~0.92)
        entry = helix_db.lookup(73, "A", "G")
        assert entry is not None
        assert entry.af > 0.01
        assert entry.is_common_haplogroup_marker is True

    def test_pathogenic_variant_is_rare(self, helix_db: HelixMTdbDatabase):
        # m.3243A>G (MELAS) has low AF but above the strict is_rare threshold (0.0001)
        entry = helix_db.lookup(3243, "A", "G")
        assert entry is not None
        assert entry.af < 0.001
        # AF=0.0002 is above the is_rare cutoff of 0.0001
        assert entry.is_rare is False

    def test_novel_position_returns_none(self, helix_db: HelixMTdbDatabase):
        result = helix_db.lookup(16570, "A", "T")
        assert result is None

    def test_get_af_returns_float(self, helix_db: HelixMTdbDatabase):
        af = helix_db.get_af(73, "A", "G")
        assert af is not None
        assert isinstance(af, float)
        assert af > 0.5

    def test_get_af_returns_none_for_unknown(self, helix_db: HelixMTdbDatabase):
        af = helix_db.get_af(16570, "A", "T")
        assert af is None

    def test_case_insensitive_lookup(self, helix_db: HelixMTdbDatabase):
        entry = helix_db.lookup(73, "a", "g")
        assert entry is not None

    def test_haplogroup_marker_threshold(self, helix_db: HelixMTdbDatabase):
        # Position 263 A>G is nearly universal (AF ~0.98)
        entry = helix_db.lookup(263, "A", "G")
        assert entry is not None
        assert entry.is_common_haplogroup_marker is True
        assert entry.is_rare is False
