"""Unit tests for the MITOMAP database module."""

from __future__ import annotations

import pytest

from vartriage.mito.mitomap import MitomapDatabase


@pytest.fixture
def mitomap_db() -> MitomapDatabase:
    """Load the bundled MITOMAP database."""
    return MitomapDatabase()


class TestMitomapDatabaseLoading:
    """Verify the MITOMAP database loads correctly."""

    def test_loads_entries(self, mitomap_db: MitomapDatabase):
        assert mitomap_db.size > 0

    def test_loads_at_least_50_entries(self, mitomap_db: MitomapDatabase):
        assert mitomap_db.size >= 50


class TestMitomapLookup:
    """Verify MITOMAP pathogenic variant lookups."""

    def test_m3243ag_returns_melas(self, mitomap_db: MitomapDatabase):
        entry = mitomap_db.lookup(3243, "A", "G")
        assert entry is not None
        assert "MELAS" in entry.disease
        assert entry.status == "Cfrm"
        assert entry.locus == "MT-TL1"

    def test_m3243ag_is_confirmed(self, mitomap_db: MitomapDatabase):
        entry = mitomap_db.lookup(3243, "A", "G")
        assert entry is not None
        assert entry.is_confirmed is True

    def test_m11778ga_returns_lhon(self, mitomap_db: MitomapDatabase):
        entry = mitomap_db.lookup(11778, "G", "A")
        assert entry is not None
        assert "LHON" in entry.disease
        assert entry.status == "Cfrm"
        assert entry.locus == "MT-ND4"

    def test_m8993tg_returns_narp(self, mitomap_db: MitomapDatabase):
        entry = mitomap_db.lookup(8993, "T", "G")
        assert entry is not None
        assert "NARP" in entry.disease
        assert entry.status == "Cfrm"

    def test_unknown_position_returns_none(self, mitomap_db: MitomapDatabase):
        result = mitomap_db.lookup(9999, "A", "T")
        assert result is None

    def test_case_insensitive_lookup(self, mitomap_db: MitomapDatabase):
        entry = mitomap_db.lookup(3243, "a", "g")
        assert entry is not None
        assert "MELAS" in entry.disease

    def test_reported_status_not_confirmed(self, mitomap_db: MitomapDatabase):
        # m.3291T>C is Reported, not Cfrm
        entry = mitomap_db.lookup(3291, "T", "C")
        assert entry is not None
        assert entry.status == "Reported"
        assert entry.is_confirmed is False
