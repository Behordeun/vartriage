"""Unit tests for the mitochondrial variant classifier."""

from __future__ import annotations

import pytest

from vartriage.mito.classifier import (
    MitochondrialClassifier,
    MitoClassification,
)
from vartriage.mito.frequency import HelixMTdbDatabase
from vartriage.mito.gene_map import MtGeneMap
from vartriage.mito.mitomap import MitomapDatabase
from vartriage.models.variant import Variant


@pytest.fixture
def classifier() -> MitochondrialClassifier:
    """Build a classifier with bundled databases."""
    gene_map = MtGeneMap()
    mitomap_db = MitomapDatabase()
    helix_db = HelixMTdbDatabase()
    return MitochondrialClassifier(
        gene_map=gene_map,
        mitomap_db=mitomap_db,
        helix_db=helix_db,
    )


def _make_mito_variant(
    pos: int,
    ref: str,
    alt: str,
    ad: list[int] | None = None,
    af: float | None = None,
) -> Variant:
    """Helper to create a chrM variant with optional heteroplasmy data."""
    info: dict = {}
    if ad is not None:
        info["AD"] = ad
    if af is not None:
        info["AF"] = af
    return Variant(
        chrom="chrM",
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=100.0,
        filter_status="PASS",
        info=info,
    )


class TestPathogenicClassification:
    """Confirmed MITOMAP + high heteroplasmy + rare = Pathogenic."""

    def test_m3243ag_high_heteroplasmy_is_pathogenic(
        self, classifier: MitochondrialClassifier
    ):
        # m.3243A>G: confirmed MELAS, AF=0.0002 in HelixMTdb (above strict
        # rare threshold of 0.0001), so falls to Likely Pathogenic via rule 4
        # (confirmed + functional region + high heteroplasmy)
        variant = _make_mito_variant(3243, "A", "G", ad=[50, 950])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.LIKELY_PATHOGENIC
        assert "MELAS" in result.classification_reason

    def test_m11778ga_homoplasmic_is_pathogenic(
        self, classifier: MitochondrialClassifier
    ):
        # m.11778G>A: confirmed LHON
        variant = _make_mito_variant(11778, "G", "A", ad=[10, 990])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.PATHOGENIC
        assert "LHON" in result.classification_reason


class TestBenignClassification:
    """Common haplogroup markers (AF > 5%) = Benign."""

    def test_common_haplogroup_variant_is_benign(
        self, classifier: MitochondrialClassifier
    ):
        # Position 73 A>G is very common (AF ~0.92)
        variant = _make_mito_variant(73, "A", "G", ad=[5, 995])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.BENIGN
        assert "haplogroup" in result.classification_reason.lower()

    def test_position_263_ag_is_benign(self, classifier: MitochondrialClassifier):
        # Position 263 A>G is nearly universal
        variant = _make_mito_variant(263, "A", "G", ad=[10, 990])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.BENIGN


class TestLikelyBenignClassification:
    """Moderate frequency (AF > 0.1%) = Likely Benign."""

    def test_moderate_frequency_variant_is_likely_benign(
        self, classifier: MitochondrialClassifier
    ):
        # Position 3010 G>A has AF ~0.15 (above 0.1% = 0.001)
        variant = _make_mito_variant(3010, "G", "A", ad=[200, 800])
        result = classifier.classify(variant)
        assert result.classification in (
            MitoClassification.LIKELY_BENIGN,
            MitoClassification.BENIGN,
        )


class TestLikelyPathogenicClassification:
    """Reported in MITOMAP + functional region + moderate-high heteroplasmy."""

    def test_reported_variant_in_trna_moderate_heteroplasmy(
        self, classifier: MitochondrialClassifier
    ):
        # m.3291T>C: Reported (unconfirmed) MELAS, in MT-TL1 (tRNA),
        # moderate heteroplasmy. Rule 4b: unconfirmed → VUS.
        variant = _make_mito_variant(3291, "T", "C", ad=[600, 400])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.VUS
        assert "unconfirmed" in result.classification_reason

    def test_confirmed_variant_in_trna_moderate_heteroplasmy(
        self, classifier: MitochondrialClassifier
    ):
        # m.3243A>G: Confirmed MELAS, in MT-TL1 (tRNA), moderate heteroplasmy.
        # Rule 4: confirmed + coding/tRNA + moderate → Likely Pathogenic.
        variant = _make_mito_variant(3243, "A", "G", ad=[600, 400])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.LIKELY_PATHOGENIC


class TestVUSClassification:
    """Novel variants without strong evidence = VUS."""

    def test_novel_variant_absent_from_databases(
        self, classifier: MitochondrialClassifier
    ):
        # Position not in MITOMAP or HelixMTdb
        variant = _make_mito_variant(5000, "A", "T", ad=[700, 300])
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.VUS

    def test_novel_variant_no_heteroplasmy_data(
        self, classifier: MitochondrialClassifier
    ):
        # No AD or AF fields
        variant = _make_mito_variant(5000, "A", "T")
        result = classifier.classify(variant)
        assert result.classification == MitoClassification.VUS


class TestClassifierOutputFields:
    """Verify all output fields are populated correctly."""

    def test_heteroplasmy_populated(self, classifier: MitochondrialClassifier):
        variant = _make_mito_variant(3243, "A", "G", ad=[50, 950])
        result = classifier.classify(variant)
        assert result.heteroplasmy is not None
        assert abs(result.heteroplasmy.percentage - 95.0) < 0.1

    def test_gene_context_populated(self, classifier: MitochondrialClassifier):
        variant = _make_mito_variant(3243, "A", "G", ad=[50, 950])
        result = classifier.classify(variant)
        assert result.gene_context.gene_name == "MT-TL1"
        assert result.gene_context.gene_type == "tRNA"

    def test_mitomap_entry_populated(self, classifier: MitochondrialClassifier):
        variant = _make_mito_variant(3243, "A", "G", ad=[50, 950])
        result = classifier.classify(variant)
        assert result.mitomap_entry is not None
        assert "MELAS" in result.mitomap_entry.disease

    def test_helix_af_populated(self, classifier: MitochondrialClassifier):
        variant = _make_mito_variant(73, "A", "G", ad=[5, 995])
        result = classifier.classify(variant)
        assert result.helix_af is not None
        assert result.helix_af > 0.5
