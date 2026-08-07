"""Integration tests for the mitochondrial pipeline.

Tests the full path from Variant objects through the MitochondrialPipeline,
verifying correct classification of known pathogenic, benign, and novel variants.
"""

from __future__ import annotations

import pytest

from vartriage.mito.classifier import MitoClassification
from vartriage.mito.config import MitoConfig
from vartriage.mito.pipeline import MitochondrialPipeline
from vartriage.models.variant import Variant


def _make_variant(
    pos: int,
    ref: str,
    alt: str,
    chrom: str = "chrM",
    ad: list[int] | None = None,
    af: float | None = None,
) -> Variant:
    """Build a Variant with optional heteroplasmy FORMAT fields."""
    info: dict = {}
    if ad is not None:
        info["AD"] = ad
    if af is not None:
        info["AF"] = af
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=200.0,
        filter_status="PASS",
        info=info,
    )


@pytest.fixture
def mito_pipeline() -> MitochondrialPipeline:
    """Build the mitochondrial pipeline with default config."""
    return MitochondrialPipeline(MitoConfig())


class TestKnownPathogenicVariant:
    """m.11778G>A (LHON) at high heteroplasmy should classify as Pathogenic."""

    def test_full_pipeline_classification(self, mito_pipeline: MitochondrialPipeline):
        # m.11778G>A: confirmed LHON, AF=0.00005 in HelixMTdb (rare)
        variant = _make_variant(11778, "G", "A", ad=[20, 980])
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 1
        result = results[0]
        assert result.classification == MitoClassification.PATHOGENIC
        assert result.heteroplasmy is not None
        assert result.heteroplasmy.category == "homoplasmic"
        assert result.mitomap_entry is not None
        assert "LHON" in result.mitomap_entry.disease
        assert result.gene_context.gene_name == "MT-ND4"

    def test_m8344ag_merrf_pathogenic(self, mito_pipeline: MitochondrialPipeline):
        # m.8344A>G: confirmed MERRF, AF=0.00004 (rare)
        variant = _make_variant(8344, "A", "G", ad=[30, 970])
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 1
        assert results[0].classification == MitoClassification.PATHOGENIC
        assert "MERRF" in results[0].mitomap_entry.disease


class TestCommonHaplogroupVariant:
    """Common haplogroup-defining polymorphisms should classify as Benign."""

    def test_position_73_ag_benign(self, mito_pipeline: MitochondrialPipeline):
        # Position 73 A>G is near-universal (AF ~0.92)
        variant = _make_variant(73, "A", "G", ad=[5, 995])
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 1
        assert results[0].classification == MitoClassification.BENIGN

    def test_position_750_ag_benign(self, mito_pipeline: MitochondrialPipeline):
        # Position 750 A>G: AF ~0.99
        variant = _make_variant(750, "A", "G", ad=[10, 990])
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 1
        assert results[0].classification == MitoClassification.BENIGN


class TestHeteroplasmyFiltering:
    """Variants below min_heteroplasmy threshold should be filtered out."""

    def test_sub_threshold_variant_filtered(self):
        config = MitoConfig(min_heteroplasmy=5.0)
        pipeline = MitochondrialPipeline(config)

        # 0.5% heteroplasmy, below the 5% threshold
        variant = _make_variant(3243, "A", "G", ad=[995, 5])
        results = pipeline.run(iter([variant]))

        assert len(results) == 0

    def test_above_threshold_variant_retained(self):
        config = MitoConfig(min_heteroplasmy=5.0)
        pipeline = MitochondrialPipeline(config)

        # 10% heteroplasmy, above the 5% threshold
        variant = _make_variant(3243, "A", "G", ad=[900, 100])
        results = pipeline.run(iter([variant]))

        assert len(results) == 1

    def test_default_threshold_filters_below_1_percent(
        self, mito_pipeline: MitochondrialPipeline
    ):
        # 0.3% heteroplasmy, below the default 1% threshold
        variant = _make_variant(5000, "A", "T", ad=[997, 3])
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 0


class TestMultipleVariants:
    """Verify batch processing of multiple chrM variants."""

    def test_mixed_classifications(self, mito_pipeline: MitochondrialPipeline):
        variants = [
            # Pathogenic: confirmed LHON, rare, high heteroplasmy
            _make_variant(11778, "G", "A", ad=[20, 980]),
            # Benign: common haplogroup marker
            _make_variant(73, "A", "G", ad=[5, 995]),
            # Novel VUS
            _make_variant(5000, "A", "T", ad=[500, 500]),
        ]
        results = mito_pipeline.run(iter(variants))

        assert len(results) == 3
        classifications = [r.classification for r in results]
        assert MitoClassification.PATHOGENIC in classifications
        assert MitoClassification.BENIGN in classifications
        assert MitoClassification.VUS in classifications

    def test_no_variants_returns_empty(self, mito_pipeline: MitochondrialPipeline):
        results = mito_pipeline.run(iter([]))
        assert results == []


class TestMitoConfigDisabled:
    """Verify that skip-mito config prevents processing."""

    def test_disabled_config_validation(self):
        config = MitoConfig(enabled=False)
        assert config.enabled is False

    def test_config_validation_rejects_invalid_threshold(self):
        with pytest.raises(ValueError, match="min_heteroplasmy"):
            MitoConfig(min_heteroplasmy=200.0)

    def test_config_validation_rejects_negative_threshold(self):
        with pytest.raises(ValueError, match="min_heteroplasmy"):
            MitoConfig(min_heteroplasmy=-5.0)


class TestNoHeteroplasmyData:
    """Variants without AD/AF fields should still be classified."""

    def test_confirmed_pathogenic_without_heteroplasmy(
        self, mito_pipeline: MitochondrialPipeline
    ):
        # m.11778G>A without heteroplasmy data: should be Likely Pathogenic
        # (confirmed + rare, but no heteroplasmy measurement)
        variant = _make_variant(11778, "G", "A")
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 1
        result = results[0]
        assert result.heteroplasmy is None
        assert result.classification == MitoClassification.LIKELY_PATHOGENIC

    def test_novel_without_heteroplasmy_is_vus(
        self, mito_pipeline: MitochondrialPipeline
    ):
        variant = _make_variant(5000, "A", "T")
        results = mito_pipeline.run(iter([variant]))

        assert len(results) == 1
        assert results[0].classification == MitoClassification.VUS
