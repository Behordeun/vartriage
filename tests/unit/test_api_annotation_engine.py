"""Unit tests for APIAnnotationEngine and APIScoreProvider."""

from __future__ import annotations

from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock, patch

import pytest

httpx = pytest.importorskip("httpx")

from vartriage.api._cache import ResponseCache
from vartriage.api.annotation_engine import APIAnnotationEngine
from vartriage.api.config import APIConfig
from vartriage.api.score_provider import APIScoreProvider
from vartriage.api.vep_client import VEPAnnotation
from vartriage.models.variant import (ClinVarAssertion, FunctionalConsequence,
                                      Variant)

# --- Helpers ---


def _make_variant(
    chrom: str = "chr22", pos: int = 17818804, ref: str = "G", alt: str = "A"
) -> Variant:
    return Variant(
        chrom=chrom, pos=pos, id=None, ref=ref, alt=alt, qual=30.0, filter_status="PASS"
    )


def _make_vep_annotation(
    consequence: FunctionalConsequence = FunctionalConsequence.MISSENSE,
    gene_name: str = "BRCA1",
    allele_frequency: Optional[float] = 0.002,
    cadd_phred: Optional[float] = 24.5,
) -> VEPAnnotation:
    return VEPAnnotation(
        consequence=consequence,
        gene_name=gene_name,
        allele_frequency=allele_frequency,
        cadd_phred=cadd_phred,
        transcript_id="ENST00000001",
        hgvsc="c.123A>G",
        hgvsp="p.Lys41Arg",
        ensembl_release="112",
    )


@pytest.fixture
def api_config(tmp_path: Path) -> APIConfig:
    return APIConfig(
        mode="api",
        genome_build="grch38",
        cache_path=tmp_path / "test_api.db",
        cache_ttl_days=1,
        vep_batch_size=5,
        max_retries=1,
        connect_timeout=5.0,
        read_timeout=10.0,
    )


# --- APIAnnotationEngine Tests ---


class TestAPIAnnotationEngineAnnotation:
    """Core annotation flow: VEP + ClinVar composed into AnnotatedVariant."""

    def test_annotates_single_variant_with_vep_and_clinvar(
        self, api_config: APIConfig
    ) -> None:
        engine = APIAnnotationEngine(api_config)

        # Mock the VEP and ClinVar clients
        mock_vep_result = _make_vep_annotation(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="MICAL3",
            allele_frequency=0.00215,
            cadd_phred=24.5,
        )
        engine._vep.annotate_batch = MagicMock(return_value=[mock_vep_result])
        engine._clinvar.lookup_batch = MagicMock(return_value=[ClinVarAssertion.BENIGN])

        variants = [_make_variant()]
        results = list(engine.annotate(iter(variants)))

        assert len(results) == 1
        ann = results[0]
        assert ann.consequence == FunctionalConsequence.MISSENSE
        assert ann.gene_name == "MICAL3"
        assert ann.allele_frequency == pytest.approx(0.00215)
        assert ann.clinvar_assertion == ClinVarAssertion.BENIGN
        assert ann.frequency_unknown is False
        assert ann.clinvar_unknown is False

    def test_annotates_multiple_variants_in_batch(self, api_config: APIConfig) -> None:
        engine = APIAnnotationEngine(api_config)

        vep_results = [
            _make_vep_annotation(gene_name="GENE1"),
            _make_vep_annotation(
                consequence=FunctionalConsequence.FRAMESHIFT,
                gene_name="GENE2",
                allele_frequency=None,
            ),
        ]
        clinvar_results = [ClinVarAssertion.PATHOGENIC, None]

        engine._vep.annotate_batch = MagicMock(return_value=vep_results)
        engine._clinvar.lookup_batch = MagicMock(return_value=clinvar_results)

        variants = [_make_variant(pos=100), _make_variant(pos=200)]
        results = list(engine.annotate(iter(variants)))

        assert len(results) == 2
        assert results[0].gene_name == "GENE1"
        assert results[0].clinvar_assertion == ClinVarAssertion.PATHOGENIC
        assert results[1].consequence == FunctionalConsequence.FRAMESHIFT
        assert results[1].clinvar_assertion is None
        assert results[1].clinvar_unknown is True

    def test_vep_failure_produces_intergenic_with_warning(
        self, api_config: APIConfig
    ) -> None:
        engine = APIAnnotationEngine(api_config)

        # VEP returns None (failed)
        engine._vep.annotate_batch = MagicMock(return_value=[None])
        engine._clinvar.lookup_batch = MagicMock(return_value=[None])

        variants = [_make_variant()]
        results = list(engine.annotate(iter(variants)))

        assert results[0].consequence == FunctionalConsequence.INTERGENIC
        assert results[0].gene_name is None

        # Should have a VEP warning
        vep_warnings = [w for w in engine.warnings if w.source == "VEP"]
        assert len(vep_warnings) == 1
        assert vep_warnings[0].reason == "api_error"

    def test_missing_frequency_tracked_in_warnings(self, api_config: APIConfig) -> None:
        engine = APIAnnotationEngine(api_config)

        # VEP returns annotation but no frequency
        vep_ann = _make_vep_annotation(allele_frequency=None)
        engine._vep.annotate_batch = MagicMock(return_value=[vep_ann])
        engine._clinvar.lookup_batch = MagicMock(return_value=[None])

        variants = [_make_variant()]
        results = list(engine.annotate(iter(variants)))

        assert results[0].frequency_unknown is True
        gnomad_warnings = [w for w in engine.warnings if w.source == "gnomAD"]
        assert len(gnomad_warnings) == 1

    def test_clinvar_miss_tracked_in_warnings(self, api_config: APIConfig) -> None:
        engine = APIAnnotationEngine(api_config)

        engine._vep.annotate_batch = MagicMock(return_value=[_make_vep_annotation()])
        engine._clinvar.lookup_batch = MagicMock(return_value=[None])

        variants = [_make_variant()]
        list(engine.annotate(iter(variants)))

        clinvar_warnings = [w for w in engine.warnings if w.source == "ClinVar"]
        assert len(clinvar_warnings) == 1
        assert clinvar_warnings[0].reason == "not_found"

    def test_batches_variants_by_config_batch_size(self, api_config: APIConfig) -> None:
        """With batch_size=5, 7 variants should produce 2 VEP calls."""
        engine = APIAnnotationEngine(api_config)

        call_args: list[int] = []

        def mock_vep_batch(keys: list) -> list:
            call_args.append(len(keys))
            return [_make_vep_annotation() for _ in keys]

        engine._vep.annotate_batch = MagicMock(side_effect=mock_vep_batch)
        engine._clinvar.lookup_batch = MagicMock(
            side_effect=lambda keys: [None] * len(keys)
        )

        variants = [_make_variant(pos=i) for i in range(7)]
        list(engine.annotate(iter(variants)))

        assert call_args == [5, 2]  # batch_size=5 splits into 5+2


# --- APIScoreProvider Tests ---


class TestAPIScoreProviderCADD:
    """CADD score hierarchy: VEP plugin first, standalone API fallback."""

    def test_uses_vep_cadd_when_available(self, api_config: APIConfig) -> None:
        provider = APIScoreProvider(api_config)

        vep_annotations = [
            _make_vep_annotation(cadd_phred=32.0),
            _make_vep_annotation(cadd_phred=18.5),
        ]
        keys = [("chr1", 100, "A", "T"), ("chr1", 200, "G", "C")]

        results = provider.lookup_cadd_batch(keys, vep_annotations=vep_annotations)

        assert results[0] == pytest.approx(32.0)
        assert results[1] == pytest.approx(18.5)

    def test_falls_back_to_api_when_vep_cadd_missing(
        self, api_config: APIConfig
    ) -> None:
        provider = APIScoreProvider(api_config)

        # VEP has CADD for first, None for second
        vep_annotations = [
            _make_vep_annotation(cadd_phred=25.0),
            _make_vep_annotation(cadd_phred=None),
        ]

        # Mock the CADD API to return a score for the second variant
        provider._cadd.lookup_batch = MagicMock(return_value=[20.0])

        keys = [("chr1", 100, "A", "T"), ("chr1", 200, "G", "C")]
        results = provider.lookup_cadd_batch(keys, vep_annotations=vep_annotations)

        assert results[0] == pytest.approx(25.0)  # from VEP
        assert results[1] == pytest.approx(20.0)  # from API fallback

        # Only one API call should have been made (for the missing one)
        provider._cadd.lookup_batch.assert_called_once()
        call_keys = provider._cadd.lookup_batch.call_args[0][0]
        assert len(call_keys) == 1
        assert call_keys[0] == ("chr1", 200, "G", "C")

    def test_no_vep_annotations_queries_all_via_api(
        self, api_config: APIConfig
    ) -> None:
        provider = APIScoreProvider(api_config)
        provider._cadd.lookup_batch = MagicMock(return_value=[15.0, 22.0])

        keys = [("chr1", 100, "A", "T"), ("chr1", 200, "G", "C")]
        results = provider.lookup_cadd_batch(keys, vep_annotations=None)

        assert results[0] == pytest.approx(15.0)
        assert results[1] == pytest.approx(22.0)

    def test_api_fallback_returns_none_on_failure(self, api_config: APIConfig) -> None:
        provider = APIScoreProvider(api_config)
        vep_annotations = [_make_vep_annotation(cadd_phred=None)]
        provider._cadd.lookup_batch = MagicMock(return_value=[None])

        keys = [("chr1", 100, "A", "T")]
        results = provider.lookup_cadd_batch(keys, vep_annotations=vep_annotations)

        assert results[0] is None


class TestAPIScoreProviderSpliceAI:
    """SpliceAI lookups with consequence filtering."""

    def test_delegates_to_spliceai_client(self, api_config: APIConfig) -> None:
        provider = APIScoreProvider(api_config)
        provider._spliceai.lookup_batch = MagicMock(return_value=[0.85, None])

        keys = [("chr1", 100, "A", "T"), ("chr2", 200, "G", "C")]
        consequences = [
            FunctionalConsequence.SPLICE_SITE,
            FunctionalConsequence.INTERGENIC,
        ]

        results = provider.lookup_spliceai_batch(keys, consequences=consequences)

        assert results[0] == pytest.approx(0.85)
        assert results[1] is None


class TestAPIScoreProviderREVEL:
    """REVEL limitation in API mode."""

    def test_revel_returns_all_none(self, api_config: APIConfig) -> None:
        provider = APIScoreProvider(api_config)
        keys = [("chr1", 100, "A", "T"), ("chr2", 200, "G", "C")]

        results = provider.lookup_revel_batch(keys)

        assert results == [None, None]

    def test_revel_logs_warning_once(self, api_config: APIConfig) -> None:
        provider = APIScoreProvider(api_config)
        keys = [("chr1", 100, "A", "T")]

        with patch("vartriage.api.score_provider.logger") as mock_logger:
            provider.lookup_revel_batch(keys)
            provider.lookup_revel_batch(keys)

            # Warning logged only on first call
            assert mock_logger.warning.call_count == 1
            assert "REVEL" in mock_logger.warning.call_args[0][0]
