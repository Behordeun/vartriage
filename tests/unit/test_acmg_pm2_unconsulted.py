"""Tests ensuring unconsulted frequency sources never trigger ACMG PM2.

In clinical genomics and ACMG/AMP variant interpretation, PM2 encodes
evidence of rarity or confirmed absence from population cohorts (e.g. gnomAD).
The absence of evidence (an unconsulted reference database, a network failure,
or an unannotated run) must not be conflated with evidence of absence.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from vartriage.annotation.engine import AnnotationEngine
from vartriage.annotation.frequency import DictFrequencyDatabase
from vartriage.annotation.frequency_polars import (
    POLARS_AVAILABLE,
    PolarsFrequencyDatabase,
)
from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase
from vartriage.api.annotation_engine import APIAnnotationEngine
from vartriage.api.vep_client import VEPAnnotation
from vartriage.classification.acmg import ACMGClassifier
from vartriage.models.config import AnnotationConfig
from vartriage.models.resolution import ResolutionState
from vartriage.models.variant import (
    AnnotatedVariant,
    EvidenceTag,
    FunctionalConsequence,
    PopulationFrequencies,
    ScoredVariant,
    Variant,
)
from vartriage.pipeline import Pipeline


def _make_variant(
    chrom: str = "chr1", pos: int = 1000, ref: str = "A", alt: str = "G"
) -> Variant:
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=None,
        filter_status="PASS",
        info={},
    )


def _make_scored(
    *,
    allele_frequency: float | None = None,
    frequency_unknown: bool = False,
    frequency_state: ResolutionState | None = None,
    population_frequencies: PopulationFrequencies | None = None,
) -> ScoredVariant:
    annotated = AnnotatedVariant(
        variant=_make_variant(),
        consequence=FunctionalConsequence.MISSENSE,
        gene_name="TESTG",
        allele_frequency=allele_frequency,
        frequency_unknown=frequency_unknown,
        frequency_state=frequency_state,
        population_frequencies=population_frequencies,
    )
    return ScoredVariant(annotated=annotated)


class TestPM2UnconsultedDefused:
    """Verify PM2 evaluates only on observed frequency or confirmed absence."""

    def test_unconsulted_frequency_never_fires_pm2(self) -> None:
        """Unannotated/unconsulted frequency must never fire PM2 and must log missing source."""
        scored = _make_scored(
            allele_frequency=None,
            frequency_unknown=False,
            frequency_state=None,
        )
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm2(scored, tags, missing)

        assert EvidenceTag.PM2 not in tags
        assert "gnomAD" in missing

    def test_lookup_failed_frequency_never_fires_pm2(self) -> None:
        """Transient lookup failure or circuit breaker open must not fire PM2."""
        scored = _make_scored(
            allele_frequency=None,
            frequency_unknown=False,
            frequency_state=ResolutionState.LOOKUP_FAILED,
        )
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm2(scored, tags, missing)

        assert EvidenceTag.PM2 not in tags
        assert "gnomAD" in missing

    def test_confirmed_absent_state_fires_pm2(self) -> None:
        """Confirmed absence from population database legitimately fires PM2."""
        scored = _make_scored(
            allele_frequency=None,
            frequency_unknown=True,
            frequency_state=ResolutionState.CONFIRMED_ABSENT,
        )
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm2(scored, tags, missing)

        assert EvidenceTag.PM2 in tags
        assert "gnomAD" not in missing

    def test_observed_rare_frequency_fires_pm2(self) -> None:
        """Observed AF below threshold fires PM2 without missing source warning."""
        scored = _make_scored(
            allele_frequency=0.000005,
            frequency_unknown=False,
            frequency_state=ResolutionState.FOUND,
        )
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm2(scored, tags, missing)

        assert EvidenceTag.PM2 in tags
        assert "gnomAD" not in missing

    def test_observed_common_frequency_does_not_fire_pm2(self) -> None:
        """Common allele frequency (> threshold) rejects PM2."""
        scored = _make_scored(
            allele_frequency=0.02,
            frequency_unknown=False,
            frequency_state=ResolutionState.FOUND,
        )
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()

        classifier._evaluate_pm2(scored, tags, missing)

        assert EvidenceTag.PM2 not in tags
        assert "gnomAD" not in missing


class TestAnnotationEngineUnconsultedBehavior:
    """Verify AnnotationEngine sets unconsulted state when gnomad_path is None."""

    def test_omitted_gnomad_yields_unconsulted_state(self, tmp_path: Path) -> None:
        gtf_file = tmp_path / "annotations.gtf"
        gtf_file.write_text("")
        config = AnnotationConfig(
            gene_annotation_path=gtf_file,
            gnomad_path=None,
            clinvar_path=None,
        )
        engine = AnnotationEngine(config)

        variant = _make_variant("chr1", 1000, "A", "G")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        res = results[0]
        assert res.allele_frequency is None
        assert res.frequency_unknown is False
        assert res.frequency_state is None

        # Confirmed that no spurious gnomAD not_found warning is emitted
        gnomad_warnings = [w for w in engine.warnings if w.source == "gnomAD"]
        assert len(gnomad_warnings) == 0


class TestPipelinePassthroughUnconsulted:
    """Verify Pipeline._passthrough_annotation produces unconsulted frequency status."""

    def test_passthrough_preserves_unconsulted_status(self) -> None:
        pipeline = Pipeline.__new__(Pipeline)
        variants = [_make_variant()]
        results = list(pipeline._passthrough_annotation(iter(variants)))

        assert len(results) == 1
        assert results[0].frequency_unknown is False
        assert results[0].frequency_state is None


class TestLegacyBackendFallback:
    """Verify AnnotationEngine maps legacy lookup_batch None to LOOKUP_FAILED."""

    def test_legacy_backend_none_maps_to_lookup_failed(self, tmp_path: Path) -> None:
        gtf_file = tmp_path / "annotations.gtf"
        gtf_file.write_text("")
        config = AnnotationConfig(
            gene_annotation_path=gtf_file,
            gnomad_path=None,
            clinvar_path=None,
        )
        engine = AnnotationEngine(config)

        # Mock a legacy backend that only implements lookup_batch
        legacy_mock = MagicMock(spec=["lookup_batch"])
        legacy_mock.lookup_batch.return_value = [None]
        engine._frequency_db = legacy_mock

        variant = _make_variant("chr1", 1000, "A", "G")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        res = results[0]
        assert res.allele_frequency is None
        assert res.frequency_unknown is False
        assert res.frequency_state is ResolutionState.LOOKUP_FAILED

        # Classification must NOT fire PM2 and must add gnomAD to missing sources
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()
        classifier._evaluate_pm2(ScoredVariant(annotated=res), tags, missing)

        assert EvidenceTag.PM2 not in tags
        assert "gnomAD" in missing


class TestVEPNoFrequencyFallback:
    """Verify APIAnnotationEngine maps VEP missing frequency to LOOKUP_FAILED."""

    def test_vep_missing_frequency_maps_to_lookup_failed(self) -> None:
        engine = APIAnnotationEngine.__new__(APIAnnotationEngine)
        engine._warnings = []
        engine._batch_size = 100

        mock_vep_client = MagicMock()
        mock_vep_ann = VEPAnnotation(
            consequence=FunctionalConsequence.MISSENSE,
            gene_name="TESTG",
            allele_frequency=None,  # No gnomAD frequency reported by VEP
            cadd_phred=None,
            transcript_id="ENST00000001",
            hgvsc="c.123A>G",
            hgvsp="p.Lys41Arg",
            ensembl_release="112",
        )
        mock_vep_client.annotate_batch.return_value = [mock_vep_ann]
        engine._vep = mock_vep_client

        mock_clinvar_client = MagicMock()
        mock_clinvar_client.lookup_batch.return_value = [None]
        engine._clinvar = mock_clinvar_client

        variant = _make_variant("chr1", 1000, "A", "G")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        res = results[0]
        assert res.allele_frequency is None
        assert res.frequency_unknown is False
        assert res.frequency_state is ResolutionState.LOOKUP_FAILED

        # Verify PM2 does not fire
        classifier = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()
        classifier._evaluate_pm2(ScoredVariant(annotated=res), tags, missing)

        assert EvidenceTag.PM2 not in tags
        assert "gnomAD" in missing


class TestLocalDatabasesContigScope:
    """Verify local backends produce explicit resolution states."""

    def test_dict_frequency_database_resolutions(self, tmp_path: Path) -> None:
        db = DictFrequencyDatabase()
        unloaded_res = db.lookup_batch_resolutions([("chr1", 1000, "A", "G")])
        assert unloaded_res[0].state is ResolutionState.LOOKUP_FAILED

        tsv_path = tmp_path / "gnomad_test.tsv"
        tsv_path.write_text("chrom\tpos\tref\talt\taf\nchr1\t1000\tA\tG\t0.001\n")
        db.load(tsv_path)

        # Variant 1: present -> FOUND
        # Variant 2: absent -> CONFIRMED_ABSENT
        variants = [
            ("chr1", 1000, "A", "G"),
            ("chr1", 2000, "C", "T"),
        ]
        resolutions = db.lookup_batch_resolutions(variants)

        assert len(resolutions) == 2
        assert resolutions[0].state is ResolutionState.FOUND
        assert resolutions[0].value == 0.001
        assert resolutions[1].state is ResolutionState.CONFIRMED_ABSENT

    @pytest.mark.skipif(not POLARS_AVAILABLE, reason="polars not installed")
    def test_polars_frequency_database_resolutions(self, tmp_path: Path) -> None:
        db = PolarsFrequencyDatabase()
        unloaded_res = db.lookup_batch_resolutions([("chr1", 1000, "A", "G")])
        assert unloaded_res[0].state is ResolutionState.LOOKUP_FAILED

        tsv_path = tmp_path / "gnomad_test.tsv"
        tsv_path.write_text("chrom\tpos\tref\talt\taf\nchr1\t1000\tA\tG\t0.001\n")
        db.load(tsv_path)

        variants = [
            ("chr1", 1000, "A", "G"),
            ("chr1", 2000, "C", "T"),
        ]
        resolutions = db.lookup_batch_resolutions(variants)

        assert len(resolutions) == 2
        assert resolutions[0].state is ResolutionState.FOUND
        assert resolutions[0].value == 0.001
        assert resolutions[1].state is ResolutionState.CONFIRMED_ABSENT

    def test_tabix_frequency_database_contig_scope(self) -> None:
        db = TabixFrequencyDatabase()
        unloaded_res = db.lookup_batch_resolutions([("chr1", 1000, "A", "G")])
        assert unloaded_res[0].state is ResolutionState.LOOKUP_FAILED

        mock_tabix = MagicMock()
        mock_tabix.contigs = ["chr1", "chrX"]
        mock_tabix.fetch.side_effect = lambda chrom, start, end: (
            iter(["chr1\t1000\t.\tA\tG\t.\tPASS\tAF=0.0005"])
            if chrom == "chr1" and start == 999
            else iter([])
        )
        db._tabix = mock_tabix

        variants = [
            ("chr1", 1000, "A", "G"),  # FOUND
            ("chr1", 2000, "C", "T"),  # CONFIRMED_ABSENT
            ("chr2", 3000, "G", "A"),  # LOOKUP_FAILED (not in contigs)
        ]
        resolutions = db.lookup_batch_resolutions(variants)

        assert len(resolutions) == 3
        assert resolutions[0].state is ResolutionState.FOUND
        assert resolutions[0].value == 0.0005

        assert resolutions[1].state is ResolutionState.CONFIRMED_ABSENT

        assert resolutions[2].state is ResolutionState.LOOKUP_FAILED
        assert "chr2" in (resolutions[2].reason or "")


class TestPipelineReattachmentPreservesFrequencyState:
    """Verify Pipeline._reattach_annotations preserves frequency_state."""

    def test_reattach_forwards_frequency_state(self) -> None:
        pipeline = Pipeline.__new__(Pipeline)
        variant = _make_variant("chr1", 1000, "A", "G")
        annotated = AnnotatedVariant(
            variant=variant,
            consequence=FunctionalConsequence.MISSENSE,
            allele_frequency=None,
            frequency_unknown=False,
            frequency_state=ResolutionState.LOOKUP_FAILED,
        )

        reattached = pipeline._reattach_annotations([variant], [annotated])
        assert len(reattached) == 1
        assert reattached[0].frequency_state is ResolutionState.LOOKUP_FAILED
        assert reattached[0].frequency_unknown is False
