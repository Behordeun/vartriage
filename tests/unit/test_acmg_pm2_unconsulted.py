"""Tests ensuring unconsulted frequency sources never trigger ACMG PM2.

In clinical genomics and ACMG/AMP variant interpretation, PM2 encodes
evidence of rarity or confirmed absence from population cohorts (e.g. gnomAD).
The absence of evidence (an unconsulted reference database, a network failure,
or an unannotated run) must not be conflated with evidence of absence.
"""

from __future__ import annotations

from pathlib import Path

from vartriage.annotation.engine import AnnotationEngine
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
