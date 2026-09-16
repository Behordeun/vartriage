"""The annotation engine populates per-population allele frequencies.

A frequency backend that can return ancestry-group frequencies should have
them threaded onto the annotated variant, so BA1/BS1/PM2 can reason about a
variant that is common in one population but globally rare, instead of
seeing only a single global number.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.annotation.engine import AnnotationEngine
from vartriage.models.config import AnnotationConfig
from vartriage.models.variant import Variant


class _GlobalOnlyDB:
    def load(self, reference_path):  # noqa: ANN001, ANN201
        pass

    def lookup_batch(self, variants):  # noqa: ANN001, ANN201
        return [0.001 for _ in variants]


class _PerPopulationDB(_GlobalOnlyDB):
    def lookup_batch_populations(self, variants):  # noqa: ANN001, ANN201
        # A variant common in AFR but globally rare.
        return [{"AF": 0.001, "AF_afr": 0.08, "AF_nfe": 0.0005} for _ in variants]


@pytest.fixture
def gene_annotation(tmp_path: Path) -> Path:
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("")
    return gtf


def _engine_with_db(db, gene_annotation: Path) -> AnnotationEngine:  # noqa: ANN001
    engine = AnnotationEngine(AnnotationConfig(gene_annotation_path=gene_annotation))
    engine.set_frequency_db(db)
    return engine


def _variant() -> Variant:
    return Variant(
        chrom="chr1",
        pos=100,
        id=None,
        ref="A",
        alt="G",
        qual=None,
        filter_status="PASS",
        info={},
    )


class TestPerPopulationWiring:
    def test_per_population_frequencies_are_populated(self, gene_annotation) -> None:  # noqa: ANN001
        engine = _engine_with_db(_PerPopulationDB(), gene_annotation)
        annotated = list(engine.annotate(iter([_variant()])))[0]
        pf = annotated.population_frequencies
        assert pf is not None
        assert pf.afr == 0.08
        assert pf.nfe == 0.0005
        assert pf.global_af == 0.001

    def test_global_only_backend_leaves_population_frequencies_none(
        self, gene_annotation
    ) -> None:  # noqa: ANN001
        engine = _engine_with_db(_GlobalOnlyDB(), gene_annotation)
        annotated = list(engine.annotate(iter([_variant()])))[0]
        assert annotated.population_frequencies is None
        assert annotated.allele_frequency == 0.001

    def test_max_population_af_reflects_the_ancestry_spike(
        self, gene_annotation
    ) -> None:  # noqa: ANN001
        engine = _engine_with_db(_PerPopulationDB(), gene_annotation)
        annotated = list(engine.annotate(iter([_variant()])))[0]
        pf = annotated.population_frequencies
        assert pf is not None
        assert pf.max_population_af == 0.08
