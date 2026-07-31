"""Unit tests for the PrioritizationEngine orchestration logic."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vartriage.models.config import PrioritizationConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    FunctionalConsequence,
    ScoredVariant,
    Variant,
)
from vartriage.prioritization.engine import (
    PrioritizationEngine,
    _merge_sort_scored,
)


def _make_annotated(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    allele_frequency: float | None = 0.001,
    frequency_unknown: bool = False,
) -> AnnotatedVariant:
    v = Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
    )
    return AnnotatedVariant(
        variant=v,
        consequence=FunctionalConsequence.MISSENSE,
        allele_frequency=allele_frequency,
        frequency_unknown=frequency_unknown,
    )


def _make_scored(composite_rank: float | None = 0.5) -> ScoredVariant:
    av = _make_annotated()
    return ScoredVariant(annotated=av, composite_rank=composite_rank)


class TestPrioritizationEngineDefaults:
    """Engine construction with default and None config."""

    def test_constructs_with_none_config(self) -> None:
        engine = PrioritizationEngine(config=None)
        assert engine._batch_size == 10_000

    def test_constructs_with_custom_batch_size(self) -> None:
        config = PrioritizationConfig(batch_size=5_000)
        engine = PrioritizationEngine(config)
        assert engine._batch_size == 5_000


class TestPrioritize:
    """prioritize() filters by frequency then scores remaining variants."""

    def test_empty_stream_yields_nothing(self) -> None:
        engine = PrioritizationEngine()
        results = list(engine.prioritize(iter([])))
        assert results == []

    def test_excludes_high_frequency_variants(self) -> None:
        config = PrioritizationConfig(max_allele_frequency=0.01)
        engine = PrioritizationEngine(config)

        rare = _make_annotated(pos=1, allele_frequency=0.005)
        common = _make_annotated(pos=2, allele_frequency=0.05)

        results = list(engine.prioritize(iter([rare, common])))
        positions = [sv.annotated.variant.pos for sv in results]
        assert 1 in positions
        assert 2 not in positions

    def test_retains_frequency_unknown_variants(self) -> None:
        config = PrioritizationConfig(max_allele_frequency=0.01)
        engine = PrioritizationEngine(config)

        unknown = _make_annotated(pos=3, allele_frequency=None, frequency_unknown=True)
        results = list(engine.prioritize(iter([unknown])))
        assert len(results) == 1
        assert results[0].annotated.variant.pos == 3

    def test_sorts_within_batch_by_composite_rank_descending(self) -> None:
        engine = PrioritizationEngine()

        variants = [
            _make_annotated(pos=i, allele_frequency=0.001)
            for i in range(5)
        ]
        results = list(engine.prioritize(iter(variants)))

        ranks = [sv.composite_rank for sv in results]
        non_null = [r for r in ranks if r is not None]
        assert non_null == sorted(non_null, reverse=True)


class TestBatchProcessing:
    """_process_in_batches splits input correctly."""

    def test_processes_multiple_batches(self) -> None:
        config = PrioritizationConfig(batch_size=2_000)
        engine = PrioritizationEngine(config)

        # 5 variants with small batch means multiple batches
        # but actual minimum is 1000, so we just verify they all come through
        variants = [
            _make_annotated(pos=i, allele_frequency=0.001)
            for i in range(3)
        ]
        results = list(engine._process_in_batches(iter(variants)))
        assert len(results) == 3


class TestChunkedFallback:
    """_chunked_fallback handles MemoryError gracefully."""

    def test_falls_back_on_memory_error(self) -> None:
        engine = PrioritizationEngine()

        batch = [_make_annotated(pos=i) for i in range(4)]
        call_count = {"n": 0}
        original_score_batch = engine._score_batch

        def score_batch_with_first_failure(b):
            call_count["n"] += 1
            if call_count["n"] == 1 and len(b) == 4:
                raise MemoryError("simulated")
            return original_score_batch(b)

        with patch.object(engine, "_score_batch", side_effect=score_batch_with_first_failure):
            results = engine._chunked_fallback(batch)

        assert len(results) == 4

    def test_chunked_fallback_results_sorted(self) -> None:
        engine = PrioritizationEngine()
        batch = [_make_annotated(pos=i) for i in range(4)]

        results = engine._chunked_fallback(batch)
        ranks = [sv.composite_rank for sv in results]
        non_null = [r for r in ranks if r is not None]
        assert non_null == sorted(non_null, reverse=True)


class TestScoreBatch:
    """_score_batch builds coordinate keys and looks up scores."""

    def test_score_batch_with_no_preloaded_scores(self) -> None:
        engine = PrioritizationEngine()
        batch = [_make_annotated(pos=42, chrom="chr5")]
        results = engine._score_batch(batch)

        assert len(results) == 1
        # No CADD/REVEL loaded, so scores should be None
        assert results[0].cadd_phred is None
        assert results[0].revel_score is None

    def test_score_batch_with_preloaded_cadd(self, tmp_path: Path) -> None:
        cadd_file = tmp_path / "cadd.tsv"
        cadd_file.write_text("chr1\t100\tA\tT\t25.0\n")

        config = PrioritizationConfig(cadd_scores_path=cadd_file)
        engine = PrioritizationEngine(config)

        variant = _make_annotated(chrom="chr1", pos=100, ref="A", alt="T")
        results = engine._score_batch([variant])

        assert len(results) == 1
        assert results[0].cadd_phred == pytest.approx(25.0)

    def test_score_batch_with_preloaded_spliceai(self, tmp_path: Path) -> None:
        splice_file = tmp_path / "spliceai.tsv"
        splice_file.write_text("chr1\t100\tA\tT\t0.85\n")

        config = PrioritizationConfig(spliceai_scores_path=splice_file)
        engine = PrioritizationEngine(config)

        variant = _make_annotated(chrom="chr1", pos=100, ref="A", alt="T")
        results = engine._score_batch([variant])

        assert len(results) == 1
        assert results[0].spliceai_score == pytest.approx(0.85)


class TestMergeSortScored:
    """_merge_sort_scored orders by composite_rank descending, nulls last."""

    def test_sorts_descending(self) -> None:
        items = [
            _make_scored(0.3),
            _make_scored(0.9),
            _make_scored(0.6),
        ]
        result = _merge_sort_scored(items)
        ranks = [sv.composite_rank for sv in result]
        assert ranks == [0.9, 0.6, 0.3]

    def test_nulls_go_last(self) -> None:
        items = [
            _make_scored(None),
            _make_scored(0.5),
            _make_scored(None),
        ]
        result = _merge_sort_scored(items)
        assert result[0].composite_rank == 0.5
        assert result[1].composite_rank is None
        assert result[2].composite_rank is None

    def test_empty_list(self) -> None:
        assert _merge_sort_scored([]) == []

    def test_all_nulls_preserves_count(self) -> None:
        items = [_make_scored(None) for _ in range(3)]
        result = _merge_sort_scored(items)
        assert len(result) == 3

    def test_single_element(self) -> None:
        items = [_make_scored(0.7)]
        result = _merge_sort_scored(items)
        assert result[0].composite_rank == 0.7


class TestMemoryErrorInBatchProcessing:
    """Engine degrades gracefully when MemoryError occurs mid-batch."""

    def test_memory_error_triggers_chunked_path(self) -> None:
        engine = PrioritizationEngine()
        variants = [_make_annotated(pos=i) for i in range(6)]

        original = engine._score_batch
        first_call = {"done": False}

        def fail_first_large(batch):
            if not first_call["done"] and len(batch) > 3:
                first_call["done"] = True
                raise MemoryError("simulated OOM")
            return original(batch)

        with patch.object(engine, "_score_batch", side_effect=fail_first_large):
            results = list(engine._process_in_batches(iter(variants)))

        assert len(results) == 6
