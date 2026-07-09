"""Unit tests for batch iterator utilities."""

from __future__ import annotations

import pytest

from vartriage._internal.batch import (
    _MIN_CHUNK_SIZE,
    batched,
    process_with_memory_fallback,
)


class TestBatched:
    """Tests for the batched() chunking utility."""

    def test_basic_chunking(self) -> None:
        result = list(batched(range(7), 3))
        assert result == [[0, 1, 2], [3, 4, 5], [6]]

    def test_exact_multiple(self) -> None:
        result = list(batched(range(6), 3))
        assert result == [[0, 1, 2], [3, 4, 5]]

    def test_empty_iterable(self) -> None:
        result = list(batched([], 5))
        assert result == []

    def test_batch_size_larger_than_iterable(self) -> None:
        result = list(batched([1, 2, 3], 10))
        assert result == [[1, 2, 3]]

    def test_batch_size_one(self) -> None:
        result = list(batched([10, 20, 30], 1))
        assert result == [[10], [20], [30]]

    def test_single_item(self) -> None:
        result = list(batched([42], 5))
        assert result == [[42]]

    def test_works_with_iterator_input(self) -> None:
        gen = (x * 2 for x in range(5))
        result = list(batched(gen, 2))
        assert result == [[0, 2], [4, 6], [8]]

    def test_invalid_batch_size_zero(self) -> None:
        with pytest.raises(ValueError, match="batch_size must be >= 1"):
            list(batched(range(10), 0))

    def test_invalid_batch_size_negative(self) -> None:
        with pytest.raises(ValueError, match="batch_size must be >= 1"):
            list(batched(range(10), -3))


class TestProcessWithMemoryFallback:
    """Tests for the MemoryError recovery logic."""

    def test_no_memory_error_processes_normally(self) -> None:
        items = [1, 2, 3, 4, 5]
        result = process_with_memory_fallback(
            items, lambda xs: [x * 2 for x in xs]
        )
        assert result == [2, 4, 6, 8, 10]

    def test_memory_error_triggers_chunked_fallback(self) -> None:
        call_count = {"n": 0}

        def processor(xs: list[int]) -> list[int]:
            call_count["n"] += 1
            if call_count["n"] == 1 and len(xs) > 5:
                raise MemoryError("simulated OOM")
            return [x * 10 for x in xs]

        items = list(range(10))
        result = process_with_memory_fallback(
            items, processor, initial_chunk_size=5
        )
        assert sorted(result) == [x * 10 for x in range(10)]

    def test_empty_items(self) -> None:
        result = process_with_memory_fallback([], lambda xs: xs)
        assert result == []

    def test_halves_chunk_size_on_repeated_memory_error(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """When chunk_size is above the minimum and fails, it halves."""
        import vartriage._internal.batch as batch_mod

        monkeypatch.setattr(batch_mod, "_MIN_CHUNK_SIZE", 3)

        def processor(xs: list[int]) -> list[int]:
            if len(xs) > 4:
                raise MemoryError("too large")
            return [x + 1 for x in xs]

        items = list(range(12))
        # initial_chunk_size=8, fails (>4), halves to max(3, 4)=4, succeeds
        result = process_with_memory_fallback(
            items, processor, initial_chunk_size=8
        )
        assert sorted(result) == [x + 1 for x in range(12)]

    def test_reraises_at_minimum_chunk_size(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import vartriage._internal.batch as batch_mod

        monkeypatch.setattr(batch_mod, "_MIN_CHUNK_SIZE", 3)

        def always_oom(xs: list[int]) -> list[int]:
            raise MemoryError("cannot allocate")

        items = list(range(10))
        with pytest.raises(MemoryError, match="cannot allocate"):
            process_with_memory_fallback(
                items, always_oom, initial_chunk_size=3
            )
