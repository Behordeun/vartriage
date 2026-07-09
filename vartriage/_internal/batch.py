"""Batch iterator utilities for memory-bounded pipeline processing.

Provides chunking helpers used throughout the pipeline to process variants
in configurable batches, keeping peak memory usage bounded regardless of
input size. Includes MemoryError recovery logic that halves the chunk size
(minimum 500,000) on failure and retries.
"""

from __future__ import annotations

import logging
from itertools import islice
from typing import Callable, Iterable, Iterator, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")
R = TypeVar("R")

_MIN_CHUNK_SIZE: int = 500_000


def batched(iterable: Iterable[T], batch_size: int) -> Iterator[list[T]]:
    """Yield successive lists of at most ``batch_size`` items from an iterable.

    Consumes the iterable lazily — only materializes one batch at a time.
    Yields nothing for an empty iterable. The final batch may contain fewer
    than ``batch_size`` items.

    Parameters
    ----------
    iterable : Iterable[T]
        Source of items to chunk. Can be an iterator or any iterable.
    batch_size : int
        Maximum number of items per batch. Must be at least 1.

    Yields
    ------
    list[T]
        A list of up to ``batch_size`` items from the source iterable.

    Raises
    ------
    ValueError
        If ``batch_size`` is less than 1.

    Examples
    --------
    >>> list(batched(range(7), 3))
    [[0, 1, 2], [3, 4, 5], [6]]

    >>> list(batched([], 10))
    []
    """
    if batch_size < 1:
        raise ValueError(f"batch_size must be >= 1, got {batch_size}")

    it = iter(iterable)
    while True:
        batch = list(islice(it, batch_size))
        if not batch:
            return
        yield batch


def process_with_memory_fallback(
    items: list[T],
    processor: Callable[[list[T]], list[R]],
    initial_chunk_size: int | None = None,
) -> list[R]:
    """Process a list of items with automatic MemoryError recovery.

    Attempts to process the full list in one pass. On MemoryError, splits
    the work into smaller chunks, halving the size on each retry down to a
    minimum of ``_MIN_CHUNK_SIZE`` (500,000). This prevents the pipeline
    from failing entirely when vectorized operations hit memory limits.

    Parameters
    ----------
    items : list[T]
        The items to process.
    processor : Callable[[list[T]], list[R]]
        A function that processes a list of items and returns results.
        Must be safe to call on sublists of the original.
    initial_chunk_size : int, optional
        Starting chunk size for fallback processing. Defaults to half
        the length of ``items`` (minimum ``_MIN_CHUNK_SIZE``).

    Returns
    -------
    list[R]
        Combined results from processing all chunks.

    Raises
    ------
    MemoryError
        Re-raised if processing fails even at the minimum chunk size.
    """
    try:
        return processor(items)
    except MemoryError:
        logger.warning(
            "MemoryError processing %d items. Falling back to chunked mode.",
            len(items),
        )

    chunk_size = initial_chunk_size or max(_MIN_CHUNK_SIZE, len(items) // 2)
    results: list[R] = []

    start = 0
    while start < len(items):
        chunk = items[start:start + chunk_size]
        try:
            results.extend(processor(chunk))
            start += chunk_size
        except MemoryError:
            if chunk_size <= _MIN_CHUNK_SIZE:
                logger.error(
                    "MemoryError persists at minimum chunk size (%d). "
                    "Re-raising.",
                    _MIN_CHUNK_SIZE,
                )
                raise
            chunk_size = max(_MIN_CHUNK_SIZE, chunk_size // 2)
            logger.warning(
                "MemoryError at chunk_size=%d, reducing to %d.",
                chunk_size * 2,
                chunk_size,
            )

    return results
