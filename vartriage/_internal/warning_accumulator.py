"""Warning accumulator for tracking MissingDataWarning events across a pipeline run.

Collects warnings as they are emitted from annotation and prioritization stages,
tracks counts by source, and emits a summary warning when the total count
exceeds a configurable threshold.
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

from vartriage.models.config import MissingDataConfig
from vartriage.models.warnings import MissingDataWarning

_CONNECTION_FAILURE_REASONS = frozenset({
    "connection_error",
    "timeout",
    "connection_timeout",
    "read_timeout",
    "dns_error",
})

_NOT_FOUND_REASONS = frozenset({
    "not_found",
    None,
})


def is_connection_failure(reason: Optional[str]) -> bool:
    """Determine whether a warning reason indicates a connection/timeout failure.

    Parameters
    ----------
    reason : str or None
        The reason field from a MissingDataWarning.

    Returns
    -------
    bool
        True if the reason represents an infrastructure failure (connection
        error, timeout, etc.) rather than a normal "not found" result.
    """
    if reason is None:
        return False
    return reason in _CONNECTION_FAILURE_REASONS


@dataclass
class MissingDataSummaryWarning(UserWarning):
    """Summary warning emitted when MissingDataWarning count exceeds the threshold.

    Attributes
    ----------
    total_count : int
        Total number of MissingDataWarning events accumulated.
    sources : frozenset[str]
        Names of reference sources that contributed to the count.
    not_found_count : int
        Number of warnings caused by variants simply not found in a database.
    connection_failure_count : int
        Number of warnings caused by connection/timeout failures.
    """

    total_count: int = 0
    sources: frozenset[str] = field(default_factory=frozenset)
    not_found_count: int = 0
    connection_failure_count: int = 0

    def __str__(self) -> str:
        source_list = ", ".join(sorted(self.sources))
        parts = [
            f"Missing data threshold exceeded: {self.total_count} warnings accumulated",
            f"Contributing sources: {source_list}",
            f"Not found: {self.not_found_count}",
            f"Connection/timeout failures: {self.connection_failure_count}",
        ]
        return ". ".join(parts)


class WarningAccumulator:
    """Accumulates MissingDataWarning instances and emits a summary when threshold is exceeded.

    Shared across pipeline stages as a context object. Each stage appends
    warnings, and the accumulator monitors the total count against the
    configured threshold.

    Parameters
    ----------
    config : MissingDataConfig
        Configuration specifying the warning threshold (default 1000).

    Examples
    --------
    >>> from vartriage.models.config import MissingDataConfig
    >>> acc = WarningAccumulator(MissingDataConfig(warning_threshold=2))
    >>> acc.add(MissingDataWarning("chr1", 100, "A", "T", "gnomAD", "not_found"))
    >>> acc.add(MissingDataWarning("chr1", 200, "G", "C", "ClinVar", "not_found"))
    >>> acc.total_count
    2
    >>> acc.threshold_exceeded
    True
    """

    def __init__(self, config: Optional[MissingDataConfig] = None) -> None:
        self._config = config or MissingDataConfig()
        self._warnings: list[MissingDataWarning] = []
        self._count_by_source: defaultdict[str, int] = defaultdict(int)
        self._not_found_count: int = 0
        self._connection_failure_count: int = 0
        self._summary_emitted: bool = False

    @property
    def threshold(self) -> int:
        """The configured warning threshold.

        Returns
        -------
        int
            Maximum warnings before a summary is emitted.
        """
        return self._config.warning_threshold

    @property
    def total_count(self) -> int:
        """Total number of MissingDataWarning events accumulated.

        Returns
        -------
        int
            Current total count.
        """
        return len(self._warnings)

    @property
    def threshold_exceeded(self) -> bool:
        """Whether the total count has exceeded the configured threshold.

        Returns
        -------
        bool
            True if total_count > threshold.
        """
        return self.total_count > self.threshold

    @property
    def sources(self) -> frozenset[str]:
        """Names of all reference sources that contributed warnings.

        Returns
        -------
        frozenset[str]
            Set of source names (e.g., "gnomAD", "ClinVar", "REVEL").
        """
        return frozenset(self._count_by_source.keys())

    @property
    def count_by_source(self) -> dict[str, int]:
        """Warning counts broken down by reference source.

        Returns
        -------
        dict[str, int]
            Mapping of source name to count.
        """
        return dict(self._count_by_source)

    @property
    def not_found_count(self) -> int:
        """Number of warnings where the variant was simply not found in the database.

        Returns
        -------
        int
            Count of "not_found" warnings.
        """
        return self._not_found_count

    @property
    def connection_failure_count(self) -> int:
        """Number of warnings caused by connection errors or timeouts.

        Returns
        -------
        int
            Count of connection/timeout failure warnings.
        """
        return self._connection_failure_count

    @property
    def warnings_list(self) -> list[MissingDataWarning]:
        """All accumulated MissingDataWarning instances.

        Returns
        -------
        list[MissingDataWarning]
            Full list of warnings in accumulation order.
        """
        return list(self._warnings)

    @property
    def summary_emitted(self) -> bool:
        """Whether the summary warning has already been emitted.

        Returns
        -------
        bool
            True if emit_summary_if_exceeded has fired.
        """
        return self._summary_emitted

    def add(self, warning: MissingDataWarning) -> None:
        """Add a MissingDataWarning to the accumulator.

        Tracks the warning, updates per-source counts, classifies the
        reason, and checks whether the threshold has been exceeded. If
        the threshold is exceeded for the first time, emits a summary
        warning via Python's warnings module.

        Parameters
        ----------
        warning : MissingDataWarning
            The warning instance to accumulate.
        """
        self._warnings.append(warning)
        self._count_by_source[warning.source] += 1

        if is_connection_failure(warning.reason):
            self._connection_failure_count += 1
        else:
            self._not_found_count += 1

        if self.threshold_exceeded and not self._summary_emitted:
            self._emit_summary()

    def add_batch(self, batch_warnings: list[MissingDataWarning]) -> None:
        """Add multiple warnings at once.

        Parameters
        ----------
        batch_warnings : list[MissingDataWarning]
            Batch of warnings to accumulate.
        """
        for w in batch_warnings:
            self.add(w)

    def build_summary(self) -> MissingDataSummaryWarning:
        """Build a summary warning object from current state.

        Returns
        -------
        MissingDataSummaryWarning
            Summary containing total count, sources, and failure breakdown.
        """
        return MissingDataSummaryWarning(
            total_count=self.total_count,
            sources=self.sources,
            not_found_count=self._not_found_count,
            connection_failure_count=self._connection_failure_count,
        )

    def reset(self) -> None:
        """Reset all accumulated state.

        Clears all warnings, counts, and the summary-emitted flag.
        """
        self._warnings.clear()
        self._count_by_source.clear()
        self._not_found_count = 0
        self._connection_failure_count = 0
        self._summary_emitted = False

    def _emit_summary(self) -> None:
        """Emit the summary warning via Python's warnings module."""
        summary = self.build_summary()
        warnings.warn(summary, stacklevel=3)
        self._summary_emitted = True
