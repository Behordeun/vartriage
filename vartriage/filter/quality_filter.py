"""Quality-based variant filtering.

Excludes variants that fail machine-level quality controls based on the
FILTER field value and the Phred-scaled QUAL score.
"""

from __future__ import annotations

import warnings
from typing import Iterator

from vartriage.models.config import QualityFilterConfig
from vartriage.models.variant import Variant
from vartriage.models.warnings import MissingDataWarning

_PASSING_FILTER_VALUES = frozenset({"PASS", "."})


class QualityFilter:
    """Exclude variants that fail quality controls.

    Applies two sequential checks to each variant:

    1. The FILTER field must be ``"PASS"`` or ``"."`` (missing/not-applied).
    2. The QUAL score must be present and at least ``min_qual``.

    Variants that fail either check are silently dropped from the output
    stream, except when QUAL is missing. In that case a
    ``MissingDataWarning`` is emitted via Python's ``warnings`` module before
    the variant is excluded.

    Parameters
    ----------
    config : QualityFilterConfig
        Filtering parameters. The ``min_qual`` field specifies the minimum
        acceptable QUAL score (default 20.0, valid range [0, 1_000_000]).

    Raises
    ------
    ValueError
        If ``config.min_qual`` is outside [0, 1_000_000]. This is enforced by
        ``QualityFilterConfig.__post_init__`` at config construction time.

    Examples
    --------
    >>> from vartriage.models.config import QualityFilterConfig
    >>> from vartriage.models.variant import Variant
    >>> cfg = QualityFilterConfig(min_qual=30.0)
    >>> qf = QualityFilter(cfg)
    >>> variants = [
    ...     Variant(chrom="chr1", pos=100, id=None, ref="A", alt="T",
    ...             qual=50.0, filter_status="PASS"),
    ...     Variant(chrom="chr1", pos=200, id=None, ref="G", alt="C",
    ...             qual=10.0, filter_status="PASS"),
    ... ]
    >>> list(qf.apply(iter(variants)))  # only first variant passes
    [Variant(chrom='chr1', pos=100, ...)]
    """

    def __init__(self, config: QualityFilterConfig | None = None) -> None:
        if config is None:
            config = QualityFilterConfig()
        self._min_qual = config.min_qual

    def apply(self, variants: Iterator[Variant]) -> Iterator[Variant]:
        """Filter variants by FILTER field and QUAL score.

        Parameters
        ----------
        variants : Iterator[Variant]
            Input stream of variant records. May be empty.

        Yields
        ------
        Variant
            Variants that pass both the FILTER and QUAL checks, in the same
            relative order they appeared in the input.

        Warns
        -----
        MissingDataWarning
            When a variant has ``qual=None``. The warning message includes
            the chromosome and position to aid debugging.
        """
        min_qual = self._min_qual

        for variant in variants:
            if variant.filter_status not in _PASSING_FILTER_VALUES:
                continue

            if variant.qual is None:
                reason = (
                    f"Missing QUAL score for variant at "
                    f"{variant.chrom}:{variant.pos}"
                )
                warning_data = MissingDataWarning(
                    chrom=variant.chrom,
                    pos=variant.pos,
                    ref=variant.ref,
                    alt=variant.alt,
                    source="QUAL",
                    reason=reason,
                )
                warnings.warn(
                    UserWarning(warning_data),
                    stacklevel=2,
                )
                continue

            if variant.qual < min_qual:
                continue

            yield variant
