"""Allele frequency-based variant filtering.

Excludes variants whose population allele frequency exceeds a configurable
threshold, while retaining variants with unknown frequency data to prevent
silent loss of potentially pathogenic candidates.
"""

from __future__ import annotations

from typing import Iterator

from vartriage.models.config import PrioritizationConfig
from vartriage.models.variant import AnnotatedVariant


class FrequencyFilter:
    """Exclude variants with allele frequency above a maximum threshold.

    Variants are retained when any of the following hold:

    1. Their ``frequency_unknown`` flag is True (absent from gnomAD).
    2. Their ``allele_frequency`` is None (implies unknown).
    3. Their ``allele_frequency`` is less than or equal to the configured
       maximum threshold.

    The input ordering of passing variants is preserved in the output.

    Parameters
    ----------
    config : PrioritizationConfig, optional
        Configuration containing the ``max_allele_frequency`` threshold.
        When None, a default config with ``max_allele_frequency=0.01`` is
        used.

    Raises
    ------
    ValueError
        If ``config.max_allele_frequency`` is outside [0.0, 1.0]. This is
        enforced by ``PrioritizationConfig.__post_init__`` at config
        construction time.

    Examples
    --------
    >>> from vartriage.models.config import PrioritizationConfig
    >>> from vartriage.models.variant import (
    ...     AnnotatedVariant, FunctionalConsequence, Variant,
    ... )
    >>> cfg = PrioritizationConfig(max_allele_frequency=0.01)
    >>> ff = FrequencyFilter(cfg)
    >>> v = Variant(chrom="chr1", pos=100, id=None, ref="A", alt="T",
    ...            qual=30.0, filter_status="PASS")
    >>> rare = AnnotatedVariant(variant=v,
    ...     consequence=FunctionalConsequence.MISSENSE,
    ...     allele_frequency=0.005)
    >>> list(ff.apply(iter([rare])))
    [AnnotatedVariant(...)]
    """

    def __init__(self, config: PrioritizationConfig | None = None) -> None:
        if config is None:
            config = PrioritizationConfig()
        self._max_af = config.max_allele_frequency

    def apply(self, variants: Iterator[AnnotatedVariant]) -> Iterator[AnnotatedVariant]:
        """Filter variants by allele frequency threshold.

        Parameters
        ----------
        variants : Iterator[AnnotatedVariant]
            Input stream of annotated variant records. May be empty.

        Yields
        ------
        AnnotatedVariant
            Variants that pass the frequency filter, in the same relative
            order they appeared in the input.
        """
        max_af = self._max_af

        for variant in variants:
            if variant.frequency_unknown:
                yield variant
                continue

            if variant.allele_frequency is None:
                yield variant
                continue

            if variant.allele_frequency <= max_af:
                yield variant
