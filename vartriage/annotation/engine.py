"""Annotation engine with backend auto-detection.

Composes consequence assignment, gnomAD frequency lookups, and ClinVar
lookups into a batch-processing pipeline. Picks the fastest available
backend (pyranges, polars) at init time, falling back to pure-Python
when the optional deps aren't installed.
"""

from __future__ import annotations

import logging
from itertools import islice
from pathlib import Path
from typing import Iterator, Optional

from vartriage.models.config import AnnotationConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    ClinVarAssertion,
    Variant,
)
from vartriage.models.warnings import MissingDataWarning
from vartriage.protocols import (
    ClinVarDatabase,
    FrequencyDatabase,
    IntervalIndex,
)

logger = logging.getLogger(__name__)


def _pyranges_available() -> bool:
    """Check whether pyranges is importable."""
    try:
        import pyranges  # noqa: F401

        return True
    except ImportError:
        return False


def _polars_available() -> bool:
    """Check whether polars is importable."""
    try:
        import polars  # noqa: F401

        return True
    except ImportError:
        return False


class AnnotationEngine:
    """Annotates variants with consequence, frequency, and ClinVar data.

    Processes in configurable batch sizes. Picks the fastest backend
    at construction time:

    - Consequence: pyranges if available, else pure-Python sorted intervals
    - Frequency/ClinVar: polars if available, else dict-based lookups

    Parameters
    ----------
    config : AnnotationConfig
        Reference file paths and batch_size.

    Raises
    ------
    FileNotFoundError
        If a required reference file is missing.
    ValueError
        If batch_size is outside [1_000, 100_000].
    """

    def __init__(self, config: AnnotationConfig) -> None:
        self._config = config
        self._warnings: list[MissingDataWarning] = []

        # Validate file paths upfront (fail-fast)
        self._validate_paths(config)

        # Initialize consequence annotator
        self._consequence_annotator: IntervalIndex = (
            self._build_consequence_annotator(
                config.gene_annotation_path
            )
        )

        # Initialize frequency database
        self._frequency_db: FrequencyDatabase = (
            self._build_frequency_db(config.gnomad_path)
        )

        # Initialize ClinVar database (optional)
        self._clinvar_db: Optional[ClinVarDatabase] = (
            self._build_clinvar_db(config.clinvar_path)
        )

    @property
    def warnings(self) -> list[MissingDataWarning]:
        """Warnings accumulated for variants missing from references."""
        return self._warnings

    def annotate(
        self, variants: Iterator[Variant]
    ) -> Iterator[AnnotatedVariant]:
        """Annotate variants with consequence, frequency, and ClinVar data.

        Processes variants in batches of ``config.batch_size`` (default
        10,000). Each batch undergoes consequence assignment, frequency
        lookup, and ClinVar lookup via the Protocol interfaces selected
        at construction time.

        Parameters
        ----------
        variants : Iterator[Variant]
            Input stream of raw variants.

        Yields
        ------
        AnnotatedVariant
            Variants enriched with functional consequence, allele
            frequency, and ClinVar clinical significance.
        """
        batch_size = self._config.batch_size

        while True:
            batch = list(islice(variants, batch_size))
            if not batch:
                break

            yield from self._annotate_batch(batch)

    def _annotate_batch(
        self, batch: list[Variant]
    ) -> list[AnnotatedVariant]:
        """Run consequence + frequency + ClinVar on a single batch."""
        # Consequence assignment
        consequences = self._consequence_annotator.assign_batch(batch)

        # Frequency lookup
        variant_keys = [
            (v.chrom, v.pos, v.ref, v.alt) for v in batch
        ]
        frequencies = self._frequency_db.lookup_batch(variant_keys)

        # ClinVar lookup
        clinvar_assertions: list[Optional[ClinVarAssertion]] = []
        if self._clinvar_db is not None:
            clinvar_assertions = self._clinvar_db.lookup_batch(variant_keys)
        else:
            clinvar_assertions = [None] * len(batch)

        # Compose results
        results: list[AnnotatedVariant] = []
        for i, variant in enumerate(batch):
            freq = frequencies[i]
            clinvar = clinvar_assertions[i]

            frequency_unknown = freq is None
            clinvar_unknown = clinvar is None

            # Emit warnings for missing data
            if frequency_unknown:
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="gnomAD",
                        reason="not_found",
                    )
                )

            if clinvar_unknown and self._clinvar_db is not None:
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="ClinVar",
                        reason="not_found",
                    )
                )

            results.append(
                AnnotatedVariant(
                    variant=variant,
                    consequence=consequences[i],
                    allele_frequency=freq,
                    clinvar_assertion=clinvar,
                    frequency_unknown=frequency_unknown,
                    clinvar_unknown=clinvar_unknown,
                )
            )

        return results

    def _validate_paths(self, config: AnnotationConfig) -> None:
        """Fail fast if any required reference file is missing."""
        if not config.gene_annotation_path.exists():
            raise FileNotFoundError(
                f"Gene annotation file not found: "
                f"{config.gene_annotation_path}"
            )

        if not config.gnomad_path.exists():
            raise FileNotFoundError(
                f"gnomAD reference file not found: {config.gnomad_path}"
            )

        if config.clinvar_path is not None and not config.clinvar_path.exists():
            raise FileNotFoundError(
                f"ClinVar reference file not found: {config.clinvar_path}"
            )

    def _build_consequence_annotator(
        self, annotation_path: Path
    ) -> IntervalIndex:
        """Pick the best consequence annotator available.

        Tries pyranges first, falls back to the pure-Python interval tree.

        Parameters
        ----------
        annotation_path : Path
            GTF/GFF gene annotation file.

        Returns
        -------
        IntervalIndex
            Loaded consequence annotator.
        """
        if _pyranges_available():
            try:
                from vartriage.annotation.consequence_pyranges import (
                    PyRangesConsequenceAnnotator,
                )

                logger.info(
                    "Using pyranges backend for consequence annotation"
                )
                return PyRangesConsequenceAnnotator(annotation_path)
            except Exception as exc:
                logger.warning(
                    "pyranges backend failed, falling back to pure-Python: %s",
                    exc,
                )

        from vartriage.annotation.consequence import (
            ConsequenceAnnotator,
        )

        logger.info(
            "Using pure-Python backend for consequence annotation"
        )
        return ConsequenceAnnotator(annotation_path)

    def _build_frequency_db(
        self, gnomad_path: Path
    ) -> FrequencyDatabase:
        """Pick the best frequency database available.

        Tries polars first, falls back to the pure-Python dict.

        Parameters
        ----------
        gnomad_path : Path
            gnomAD reference TSV.

        Returns
        -------
        FrequencyDatabase
            Loaded frequency database.
        """
        if _polars_available():
            try:
                from vartriage.annotation.frequency_polars import (
                    PolarsFrequencyDatabase,
                )

                logger.info(
                    "Using polars backend for frequency lookup"
                )
                freq_db: FrequencyDatabase = PolarsFrequencyDatabase()
                freq_db.load(gnomad_path)
                return freq_db
            except Exception as exc:
                logger.warning(
                    "polars frequency backend failed, falling back to "
                    "pure-Python: %s",
                    exc,
                )

        from vartriage.annotation.frequency import (
            DictFrequencyDatabase,
        )

        logger.info("Using pure-Python backend for frequency lookup")
        freq_db = DictFrequencyDatabase()
        freq_db.load(gnomad_path)
        return freq_db

    def _build_clinvar_db(
        self, clinvar_path: Optional[Path]
    ) -> Optional[ClinVarDatabase]:
        """Pick the best ClinVar backend, or None if no path given.

        Parameters
        ----------
        clinvar_path : Optional[Path]
            ClinVar reference file, or None to skip.

        Returns
        -------
        Optional[ClinVarDatabase]
            Loaded ClinVar database, or None.
        """
        if clinvar_path is None:
            return None

        if _polars_available():
            try:
                from vartriage.annotation.clinvar_polars import (
                    PolarsClinVarDatabase,
                )

                logger.info(
                    "Using polars backend for ClinVar lookup"
                )
                clinvar_db: ClinVarDatabase = PolarsClinVarDatabase()
                clinvar_db.load(clinvar_path)
                return clinvar_db
            except Exception as exc:
                logger.warning(
                    "polars ClinVar backend failed, falling back to "
                    "pure-Python: %s",
                    exc,
                )

        from vartriage.annotation.clinvar import (
            DictClinVarDatabase,
        )

        logger.info("Using pure-Python backend for ClinVar lookup")
        clinvar_db = DictClinVarDatabase()
        clinvar_db.load(clinvar_path)
        return clinvar_db
