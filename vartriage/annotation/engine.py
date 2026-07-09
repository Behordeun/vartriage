"""Annotation engine orchestrator with backend auto-detection.

Composes functional consequence assignment, population frequency lookups,
and ClinVar assertion lookups into a single batch-processing pipeline.
Auto-detects available backends (pyranges, polars) at construction time
and selects the fastest available, falling back to pure-Python
implementations when optional dependencies are not installed.
"""

from __future__ import annotations

import logging
from itertools import islice
from pathlib import Path
from typing import Any, Iterator, Optional

from vartriage.models.config import AnnotationConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    ClinVarAssertion,
    Variant,
)
from vartriage.models.warnings import MissingDataWarning

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
    """Annotate variants with functional consequence and population data.

    Orchestrates consequence assignment, gnomAD frequency lookup, and
    ClinVar assertion lookup in configurable batch sizes. Auto-detects
    the fastest available backend at construction time:

    - Consequence: PyRangesIntervalIndex if pyranges installed, else
      SortedArrayIntervalIndex (pure-Python)
    - Frequency: PolarsFrequencyDatabase if polars installed, else
      DictFrequencyDatabase (pure-Python)
    - ClinVar: PolarsClinVarDatabase if polars installed, else
      DictClinVarDatabase (pure-Python)

    Parameters
    ----------
    config : AnnotationConfig
        Paths to gene annotation (GTF/GFF), gnomAD, and ClinVar reference
        files, plus batch_size for vectorized operations.

    Raises
    ------
    FileNotFoundError
        If any required reference file path does not exist.
    ReferenceFileError
        If a reference file cannot be read or parsed.
    ValueError
        If batch_size is outside [1_000, 100_000].
    """

    def __init__(self, config: AnnotationConfig) -> None:
        self._config = config
        self._warnings: list[MissingDataWarning] = []

        # Validate file paths upfront (fail-fast)
        self._validate_paths(config)

        # Initialize consequence annotator
        self._consequence_annotator = self._build_consequence_annotator(
            config.gene_annotation_path
        )

        # Initialize frequency database
        self._frequency_db = self._build_frequency_db(config.gnomad_path)

        # Initialize ClinVar database (optional)
        self._clinvar_db = self._build_clinvar_db(config.clinvar_path)

    @property
    def warnings(self) -> list[MissingDataWarning]:
        """All MissingDataWarning instances accumulated during annotation.

        Returns
        -------
        list[MissingDataWarning]
            Warnings emitted for variants missing from reference databases.
        """
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
        """Annotate a single batch of variants.

        Parameters
        ----------
        batch : list[Variant]
            Batch of variants to annotate.

        Returns
        -------
        list[AnnotatedVariant]
            Annotated variants in the same order as input.
        """
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
        """Validate all reference file paths exist at initialization.

        Parameters
        ----------
        config : AnnotationConfig
            Configuration containing file paths to validate.

        Raises
        ------
        FileNotFoundError
            If any required reference file does not exist.
        """
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

    def _build_consequence_annotator(self, annotation_path: Path) -> Any:
        """Build the best available consequence annotator.

        Attempts to use the pyranges-based backend first; falls back
        to the pure-Python sorted interval tree if pyranges is not
        installed.

        Parameters
        ----------
        annotation_path : Path
            Path to the GTF/GFF gene annotation file.

        Returns
        -------
        ConsequenceAnnotator or PyRangesConsequenceAnnotator
            The consequence annotator instance with loaded data.
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

    def _build_frequency_db(self, gnomad_path: Path) -> Any:
        """Build the best available frequency database.

        Attempts to use the polars-based backend first; falls back to
        the pure-Python dict-based implementation if polars is not
        installed.

        Parameters
        ----------
        gnomad_path : Path
            Path to the gnomAD reference TSV file.

        Returns
        -------
        DictFrequencyDatabase or PolarsFrequencyDatabase
            The frequency database instance with loaded data.
        """
        if _polars_available():
            try:
                from vartriage.annotation.frequency_polars import (
                    PolarsFrequencyDatabase,
                )

                logger.info(
                    "Using polars backend for frequency lookup"
                )
                freq_db: Any = PolarsFrequencyDatabase()
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

    def _build_clinvar_db(self, clinvar_path: Optional[Path]) -> Any:
        """Build the best available ClinVar database, or None if path is absent.

        Parameters
        ----------
        clinvar_path : Optional[Path]
            Path to the ClinVar reference file, or None to skip ClinVar.

        Returns
        -------
        DictClinVarDatabase or PolarsClinVarDatabase or None
            The ClinVar database instance with loaded data, or None.
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
                clinvar_db: Any = PolarsClinVarDatabase()
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
