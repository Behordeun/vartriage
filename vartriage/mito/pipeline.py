"""Mitochondrial variant analysis sub-pipeline.

Processes chrM/MT variants through the mtDNA-specific pathway:
heteroplasmy extraction -> gene map annotation -> MITOMAP/HelixMTdb
lookup -> mitochondrial classification. Runs separately from the
nuclear ACMG pipeline and its results merge at the report stage.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator

from vartriage.mito.classifier import (
    MitochondrialClassifier,
    MitoClassifiedVariant,
)
from vartriage.mito.config import MitoConfig
from vartriage.mito.frequency import HelixMTdbDatabase
from vartriage.mito.gene_map import MtGeneMap
from vartriage.mito.mitomap import MitomapDatabase
from vartriage.models.variant import Variant

logger = logging.getLogger(__name__)


class MitochondrialPipeline:
    """Sub-pipeline for mtDNA variant classification.

    Constructs the three databases (gene map, MITOMAP, HelixMTdb) at
    init time from config paths or package defaults, then classifies
    variants via the MitochondrialClassifier.

    Parameters
    ----------
    config
        MitoConfig with threshold settings and optional custom data paths.
    """

    def __init__(self, config: MitoConfig) -> None:
        self._config = config

        self._gene_map = MtGeneMap(data_path=config.gene_map_path)
        self._mitomap_db = MitomapDatabase(data_path=config.mitomap_path)
        self._helix_db = HelixMTdbDatabase(data_path=config.helixmtdb_path)

        self._classifier = MitochondrialClassifier(
            gene_map=self._gene_map,
            mitomap_db=self._mitomap_db,
            helix_db=self._helix_db,
        )

        logger.info(
            "MitochondrialPipeline initialized: "
            "gene_map=%d intervals, mitomap=%d entries, helixmtdb=%d entries",
            len(self._gene_map.entries),
            self._mitomap_db.size,
            self._helix_db.size,
        )

    def run(self, variants: Iterator[Variant]) -> list[MitoClassifiedVariant]:
        """Classify a stream of mitochondrial variants.

        Filters sub-threshold variants based on min_heteroplasmy config,
        then classifies the remainder.

        Parameters
        ----------
        variants
            Iterator of Variant records (all expected to be chrM/MT).

        Returns
        -------
        list[MitoClassifiedVariant]
            Classified mitochondrial variants that pass the heteroplasmy
            threshold filter.
        """
        results: list[MitoClassifiedVariant] = []
        filtered_count = 0

        for variant in variants:
            classified = self._classifier.classify(variant)

            # Filter sub-threshold variants
            if (
                classified.heteroplasmy is not None
                and classified.heteroplasmy.percentage < self._config.min_heteroplasmy
            ):
                filtered_count += 1
                continue

            results.append(classified)

        if filtered_count > 0:
            logger.info(
                "Filtered %d mtDNA variants below %.1f%% heteroplasmy threshold",
                filtered_count,
                self._config.min_heteroplasmy,
            )

        logger.info("MitochondrialPipeline classified %d variants", len(results))
        return results
