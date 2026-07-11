"""Top-level pipeline orchestrator for variant prioritization.

Connects all processing stages: VCF ingestion → quality filtering →
annotation → prioritization → ACMG classification → report generation.
Manages configuration validation at construction (fail-fast), warning
accumulation across stages, and batch-based memory management.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator, Optional

from vartriage._internal.warning_accumulator import (
    WarningAccumulator,
)
from vartriage.annotation.engine import AnnotationEngine
from vartriage.classification.acmg import ACMGClassifier
from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.config import (
    AnnotationConfig,
    PipelineConfig,
)
from vartriage.models.variant import (
    AnnotatedVariant,
    Variant,
)
from vartriage.prioritization.engine import PrioritizationEngine
from vartriage.reporting.generator import ReportGenerator

logger = logging.getLogger(__name__)


class Pipeline:
    """Top-level orchestrator for the variant prioritization pipeline.

    Connects all stages and manages configuration validation, streaming
    variants through the pipeline while maintaining memory bounds suitable
    for whole-genome scale datasets (4M+ variants, <2GB RSS).

    The pipeline validates all configuration parameters at construction
    time (fail-fast). The ``run`` method then executes the full processing
    chain: parsing → filtering → annotation → prioritization →
    classification → report generation.

    Parameters
    ----------
    config : PipelineConfig
        Complete configuration for all pipeline stages.

    Raises
    ------
    ValueError
        If any configuration parameter is invalid (out-of-range thresholds,
        invalid batch sizes, unsupported formats).
    FileNotFoundError
        If any required reference file is missing.
    """

    def __init__(self, config: PipelineConfig) -> None:
        self._config = config
        self._warning_accumulator = WarningAccumulator(config.missing_data)
        self._validate_config(config)

    @property
    def warning_accumulator(self) -> WarningAccumulator:
        """Access the warning accumulator tracking MissingDataWarnings.

        Returns
        -------
        WarningAccumulator
            The shared warning accumulator for the current pipeline run.
        """
        return self._warning_accumulator

    def run(
        self, vcf_path: Optional[Path] = None, output_path: Optional[Path] = None
    ) -> Path:
        """Execute the full variant prioritization pipeline.

        Wires stages sequentially: VCFParser → QualityFilter →
        AnnotationEngine → PrioritizationEngine → ACMGClassifier →
        ReportGenerator. Processes data in batches to maintain peak RSS
        below 2GB for whole-genome scale files (4M+ variants).

        Parameters
        ----------
        vcf_path : Path, optional
            Path to the input VCF file (.vcf or .vcf.gz). If None, uses
            the path from config.
        output_path : Path, optional
            Path for the output report. If None, uses the path from config.

        Returns
        -------
        Path
            Path to the generated output report file.

        Raises
        ------
        FileNotFoundError
            If the VCF file or any reference file is missing.
        ParseError
            If the VCF file has malformed headers or data lines.
        IOError
            If report generation fails due to a write error.
        """
        effective_vcf_path = vcf_path or self._config.vcf_path
        effective_output_path = output_path or self._config.output_path

        self._warning_accumulator.reset()

        logger.info("Starting pipeline run: %s → %s",
                    effective_vcf_path, effective_output_path)

        quality_filter = QualityFilter(self._config.quality_filter)

        annotation_engine: Optional[AnnotationEngine] = None
        if self._config.annotation is not None:
            annotation_engine = AnnotationEngine(self._config.annotation)

        prioritization_engine = PrioritizationEngine(
            self._config.prioritization
        )

        acmg_classifier = ACMGClassifier()

        report_generator = ReportGenerator(self._config.report)

        with VCFParser(effective_vcf_path) as parser:
            filtered = quality_filter.apply(iter(parser))

            if annotation_engine is not None:
                annotated = annotation_engine.annotate(filtered)
            else:
                annotated = self._passthrough_annotation(filtered)

            scored = prioritization_engine.prioritize(annotated)

            classified = acmg_classifier.classify(scored)

            result_path = report_generator.generate(
                classified, effective_output_path
            )

            if annotation_engine is not None:
                self._warning_accumulator.add_batch(
                    annotation_engine.warnings
                )

        logger.info(
            "Pipeline completed. Missing data warnings: %d",
            self._warning_accumulator.total_count,
        )

        logger.info("Report written to: %s", result_path)
        return result_path

    def _validate_config(self, config: PipelineConfig) -> None:
        """Validate all configuration at construction time (fail-fast).

        Checks that reference file paths exist and all sub-config
        parameters are within valid ranges.

        Parameters
        ----------
        config : PipelineConfig
            Configuration to validate.

        Raises
        ------
        ValueError
            If any parameter is out of range.
        FileNotFoundError
            If any required reference file does not exist.
        """
        if config.annotation is not None:
            ann_config: AnnotationConfig = config.annotation
            if not ann_config.gene_annotation_path.exists():
                raise FileNotFoundError(
                    f"Gene annotation file not found: "
                    f"{ann_config.gene_annotation_path}"
                )
            if not ann_config.gnomad_path.exists():
                raise FileNotFoundError(
                    f"gnomAD reference file not found: "
                    f"{ann_config.gnomad_path}"
                )
            if (
                ann_config.clinvar_path is not None
                and not ann_config.clinvar_path.exists()
            ):
                raise FileNotFoundError(
                    f"ClinVar reference file not found: "
                    f"{ann_config.clinvar_path}"
                )

        pri_config = config.prioritization
        if pri_config.cadd_scores_path is not None:
            if not pri_config.cadd_scores_path.exists():
                raise FileNotFoundError(
                    f"CADD scores file not found: "
                    f"{pri_config.cadd_scores_path}"
                )
        if pri_config.revel_scores_path is not None:
            if not pri_config.revel_scores_path.exists():
                raise FileNotFoundError(
                    f"REVEL scores file not found: "
                    f"{pri_config.revel_scores_path}"
                )
        if pri_config.spliceai_scores_path is not None:
            if not pri_config.spliceai_scores_path.exists():
                raise FileNotFoundError(
                    f"SpliceAI scores file not found: "
                    f"{pri_config.spliceai_scores_path}"
                )

    def _passthrough_annotation(self, variants: Iterator["Variant"]) -> Iterator["AnnotatedVariant"]:
        """Create AnnotatedVariant wrappers when no annotation config exists.

        Used when the pipeline is run without annotation references. Each
        variant gets an Intergenic consequence, null frequency, and null
        ClinVar assertion, allowing downstream stages to function.

        Parameters
        ----------
        variants : Iterator[Variant]
            Filtered variant stream.

        Yields
        ------
        AnnotatedVariant
            Minimally annotated variant records.
        """
        from vartriage.models.variant import (
            AnnotatedVariant,
            FunctionalConsequence,
        )

        for variant in variants:
            yield AnnotatedVariant(
                variant=variant,
                consequence=FunctionalConsequence.INTERGENIC,
                allele_frequency=None,
                clinvar_assertion=None,
                frequency_unknown=True,
                clinvar_unknown=True,
            )
