"""Cohort pipeline orchestrator.

Runs the standard single-sample Pipeline for each VCF in the cohort,
collects classified variants, then aggregates, computes statistics,
and generates cohort-level reports. Supports both sequential and
parallel (thread pool) sample processing.
"""

from __future__ import annotations

import logging
import tempfile
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from vartriage.cohort.aggregator import CohortAggregator
from vartriage.cohort.report import CohortReportGenerator
from vartriage.cohort.statistics import CohortStatistics
from vartriage.models.cohort import (
    CohortConfig,
    CohortSummary,
    CohortVariant,
    GeneBurden,
)
from vartriage.models.config import (
    AnnotationConfig,
    PipelineConfig,
    PrioritizationConfig,
    ReportConfig,
)
from vartriage.models.variant import ClassifiedVariant

logger = logging.getLogger(__name__)


class CohortPipeline:
    """Orchestrate multi-sample cohort analysis.

    Processes each sample VCF through the standard vartriage pipeline,
    collects classified variants, then merges them via CohortAggregator
    for cross-sample analysis.

    Parameters
    ----------
    cohort_config : CohortConfig
        Cohort-level settings (sample list, thresholds, output).
    pipeline_config : PipelineConfig | None
        Base pipeline configuration applied to each sample. When None,
        a minimal config is constructed per-sample using only the
        cohort_config's sample VCF paths with default quality/prioritization
        settings.
    annotation_config : AnnotationConfig | None
        Shared annotation config for all samples. Overrides the
        pipeline_config's annotation setting when provided.
    prioritization_config : PrioritizationConfig | None
        Shared prioritization config. Overrides pipeline_config when provided.
    """

    def __init__(
        self,
        cohort_config: CohortConfig,
        pipeline_config: Optional[PipelineConfig] = None,
        annotation_config: Optional[AnnotationConfig] = None,
        prioritization_config: Optional[PrioritizationConfig] = None,
    ) -> None:
        self._cohort_config = cohort_config
        self._base_pipeline_config = pipeline_config
        self._annotation_config = annotation_config
        self._prioritization_config = prioritization_config
        self._aggregator = CohortAggregator(cohort_config)
        self._tmp_dir: Optional[tempfile.TemporaryDirectory[str]] = None

        # Results populated after run()
        self._variants: list[CohortVariant] = []
        self._gene_burdens: list[GeneBurden] = []
        self._summary: Optional[CohortSummary] = None
        self._samples_processed: list[str] = []

    @property
    def variants(self) -> list[CohortVariant]:
        """Aggregated cohort variants (populated after run())."""
        return self._variants

    @property
    def gene_burdens(self) -> list[GeneBurden]:
        """Per-gene burden records (populated after run())."""
        return self._gene_burdens

    @property
    def summary(self) -> Optional[CohortSummary]:
        """Cohort summary statistics (populated after run())."""
        return self._summary

    def run(self) -> list[Path]:
        """Execute the full cohort analysis pipeline.

        Sequence:
        1. Process each sample VCF through the standard pipeline
        2. Aggregate variants across samples
        3. Compute cohort statistics
        4. Generate reports

        Returns
        -------
        list[Path]
            Paths to generated report files.

        Raises
        ------
        FileNotFoundError
            If any sample VCF file does not exist.
        """
        logger.info(
            "Starting cohort analysis '%s' with %d samples",
            self._cohort_config.cohort_name,
            self._cohort_config.sample_count,
        )

        # Validate all VCF files exist before processing
        for vcf_path in self._cohort_config.sample_vcfs:
            if not vcf_path.exists():
                raise FileNotFoundError(
                    f"Sample VCF not found: {vcf_path}"
                )

        # Process samples
        if self._cohort_config.parallel:
            self._process_parallel()
        else:
            self._process_sequential()

        # Aggregate
        logger.info("Aggregating variants across %d samples", len(self._samples_processed))
        self._variants = self._aggregator.aggregate()

        # Statistics
        stats = CohortStatistics(self._cohort_config, self._variants)
        self._gene_burdens = stats.compute_gene_burden()
        self._summary = stats.compute_summary(self._samples_processed)

        # Report generation
        reporter = CohortReportGenerator(self._cohort_config)
        report_paths = reporter.generate(
            self._variants, self._gene_burdens, self._summary
        )

        logger.info(
            "Cohort analysis complete: %d variants, %d shared, %d genes",
            self._summary.total_variants,
            self._summary.shared_variants,
            self._summary.genes_affected,
        )

        # Clean up temporary directory used for placeholder output paths
        if self._tmp_dir is not None:
            self._tmp_dir.cleanup()
            self._tmp_dir = None

        return report_paths

    def _process_sequential(self) -> None:
        """Process each sample VCF sequentially."""
        for vcf_path in self._cohort_config.sample_vcfs:
            sample_id = self._cohort_config.label_for(vcf_path)
            classified = self._run_single_sample(vcf_path, sample_id)
            self._aggregator.add_sample(sample_id, vcf_path, classified)
            self._samples_processed.append(sample_id)

    def _process_parallel(self) -> None:
        """Process sample VCFs concurrently using a thread pool.

        Thread-based parallelism works here because the per-sample
        pipeline is I/O-bound (VCF parsing, reference file reads).
        """
        max_workers = self._cohort_config.max_workers
        futures_map: dict[Future[list[ClassifiedVariant]], tuple[str, Path]] = {}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for vcf_path in self._cohort_config.sample_vcfs:
                sample_id = self._cohort_config.label_for(vcf_path)
                future = executor.submit(
                    self._run_single_sample, vcf_path, sample_id
                )
                futures_map[future] = (sample_id, vcf_path)

            for future in as_completed(futures_map):
                sample_id, vcf_path = futures_map[future]
                try:
                    classified: list[ClassifiedVariant] = future.result()
                    self._aggregator.add_sample(sample_id, vcf_path, classified)
                    self._samples_processed.append(sample_id)
                except Exception:
                    logger.exception(
                        "Failed to process sample '%s' (%s)",
                        sample_id,
                        vcf_path,
                    )
                    raise

    def _run_single_sample(
        self, vcf_path: Path, sample_id: str
    ) -> list[ClassifiedVariant]:
        """Run the standard pipeline on a single VCF and collect results.

        Delegates to Pipeline.run_to_classification() which executes
        all stages up through ACMG classification without writing a
        report file.

        Parameters
        ----------
        vcf_path : Path
            Path to the sample's VCF file.
        sample_id : str
            Sample identifier for logging.

        Returns
        -------
        list[ClassifiedVariant]
            All classified variants from this sample.
        """
        from vartriage.pipeline import Pipeline

        logger.info("Processing sample '%s': %s", sample_id, vcf_path)

        config = self._build_sample_config(vcf_path)
        pipeline = Pipeline(config)
        classified = list(pipeline.run_to_classification(vcf_path))

        logger.info(
            "Sample '%s' yielded %d classified variants",
            sample_id,
            len(classified),
        )
        return classified

    def _build_sample_config(self, vcf_path: Path) -> PipelineConfig:
        """Build a PipelineConfig for a single sample.

        If a base pipeline_config was provided, clones it with the
        sample's VCF path. Otherwise builds a minimal config.
        """
        # Placeholder output path (run_to_classification never writes to it)
        if self._tmp_dir is None:
            self._tmp_dir = tempfile.TemporaryDirectory(prefix="vartriage_cohort_")
        tmp_output = Path(self._tmp_dir.name) / f"{vcf_path.stem}_output.json"

        if self._base_pipeline_config is not None:
            # Clone the base config with this sample's VCF path
            base = self._base_pipeline_config
            return PipelineConfig(
                vcf_path=vcf_path,
                output_path=tmp_output,
                quality_filter=base.quality_filter,
                annotation=self._annotation_config or base.annotation,
                prioritization=self._prioritization_config or base.prioritization,
                report=ReportConfig(output_format="json"),
                missing_data=base.missing_data,
                gene_filter=base.gene_filter,
                region_filter=base.region_filter,
            )

        # Minimal config when no base is provided
        return PipelineConfig(
            vcf_path=vcf_path,
            output_path=tmp_output,
            annotation=self._annotation_config,
            prioritization=(
                self._prioritization_config or PrioritizationConfig()
            ),
            report=ReportConfig(output_format="json"),
        )
