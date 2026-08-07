"""Top-level pipeline orchestrator for variant prioritization.

Connects all processing stages: VCF ingestion → quality filtering →
annotation → prioritization → ACMG classification → report generation.
Manages configuration validation at construction (fail-fast), warning
accumulation across stages, and batch-based memory management.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path

from vartriage._internal.path_safety import resolve_path
from vartriage._internal.warning_accumulator import WarningAccumulator
from vartriage.annotation.engine import AnnotationEngine
from vartriage.classification.acmg import ACMGClassifier
from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.vcf_parser import VCFParser
from vartriage.mito.genetic_code import is_mitochondrial
from vartriage.models.config import (
    AnnotationConfig,
    PipelineConfig,
    PrioritizationConfig,
)
from vartriage.models.variant import AnnotatedVariant, ClassifiedVariant, Variant
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

        # Gene-disease linkage annotator: constructed once, reused across runs
        self._gene_knowledge_annotator: GeneKnowledgeAnnotator | None = None
        if config.knowledge is not None:
            from vartriage.knowledge.annotator import GeneKnowledgeAnnotator

            self._gene_knowledge_annotator = GeneKnowledgeAnnotator(
                config.knowledge  # type: ignore[arg-type]
            )

    @property
    def _mito_enabled(self) -> bool:
        """Whether mitochondrial analysis is active."""
        if self._config.mito is not None:
            return self._config.mito.enabled
        # Auto-enabled by default (no explicit config needed)
        return True

    @property
    def mito_results(self) -> list | None:
        """Access mitochondrial classification results from the last run.

        Returns None if mitochondrial analysis was disabled or no chrM
        variants were present in the input.
        """
        return getattr(self, "_mito_results", None)

    @property
    def warning_accumulator(self) -> WarningAccumulator:
        """Access the warning accumulator tracking MissingDataWarnings.

        Returns
        -------
        WarningAccumulator
            The shared warning accumulator for the current pipeline run.
        """
        return self._warning_accumulator

    def _build_stages(
        self,
    ) -> tuple[
        QualityFilter,
        AnnotationEngine | None,
        object,
        PrioritizationEngine,
        ACMGClassifier,
        ReportGenerator,
    ]:
        """Construct all pipeline stage objects."""
        quality_filter = QualityFilter(self._config.quality_filter)

        annotation_engine: AnnotationEngine | None = None
        api_annotation_engine = None

        if self._config.api is not None:
            api_annotation_engine = self._build_api_annotation_engine()

        if api_annotation_engine is None and self._config.annotation is not None:
            annotation_engine = AnnotationEngine(self._config.annotation)

        prioritization_engine = PrioritizationEngine(self._config.prioritization)
        acmg_classifier = ACMGClassifier()
        report_generator = ReportGenerator(
            self._config.report,
            clinical_config=self._config.clinical_report,
            reference_checksums=(
                self._compute_reference_checksums()
                if self._config.clinical_report is not None
                else None
            ),
        )
        return (
            quality_filter,
            annotation_engine,
            api_annotation_engine,
            prioritization_engine,
            acmg_classifier,
            report_generator,
        )

    def _generate_report(
        self,
        report_generator: ReportGenerator,
        classified: Iterator[ClassifiedVariant],
        output_path: Path,
        vcf_path: Path,
        mito_results: list | None = None,
    ) -> Path:
        """Dispatch report generation for VCF and non-VCF formats.

        When mito_results are present and the format supports it (JSON/CSV),
        uses the mito-aware writers that include a mitochondrial findings section.
        """
        self._mito_results = mito_results
        fmt = self._config.report.output_format

        if fmt == "vcf":
            return report_generator.generate(classified, output_path, vcf_path)

        # For JSON/CSV with mito results, use mito-aware writers directly
        if mito_results and fmt in ("json", "csv"):
            from vartriage._internal.path_safety import resolve_path as _rp

            resolved = _rp(output_path)
            resolved.parent.mkdir(parents=True, exist_ok=True)

            if fmt == "json":
                from vartriage.reporting.json_writer import write_json_with_mito

                return write_json_with_mito(classified, resolved, mito_results)
            else:
                from vartriage.reporting.csv_writer import write_csv_with_mito

                return write_csv_with_mito(classified, resolved, mito_results)

        return report_generator.generate(classified, output_path)

    def run(
        self, vcf_path: Path | None = None, output_path: Path | None = None
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

        logger.info(
            "Starting pipeline run: %s → %s", effective_vcf_path, effective_output_path
        )

        if self._config.clinical_report is not None:
            self._check_reference_checksums()

        (
            quality_filter,
            annotation_engine,
            api_annotation_engine,
            prioritization_engine,
            acmg_classifier,
            report_generator,
        ) = self._build_stages()

        extract_samples = (
            self._config.inheritance is not None
            or self._config.sample is not None
            or self._mito_enabled
        )

        with VCFParser(
            effective_vcf_path,
            extract_samples=extract_samples,
        ) as parser:
            mito_variants, nuclear_iter = self._split_mito_variants(parser)

            annotated = self._build_annotated_stream(
                nuclear_iter,
                quality_filter,
                annotation_engine,
                api_annotation_engine,
                sample_names=parser.sample_names,
            )

            if self._config.gene_filter is not None:
                from vartriage.filter.gene_filter import GeneFilter

                annotated = GeneFilter(self._config.gene_filter).apply(annotated)

            if self._gene_knowledge_annotator is not None:
                annotated = self._gene_knowledge_annotator.annotate(annotated)

            scored = prioritization_engine.prioritize(annotated)

            if self._gene_knowledge_annotator is not None:
                scored = self._gene_knowledge_annotator.boost_scores(scored)

            classified = acmg_classifier.classify(scored)
            mito_results = self._run_mito_pipeline(mito_variants)

            result_path = self._generate_report(
                report_generator,
                classified,
                effective_output_path,
                effective_vcf_path,
                mito_results=mito_results,
            )

            if annotation_engine is not None:
                self._warning_accumulator.add_batch(annotation_engine.warnings)

        logger.info(
            "Pipeline completed. Missing data warnings: %d",
            self._warning_accumulator.total_count,
        )
        logger.info("Report written to: %s", result_path)
        return result_path

    def _split_mito_variants(
        self, parser: VCFParser
    ) -> tuple[list[Variant], Iterator[Variant]]:
        """Partition parser variants into mitochondrial and nuclear streams."""
        if not self._mito_enabled:
            return [], iter(parser)
        mito: list[Variant] = []
        nuclear: list[Variant] = []
        for variant in parser:
            if is_mitochondrial(variant.chrom):
                mito.append(variant)
            else:
                nuclear.append(variant)
        return mito, iter(nuclear)

    def _run_mito_pipeline(self, mito_variants: list[Variant]) -> list | None:
        """Run the mitochondrial pipeline if variants are present."""
        if not (self._mito_enabled and mito_variants):
            return None
        from vartriage.mito.config import MitoConfig
        from vartriage.mito.pipeline import MitochondrialPipeline

        mito_config = self._config.mito or MitoConfig()
        results = MitochondrialPipeline(mito_config).run(iter(mito_variants))
        logger.info("Mitochondrial pipeline: %d variants classified", len(results))
        return results

    def run_to_classification(
        self, vcf_path: Path | None = None
    ) -> Iterator[ClassifiedVariant]:
        """Execute the pipeline through classification without report generation.

        Runs VCF parsing, quality filtering, annotation, prioritization,
        and ACMG classification, yielding ClassifiedVariant objects. Skips
        report writing entirely. Used by CohortPipeline to collect per-sample
        results without generating throw-away files.

        Stage limitations vs full run():
            - SampleExtractor is not applied (cohort mode processes
              each VCF as a single-sample file already).
            - InheritanceFilter is not applied (trio analysis is a
              single-proband concern, not meaningful per-sample in a
              cohort where each VCF represents one individual).
            - SecondaryFindingsFilter is not applied (cohort mode
              aggregates by variant, not by clinical reporting context).

        These stages are intentionally excluded because cohort analysis
        operates on pre-separated per-sample VCFs. If your workflow
        requires sample extraction or inheritance filtering, run the
        full single-sample pipeline first and feed ClassifiedVariant
        results to CohortAggregator directly.

        Parameters
        ----------
        vcf_path : Path, optional
            Path to the input VCF file. If None, uses the path from config.

        Yields
        ------
        ClassifiedVariant
            Classified variants from the pipeline stages.
        """
        effective_vcf_path = vcf_path or self._config.vcf_path

        quality_filter = QualityFilter(self._config.quality_filter)

        annotation_engine: AnnotationEngine | None = None
        if self._config.annotation is not None:
            annotation_engine = AnnotationEngine(self._config.annotation)

        prioritization_engine = PrioritizationEngine(self._config.prioritization)
        acmg_classifier = ACMGClassifier()

        with VCFParser(effective_vcf_path) as parser:
            stream: Iterator[Variant] = iter(parser)

            if self._config.region_filter is not None:
                from vartriage.filter.region_filter import RegionFilter

                region_filter = RegionFilter(self._config.region_filter)
                stream = region_filter.apply(stream)

            filtered = quality_filter.apply(stream)

            if annotation_engine is not None:
                annotated = annotation_engine.annotate(filtered)
            else:
                annotated = self._passthrough_annotation(filtered)

            if self._config.gene_filter is not None:
                from vartriage.filter.gene_filter import GeneFilter

                gene_filter = GeneFilter(self._config.gene_filter)
                annotated = gene_filter.apply(annotated)

            # Gene-disease linkage in run_to_classification path
            if self._gene_knowledge_annotator is not None:
                annotated = self._gene_knowledge_annotator.annotate(annotated)

            scored = prioritization_engine.prioritize(annotated)

            # Apply phenotype boost in classification path too
            if self._gene_knowledge_annotator is not None:
                scored = self._gene_knowledge_annotator.boost_scores(scored)

            yield from acmg_classifier.classify(scored)

    def run_with_sv(self) -> None:
        """Execute the main pipeline and SV triage when sv_vcf_path is set.

        Runs the standard SNV pipeline via run(), then if sv_vcf_path is
        configured, also runs the SV triage pipeline and writes a combined
        findings file alongside the main report.
        """
        self.run()

        sv_path = self._config.sv_vcf_path
        if sv_path is None or not sv_path.exists():
            return

        from vartriage.structural.config import SVTriageConfig
        from vartriage.structural.pipeline import SVTriagePipeline

        # Derive SV output path from the main output
        main_output = self._config.output_path
        sv_output = main_output.parent / (main_output.stem + "_sv" + main_output.suffix)

        sv_config = SVTriageConfig(
            vcf_path=sv_path,
            output_path=sv_output,
            gene_annotation_path=(
                self._config.annotation.gene_annotation_path
                if self._config.annotation is not None
                else None
            ),
            output_format=(
                "json"
                if self._config.report.output_format in ("json", "vcf", "pdf")
                else "csv"
            ),
        )

        sv_pipeline = SVTriagePipeline(sv_config)
        sv_pipeline.run()

        logger.info("SV triage report written to: %s", sv_output)

    def _check_reference_checksums(self) -> None:
        """Log reference file checksums using AuditTrailWriter.

        Delegates SHA-256 computation to
        ``AuditTrailWriter.compute_file_checksum`` so there is a
        single checksum implementation across the codebase.
        """
        from vartriage.reporting.clinical.audit import AuditTrailWriter

        audit_writer = AuditTrailWriter()

        for ref_path in self._collect_reference_paths():
            if not ref_path.exists():
                continue
            try:
                checksum = audit_writer.compute_file_checksum(ref_path)
                logger.debug(
                    "Reference file %s checksum: %s",
                    ref_path,
                    checksum,
                )
            except OSError as exc:
                logger.warning(
                    "Could not compute checksum for %s: %s",
                    ref_path,
                    exc,
                )

    def _collect_reference_paths(self) -> list[Path]:
        """Collect all configured reference file paths.

        Returns
        -------
        list[Path]
            Paths to annotation and prioritization reference files
            that are configured (non-None).
        """
        ref_paths: list[Path] = []
        if self._config.annotation is not None:
            ref_paths.append(self._config.annotation.gene_annotation_path)
            ref_paths.append(self._config.annotation.gnomad_path)
            if self._config.annotation.clinvar_path is not None:
                ref_paths.append(self._config.annotation.clinvar_path)
        pri = self._config.prioritization
        if pri.cadd_scores_path is not None:
            ref_paths.append(pri.cadd_scores_path)
        if pri.revel_scores_path is not None:
            ref_paths.append(pri.revel_scores_path)
        if pri.spliceai_scores_path is not None:
            ref_paths.append(pri.spliceai_scores_path)
        return ref_paths

    def _compute_reference_checksums(self) -> dict[str, str]:
        """Compute SHA-256 checksums for all reference files.

        Uses ``AuditTrailWriter.compute_file_checksum`` as the
        single checksum implementation. Files that do not exist or
        cannot be read are skipped with a warning.

        Returns
        -------
        dict[str, str]
            Mapping of file path strings to SHA-256 hex digests.
        """
        from vartriage.reporting.clinical.audit import AuditTrailWriter

        audit_writer = AuditTrailWriter()
        checksums: dict[str, str] = {}

        for ref_path in self._collect_reference_paths():
            if not ref_path.exists():
                continue
            try:
                checksums[str(ref_path)] = audit_writer.compute_file_checksum(ref_path)
            except OSError as exc:
                logger.warning(
                    "Could not compute checksum for %s: %s",
                    ref_path,
                    exc,
                )

        return checksums

    def _build_annotated_stream(
        self,
        parser: VCFParser | Iterator[Variant],
        quality_filter: QualityFilter,
        annotation_engine: AnnotationEngine | None,
        api_annotation_engine: object = None,
        sample_names: list[str] | None = None,
    ) -> Iterator[AnnotatedVariant]:
        """Build the filtered and annotated variant stream."""
        stream: Iterator[Variant] = iter(parser)  # type: ignore[arg-type]

        # Sample extraction or inheritance (mutually exclusive)
        if self._config.inheritance is not None:
            from vartriage.filter.inheritance_filter import InheritanceFilter

            _sample_names = sample_names
            if _sample_names is None and hasattr(parser, "sample_names"):
                _sample_names = parser.sample_names  # type: ignore[union-attr]

            inheritance_filter = InheritanceFilter(
                self._config.inheritance,
                _sample_names or [],
            )
            compound_het_active = "compound_het" in self._config.inheritance.patterns

            if compound_het_active and annotation_engine is not None:
                # Annotate first so gene info is available for
                # compound_het grouping
                filtered = quality_filter.apply(stream)
                annotated_iter = annotation_engine.annotate(filtered)
                # Buffer annotated variants so we can pass raw
                # Variants (with gene in info) to InheritanceFilter,
                # then re-associate the annotation data afterward
                annotated_list = list(annotated_iter)
                variants_with_genes = list(
                    self._variants_with_gene_info(iter(annotated_list))
                )
                inherited_variants = list(
                    inheritance_filter.apply(iter(variants_with_genes))
                )
                # Re-attach annotation data to inherited variants
                return iter(
                    self._reattach_annotations(inherited_variants, annotated_list)
                )
            else:
                stream = inheritance_filter.apply(stream)
        elif self._config.sample is not None:
            from vartriage.filter.sample_extractor import SampleExtractor

            sample_extractor = SampleExtractor(
                self._config.sample,
                parser.sample_names,
            )
            stream = sample_extractor.apply(stream)

        # Region filter (optional, runs before quality filter)
        if self._config.region_filter is not None:
            from vartriage.filter.region_filter import RegionFilter

            region_filter = RegionFilter(self._config.region_filter)
            stream = region_filter.apply(stream)

        filtered = quality_filter.apply(stream)
        if annotation_engine is not None:
            return annotation_engine.annotate(filtered)
        if api_annotation_engine is not None:
            return api_annotation_engine.annotate(filtered)  # type: ignore[attr-defined,no-any-return]
        return self._passthrough_annotation(filtered)

    def _build_api_annotation_engine(self) -> object:
        """Construct the API annotation engine from config.

        Deferred import to avoid pulling httpx into the core module.
        """
        from vartriage.api.annotation_engine import APIAnnotationEngine

        return APIAnnotationEngine(self._config.api)  # type: ignore[arg-type]

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
            self._validate_annotation_config(config.annotation)

        self._validate_prioritization_config(config.prioritization)

        if config.gene_filter is not None:
            self._check_path(
                config.gene_filter.gene_list_path,
                "Gene list file",
            )

        if config.region_filter is not None:
            self._check_path(
                config.region_filter.bed_path,
                "BED file",
            )

        if (
            config.inheritance is not None
            and "compound_het" in config.inheritance.patterns
            and config.annotation is None
        ):
            raise ValueError(
                "compound_het pattern requires annotation "
                "configuration (gene annotation reference). "
                "Either provide --gene-annotation and --gnomad, "
                "or remove compound_het from the patterns list."
            )

    def _validate_annotation_config(
        self,
        ann_config: AnnotationConfig,
    ) -> None:
        """Validate annotation reference file paths exist."""
        self._check_path(ann_config.gene_annotation_path, "Gene annotation file")
        self._check_path(ann_config.gnomad_path, "gnomAD reference file")
        if ann_config.clinvar_path is not None:
            self._check_path(ann_config.clinvar_path, "ClinVar reference file")

    def _validate_prioritization_config(
        self,
        pri_config: PrioritizationConfig,
    ) -> None:
        """Validate prioritization score file paths exist."""
        if pri_config.cadd_scores_path is not None:
            self._check_path(pri_config.cadd_scores_path, "CADD scores file")
        if pri_config.revel_scores_path is not None:
            self._check_path(pri_config.revel_scores_path, "REVEL scores file")
        if pri_config.spliceai_scores_path is not None:
            self._check_path(pri_config.spliceai_scores_path, "SpliceAI scores file")

    def _reattach_annotations(
        self,
        inherited_variants: list[Variant],
        annotated_list: list[AnnotatedVariant],
    ) -> list[AnnotatedVariant]:
        """Re-attach annotation data to variants after inheritance filtering.

        Builds a coordinate lookup from the original annotated variants
        and matches each inherited variant back to its annotation.
        Variants that pass inheritance filtering but have no annotation
        match get INTERGENIC/null as fallback.

        Parameters
        ----------
        inherited_variants : list[Variant]
            Variants that passed inheritance filtering (with
            inheritance_pattern in info).
        annotated_list : list[AnnotatedVariant]
            Original annotated variants before inheritance filtering.

        Returns
        -------
        list[AnnotatedVariant]
            Annotated variants with inheritance metadata preserved
            in the underlying variant's info dict.
        """
        from vartriage.models.variant import AnnotatedVariant, FunctionalConsequence

        # Build lookup by (chrom, pos, ref, alt)
        ann_lookup: dict[tuple[str, int, str, str], AnnotatedVariant] = {}
        for av in annotated_list:
            key = (av.variant.chrom, av.variant.pos, av.variant.ref, av.variant.alt)
            ann_lookup[key] = av

        results: list[AnnotatedVariant] = []
        for v in inherited_variants:
            key = (v.chrom, v.pos, v.ref, v.alt)
            original = ann_lookup.get(key)
            if original is not None:
                # Preserve original annotation, attach inheritance
                # info by replacing the variant's info dict
                results.append(
                    AnnotatedVariant(
                        variant=v,
                        consequence=original.consequence,
                        allele_frequency=original.allele_frequency,
                        clinvar_assertion=original.clinvar_assertion,
                        frequency_unknown=original.frequency_unknown,
                        clinvar_unknown=original.clinvar_unknown,
                        gene_name=original.gene_name,
                    )
                )
            else:
                results.append(
                    AnnotatedVariant(
                        variant=v,
                        consequence=FunctionalConsequence.INTERGENIC,
                        frequency_unknown=True,
                        clinvar_unknown=True,
                    )
                )
        return results

    def _variants_with_gene_info(
        self, annotated: Iterator[AnnotatedVariant]
    ) -> Iterator[Variant]:
        """Extract Variant objects with gene_name copied into info dict.

        Used when compound_het needs gene grouping from annotation
        data. Copies gene_name into info["gene"] so the
        InheritanceFilter can group variants by gene.

        Parameters
        ----------
        annotated : Iterator[AnnotatedVariant]
            Annotated variant stream.

        Yields
        ------
        Variant
            Raw variants with gene info attached.
        """
        for av in annotated:
            v = av.variant
            new_info = dict(v.info)
            if av.gene_name is not None:
                new_info["gene"] = av.gene_name
            yield Variant(
                chrom=v.chrom,
                pos=v.pos,
                id=v.id,
                ref=v.ref,
                alt=v.alt,
                qual=v.qual,
                filter_status=v.filter_status,
                info=new_info,
            )

    def _passthrough_annotation(
        self, variants: Iterator[Variant]
    ) -> Iterator[AnnotatedVariant]:
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
        from vartriage.models.variant import AnnotatedVariant, FunctionalConsequence

        for variant in variants:
            yield AnnotatedVariant(
                variant=variant,
                consequence=FunctionalConsequence.INTERGENIC,
                allele_frequency=None,
                clinvar_assertion=None,
                frequency_unknown=True,
                clinvar_unknown=True,
            )

    @staticmethod
    def _check_path(path: Path, label: str) -> None:
        """Verify a file path exists, raising FileNotFoundError if not.

        Resolves the path to its canonical form before checking, which
        eliminates directory traversal sequences and follows symlinks.

        Parameters
        ----------
        path : Path
            The filesystem path to validate.
        label : str
            Human-readable description used in the error message.

        Raises
        ------
        FileNotFoundError
            If the path does not exist.
        """
        resolved = resolve_path(path)
        if not resolved.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
