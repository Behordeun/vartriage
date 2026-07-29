"""Structural variant triage pipeline orchestrator.

Connects all SV processing stages: VCF parsing → quality/size filtering →
gene annotation → dosage scoring → ClinGen classification → report output.
Validates configuration at construction (fail-fast) and processes SVs in
a streaming fashion.
"""

from __future__ import annotations

import json
import csv as csv_module
import logging
from pathlib import Path
from typing import Any, Iterator, Optional

from vartriage.structural.annotator import SVAnnotator
from vartriage.structural.classifier import SVClassifier
from vartriage.structural.config import SVTriageConfig
from vartriage.structural.models import (
    ClassifiedSV,
    SVClassification,
)
from vartriage.structural.parser import SVParser
from vartriage.structural.scoring import SVScorer

logger = logging.getLogger(__name__)


class SVTriagePipeline:
    """Orchestrate structural variant triage from VCF to classified report.

    Wires stages sequentially: SVParser → SVAnnotator → SVScorer →
    SVClassifier → report output. Configuration validation happens at
    construction so errors surface before processing starts.

    Parameters
    ----------
    config : SVTriageConfig
        Complete configuration for all SV triage stages.

    Raises
    ------
    FileNotFoundError
        If the VCF or any required reference file is missing.
    ValueError
        If configuration parameters are invalid.
    """

    def __init__(self, config: SVTriageConfig) -> None:
        self._config = config
        self._validate_config()

    def _validate_config(self) -> None:
        """Fail-fast on missing files and invalid configuration."""
        if not self._config.vcf_path.exists():
            raise FileNotFoundError(
                f"Input VCF not found: {self._config.vcf_path}"
            )

        if self._config.gene_annotation_path is not None:
            if not self._config.gene_annotation_path.exists():
                raise FileNotFoundError(
                    f"Gene annotation not found: {self._config.gene_annotation_path}"
                )

        if self._config.dosage_sensitivity_path is not None:
            if not self._config.dosage_sensitivity_path.exists():
                raise FileNotFoundError(
                    f"Dosage sensitivity file not found: "
                    f"{self._config.dosage_sensitivity_path}"
                )

        if self._config.gnomad_sv_path is not None:
            if not self._config.gnomad_sv_path.exists():
                raise FileNotFoundError(
                    f"gnomAD-SV file not found: {self._config.gnomad_sv_path}"
                )

    def run(self) -> Path:
        """Execute the full SV triage pipeline.

        Returns
        -------
        Path
            Path to the generated report file.
        """
        logger.info(
            "Starting SV triage pipeline: %s → %s",
            self._config.vcf_path,
            self._config.output_path,
        )

        # Build stages
        annotator = SVAnnotator(
            gene_annotation_path=self._config.gene_annotation_path,
            dosage_sensitivity_path=self._config.dosage_sensitivity_path,
            gnomad_sv_path=self._config.gnomad_sv_path,
            reciprocal_overlap=self._config.reciprocal_overlap,
            whole_gene_threshold=self._config.whole_gene_threshold,
        )

        scorer = SVScorer(
            max_allele_frequency=self._config.max_allele_frequency,
        )

        pathogenic_regions, region_names = self._load_regions_with_names(
            self._config.pathogenic_regions_path
        )
        benign_regions = self._load_regions(self._config.benign_regions_path)

        classifier = SVClassifier(
            pathogenic_regions=pathogenic_regions,
            benign_regions=benign_regions,
            pathogenic_region_names=region_names,
        )

        # Stream through the pipeline
        with SVParser(
            self._config.vcf_path,
            min_size=self._config.min_sv_size,
            max_size=self._config.max_sv_size,
            min_quality=self._config.min_quality,
        ) as parser:
            annotated = annotator.annotate(iter(parser))
            scored = scorer.score(annotated)
            classified = classifier.classify(scored)

            results = self._collect_results(classified)

        # Write report
        self._write_report(results)

        logger.info(
            "SV triage complete: %d variants classified, output at %s",
            len(results),
            self._config.output_path,
        )

        return self._config.output_path

    def _collect_results(
        self, classified: "Iterator[ClassifiedSV]"
    ) -> list[ClassifiedSV]:
        """Collect classified SVs, applying output filtering."""
        results: list[ClassifiedSV] = []

        for sv in classified:
            if not self._config.include_benign and sv.classification in (
                SVClassification.BENIGN,
                SVClassification.LIKELY_BENIGN,
            ):
                continue
            results.append(sv)

        # Sort by pathogenicity score descending (None last)
        results.sort(
            key=lambda s: (
                s.scored.pathogenicity_score is None,
                -(s.scored.pathogenicity_score or 0.0),
            )
        )

        return results

    def _write_report(self, results: list[ClassifiedSV]) -> None:
        """Write results to the configured output format."""
        output_path = self._config.output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if self._config.output_format == "json":
            self._write_json(results, output_path)
        elif self._config.output_format == "csv":
            self._write_csv(results, output_path)

    def _write_json(self, results: list[ClassifiedSV], path: Path) -> None:
        """Serialize results as JSON."""
        records = [self._sv_to_dict(sv) for sv in results]

        with open(path, "w") as fh:
            json.dump(
                {
                    "pipeline": "structural_variant_triage",
                    "version": "0.13.0",
                    "total_variants": len(records),
                    "variants": records,
                },
                fh,
                indent=2,
            )

    def _write_csv(self, results: list[ClassifiedSV], path: Path) -> None:
        """Serialize results as CSV."""
        fieldnames = [
            "chrom",
            "start",
            "end",
            "sv_type",
            "length",
            "consequence",
            "classification",
            "pathogenicity_score",
            "evidence_score",
            "genes_affected",
            "hi_genes_affected",
            "population_frequency",
            "evidence_categories",
            "gene_symbols",
        ]

        with open(path, "w", newline="") as fh:
            writer = csv_module.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for sv in results:
                writer.writerow(self._sv_to_flat_dict(sv))

    def _sv_to_dict(self, classified: ClassifiedSV) -> dict[str, Any]:
        """Convert a ClassifiedSV to a JSON-serializable dict."""
        sv = classified.scored.annotated.sv
        annotated = classified.scored.annotated

        gene_details = []
        for overlap in annotated.gene_overlaps:
            gene_details.append({
                "symbol": overlap.gene_symbol,
                "overlap_fraction": round(overlap.overlap_fraction, 4),
                "is_whole_gene": overlap.is_whole_gene,
                "exons_affected": overlap.exons_affected,
                "total_exons": overlap.total_exons,
                "is_haploinsufficient": overlap.is_haploinsufficient,
                "is_triplosensitive": overlap.is_triplosensitive,
                "hi_score": overlap.hi_score,
                "ts_score": overlap.ts_score,
            })

        return {
            "chrom": sv.chrom,
            "start": sv.start,
            "end": sv.end,
            "sv_type": sv.sv_type.value,
            "length": sv.length,
            "id": sv.id,
            "quality": sv.qual,
            "filter": sv.filter_status,
            "copy_number": sv.copy_number,
            "consequence": annotated.consequence.value,
            "classification": classified.classification.value,
            "pathogenicity_score": (
                round(classified.scored.pathogenicity_score, 4)
                if classified.scored.pathogenicity_score is not None
                else None
            ),
            "evidence_score": round(classified.evidence_score, 4),
            "evidence_categories": sorted(
                cat.value for cat in classified.evidence_categories
            ),
            "genes_affected": annotated.genes_affected,
            "hi_genes_affected": annotated.hi_genes_affected,
            "population_frequency": annotated.population_frequency,
            "frequency_unknown": annotated.frequency_unknown,
            "gene_details": gene_details,
            "scoring": {
                "dosage_score": (
                    round(classified.scored.dosage_score, 4)
                    if classified.scored.dosage_score is not None
                    else None
                ),
                "size_score": round(classified.scored.size_score, 4),
                "frequency_score": round(classified.scored.frequency_score, 4),
            },
            "missing_data_sources": sorted(classified.missing_data_sources),
            "syndrome_name": classified.syndrome_name,
        }

    def _sv_to_flat_dict(self, classified: ClassifiedSV) -> dict[str, Any]:
        """Convert to a flat dict for CSV output."""
        sv = classified.scored.annotated.sv
        annotated = classified.scored.annotated

        gene_symbols = ",".join(
            o.gene_symbol for o in annotated.gene_overlaps
        )
        evidence_cats = ",".join(
            sorted(cat.value for cat in classified.evidence_categories)
        )

        return {
            "chrom": sv.chrom,
            "start": sv.start,
            "end": sv.end,
            "sv_type": sv.sv_type.value,
            "length": sv.length,
            "consequence": annotated.consequence.value,
            "classification": classified.classification.value,
            "pathogenicity_score": (
                round(classified.scored.pathogenicity_score, 4)
                if classified.scored.pathogenicity_score is not None
                else ""
            ),
            "evidence_score": round(classified.evidence_score, 4),
            "genes_affected": annotated.genes_affected,
            "hi_genes_affected": annotated.hi_genes_affected,
            "population_frequency": (
                annotated.population_frequency
                if annotated.population_frequency is not None
                else ""
            ),
            "evidence_categories": evidence_cats,
            "gene_symbols": gene_symbols,
        }

    def _load_regions(
        self, path: Optional[Path]
    ) -> list[tuple[str, int, int]]:
        """Load BED regions (chrom, start, end) from a file."""
        if path is None or not path.exists():
            return []

        regions: list[tuple[str, int, int]] = []

        with open(path, "r") as fh:
            for line in fh:
                if line.startswith(("#", "track")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    regions.append((chrom, start, end))
                except ValueError:
                    continue

        logger.info("Loaded %d regions from %s", len(regions), path)
        return regions

    def _load_regions_with_names(
        self, path: Optional[Path]
    ) -> tuple[list[tuple[str, int, int]], dict[tuple[str, int, int], str]]:
        """Load BED regions with optional 4th-column syndrome names."""
        if path is None or not path.exists():
            return [], {}

        regions: list[tuple[str, int, int]] = []
        names: dict[tuple[str, int, int], str] = {}

        with open(path, "r") as fh:
            for line in fh:
                if line.startswith(("#", "track")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue

                key = (chrom, start, end)
                regions.append(key)

                if len(parts) >= 4 and parts[3].strip():
                    names[key] = parts[3].strip()

        logger.info("Loaded %d regions from %s", len(regions), path)
        return regions, names
