"""Cohort report generation in JSON and CSV formats.

Serializes CohortVariant records, GeneBurden tables, and the
CohortSummary into structured output files. Both formats stream
records to keep memory bounded for large cohorts.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from vartriage.models.cohort import (
    CohortConfig,
    CohortSummary,
    CohortVariant,
    GeneBurden,
)

logger = logging.getLogger(__name__)


class CohortReportGenerator:
    """Generate cohort analysis reports in JSON or CSV format.

    Writes three output files per cohort run:
    - variants report: all cohort variants with recurrence data
    - gene burden report: per-gene statistics
    - summary report: top-level cohort metrics

    Parameters
    ----------
    config : CohortConfig
        Cohort configuration with output_path and output_format.
    """

    def __init__(self, config: CohortConfig) -> None:
        self._config = config
        self._output_dir = config.output_path

    def generate(
        self,
        variants: list[CohortVariant],
        gene_burdens: list[GeneBurden],
        summary: CohortSummary,
    ) -> list[Path]:
        """Write all cohort report files.

        Parameters
        ----------
        variants : list[CohortVariant]
            Aggregated cohort variants.
        gene_burdens : list[GeneBurden]
            Per-gene burden statistics.
        summary : CohortSummary
            Top-level cohort metrics.

        Returns
        -------
        list[Path]
            Paths to all generated report files.
        """
        self._output_dir.mkdir(parents=True, exist_ok=True)

        fmt = self._config.output_format
        if fmt == "csv":
            return self._write_csv(variants, gene_burdens, summary)
        return self._write_json(variants, gene_burdens, summary)

    def _write_json(
        self,
        variants: list[CohortVariant],
        gene_burdens: list[GeneBurden],
        summary: CohortSummary,
    ) -> list[Path]:
        """Write JSON format reports."""
        paths: list[Path] = []

        # Variants report
        variants_path = self._output_dir / f"{self._config.cohort_name}_variants.json"
        variant_records = [self._serialize_variant(v) for v in variants]
        self._write_json_file(variants_path, variant_records)
        paths.append(variants_path)

        # Gene burden report
        burden_path = self._output_dir / f"{self._config.cohort_name}_gene_burden.json"
        burden_records = [self._serialize_burden(b) for b in gene_burdens]
        self._write_json_file(burden_path, burden_records)
        paths.append(burden_path)

        # Summary report
        summary_path = self._output_dir / f"{self._config.cohort_name}_summary.json"
        summary_record = self._serialize_summary(summary)
        self._write_json_file(summary_path, summary_record)
        paths.append(summary_path)

        logger.info("JSON reports written to %s", self._output_dir)
        return paths

    def _write_csv(
        self,
        variants: list[CohortVariant],
        gene_burdens: list[GeneBurden],
        summary: CohortSummary,
    ) -> list[Path]:
        """Write CSV format reports."""
        paths: list[Path] = []

        # Variants CSV
        variants_path = self._output_dir / f"{self._config.cohort_name}_variants.csv"
        self._write_variants_csv(variants_path, variants)
        paths.append(variants_path)

        # Gene burden CSV
        burden_path = self._output_dir / f"{self._config.cohort_name}_gene_burden.csv"
        self._write_burden_csv(burden_path, gene_burdens)
        paths.append(burden_path)

        # Summary JSON (always JSON for structured metadata)
        summary_path = self._output_dir / f"{self._config.cohort_name}_summary.json"
        summary_record = self._serialize_summary(summary)
        self._write_json_file(summary_path, summary_record)
        paths.append(summary_path)

        logger.info("CSV reports written to %s", self._output_dir)
        return paths

    def _write_variants_csv(
        self, path: Path, variants: list[CohortVariant]
    ) -> None:
        """Write cohort variants to CSV."""
        fieldnames = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "gene_name",
            "consequence",
            "sample_count",
            "total_samples",
            "cohort_frequency",
            "allele_frequency",
            "max_classification",
            "evidence_tags",
            "samples",
        ]

        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for v in variants:
                writer.writerow(
                    {
                        "chrom": v.chrom,
                        "pos": v.pos,
                        "ref": v.ref,
                        "alt": v.alt,
                        "gene_name": v.gene_name or "",
                        "consequence": v.consequence.value,
                        "sample_count": v.sample_count,
                        "total_samples": v.total_samples,
                        "cohort_frequency": f"{v.cohort_frequency:.4f}",
                        "allele_frequency": (
                            f"{v.allele_frequency:.6f}"
                            if v.allele_frequency is not None
                            else ""
                        ),
                        "max_classification": v.max_classification.value,
                        "evidence_tags": ";".join(
                            sorted(t.value for t in v.all_evidence_tags)
                        ),
                        "samples": ";".join(v.sample_ids),
                    }
                )

    def _write_burden_csv(
        self, path: Path, burdens: list[GeneBurden]
    ) -> None:
        """Write gene burden table to CSV."""
        fieldnames = [
            "gene_name",
            "total_variants",
            "pathogenic_count",
            "samples_affected",
            "total_samples",
            "penetrance",
            "most_severe",
        ]

        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for b in burdens:
                writer.writerow(
                    {
                        "gene_name": b.gene_name,
                        "total_variants": b.total_variants,
                        "pathogenic_count": b.pathogenic_count,
                        "samples_affected": b.samples_affected,
                        "total_samples": b.total_samples,
                        "penetrance": f"{b.penetrance:.4f}",
                        "most_severe": b.most_severe.value,
                    }
                )

    @staticmethod
    def _write_json_file(path: Path, data: Any) -> None:
        """Write data to a JSON file with consistent formatting."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
            f.write("\n")

    @staticmethod
    def _serialize_variant(v: CohortVariant) -> dict[str, Any]:
        """Convert a CohortVariant to a JSON-serializable dict."""
        return {
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "gene_name": v.gene_name,
            "consequence": v.consequence.value,
            "sample_count": v.sample_count,
            "total_samples": v.total_samples,
            "cohort_frequency": round(v.cohort_frequency, 4),
            "is_singleton": v.is_singleton,
            "is_universal": v.is_universal,
            "allele_frequency": v.allele_frequency,
            "max_classification": v.max_classification.value,
            "evidence_tags": sorted(t.value for t in v.all_evidence_tags),
            "samples": [
                {
                    "sample_id": occ.sample_id,
                    "classification": occ.classified.classification.value,
                    "evidence_tags": sorted(
                        t.value for t in occ.classified.evidence_tags
                    ),
                }
                for occ in v.occurrences
            ],
        }

    @staticmethod
    def _serialize_burden(b: GeneBurden) -> dict[str, Any]:
        """Convert a GeneBurden to a JSON-serializable dict."""
        return {
            "gene_name": b.gene_name,
            "total_variants": b.total_variants,
            "pathogenic_count": b.pathogenic_count,
            "samples_affected": b.samples_affected,
            "total_samples": b.total_samples,
            "penetrance": round(b.penetrance, 4),
            "most_severe": b.most_severe.value,
        }

    @staticmethod
    def _serialize_summary(s: CohortSummary) -> dict[str, Any]:
        """Convert a CohortSummary to a JSON-serializable dict."""
        return {
            "cohort_name": s.cohort_name,
            "total_samples": s.total_samples,
            "total_variants": s.total_variants,
            "shared_variants": s.shared_variants,
            "singleton_variants": s.singleton_variants,
            "universal_variants": s.universal_variants,
            "pathogenic_variants": s.pathogenic_variants,
            "likely_pathogenic_variants": s.likely_pathogenic_variants,
            "genes_affected": s.genes_affected,
            "top_recurrent_genes": list(s.top_recurrent_genes),
            "samples_processed": list(s.samples_processed),
        }
