"""Audit trail writer for clinical variant reports.

Produces a JSON sidecar file alongside each clinical report capturing
the full run manifest and per-variant decision log. This supports
reproducibility and automated QC checks.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path
from typing import Sequence

from vartriage.models.config import ClinicalReportConfig
from vartriage.models.variant import (
    ClassifiedVariant,
    EvidenceTag,
)


class AuditTrailWriter:
    """Writes JSON audit sidecar files alongside clinical reports.

    The sidecar contains:
    - run_manifest: all config params, reference checksums,
      pipeline version, timestamp, Python version, OS platform
    - decision_log: per-variant record of evidence tags, scores,
      and classification
    """

    def write(
        self,
        output_path: Path,
        config: ClinicalReportConfig,
        variants: Sequence[ClassifiedVariant],
        reference_checksums: dict[str, str],
        pipeline_version: str,
        execution_timestamp: str,
    ) -> Path:
        """Write the audit JSON sidecar file.

        Parameters
        ----------
        output_path : Path
            Path to the main report file. The sidecar is written
            at ``{output_path}.audit.json``.
        config : ClinicalReportConfig
            Configuration used for this report run.
        variants : Sequence[ClassifiedVariant]
            All classified variants included in the report.
        reference_checksums : dict[str, str]
            Mapping of reference file paths to SHA-256 hex digests.
        pipeline_version : str
            Version string of the vartriage pipeline.
        execution_timestamp : str
            ISO 8601 timestamp of the execution.

        Returns
        -------
        Path
            Path to the written audit sidecar file.

        Raises
        ------
        IOError
            If the sidecar file cannot be written due to filesystem
            errors.
        """
        sidecar_path = Path(str(output_path) + ".audit.json")

        audit_data = {
            "run_manifest": self._build_run_manifest(
                config=config,
                reference_checksums=reference_checksums,
                pipeline_version=pipeline_version,
                execution_timestamp=execution_timestamp,
            ),
            "decision_log": self._build_decision_log(variants),
        }

        try:
            sidecar_path.parent.mkdir(parents=True, exist_ok=True)
            sidecar_path.write_text(
                json.dumps(audit_data, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
        except OSError as exc:
            raise IOError(
                f"Failed to write audit sidecar at "
                f"{sidecar_path}: {exc}"
            ) from exc

        return sidecar_path

    def compute_file_checksum(self, path: Path) -> str:
        """Compute SHA-256 hex digest for a file.

        Parameters
        ----------
        path : Path
            Path to the file to checksum.

        Returns
        -------
        str
            Hex-encoded SHA-256 digest prefixed with "sha256:".

        Raises
        ------
        IOError
            If the file cannot be read.
        """
        sha256 = hashlib.sha256()
        try:
            with open(path, "rb") as f:
                while True:
                    chunk = f.read(65536)
                    if not chunk:
                        break
                    sha256.update(chunk)
        except OSError as exc:
            raise IOError(
                f"Failed to read file for checksum at "
                f"{path}: {exc}"
            ) from exc

        return f"sha256:{sha256.hexdigest()}"

    def _build_run_manifest(
        self,
        config: ClinicalReportConfig,
        reference_checksums: dict[str, str],
        pipeline_version: str,
        execution_timestamp: str,
    ) -> dict:
        """Assemble the run manifest dictionary."""
        return {
            "patient_id": config.patient_id,
            "panel_name": config.panel_name,
            "output_format": config.output_format,
            "report_template": config.report_template,
            "pipeline_version": pipeline_version,
            "execution_timestamp": execution_timestamp,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "reference_files": reference_checksums,
        }

    def _build_decision_log(
        self, variants: Sequence[ClassifiedVariant]
    ) -> list[dict]:
        """Build the per-variant decision log entries."""
        entries: list[dict] = []

        for variant in variants:
            scored = variant.scored
            annotated = scored.annotated
            raw = annotated.variant

            # Determine which tags were skipped due to missing data.
            tags_skipped = self._compute_skipped_tags(variant)

            entry = {
                "chromosome": raw.chrom,
                "position": raw.pos,
                "ref": raw.ref,
                "alt": raw.alt,
                "gene_name": annotated.gene_name,
                "evidence_tags_assigned": sorted(
                    t.value for t in variant.evidence_tags
                ),
                "evidence_tags_skipped": tags_skipped,
                "scores": {
                    "cadd_phred": scored.cadd_phred,
                    "cadd_normalized": scored.cadd_normalized,
                    "revel_score": scored.revel_score,
                    "spliceai_score": scored.spliceai_score,
                    "composite_rank": scored.composite_rank,
                },
                "classification": variant.classification.value,
            }
            entries.append(entry)

        return entries

    def _compute_skipped_tags(
        self, variant: ClassifiedVariant
    ) -> dict[str, list[str]]:
        """Determine which tags were skipped due to missing data.

        Maps tag values to the list of missing source names that
        prevented the tag from being assigned.
        """
        skipped: dict[str, list[str]] = {}
        missing_sources = variant.missing_data_sources

        if not missing_sources:
            return skipped

        # Map known data sources to the tags they could support.
        source_to_tags: dict[str, list[EvidenceTag]] = {
            "REVEL": [EvidenceTag.PP3],
            "SpliceAI": [EvidenceTag.PP3],
            "ClinVar": [EvidenceTag.PP5],
            "gnomAD": [EvidenceTag.PM2],
        }

        for source in sorted(missing_sources):
            affected_tags = source_to_tags.get(source, [])
            for tag in affected_tags:
                if tag not in variant.evidence_tags:
                    tag_key = tag.value
                    if tag_key not in skipped:
                        skipped[tag_key] = []
                    skipped[tag_key].append(source)

        return skipped
