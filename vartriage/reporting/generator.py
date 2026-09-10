"""Report generation: routes output to JSON, CSV, or PDF writers.

Writes to a temp file first and does an atomic rename on success, so the
target path never contains partial output if something fails mid-write.

Accepts both iterators and sequences. JSON/CSV stream directly; PDF
materializes everything since page layout needs random access.
"""

from __future__ import annotations

import contextlib
import os
import tempfile
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

from vartriage._internal.path_safety import resolve_path
from vartriage.models.config import ClinicalReportConfig, ReportConfig
from vartriage.models.variant import ClassifiedVariant
from vartriage.reporting.csv_writer import write_csv
from vartriage.reporting.json_writer import write_json


class ReportGenerator:
    """Writes clinical reports in JSON, CSV, PDF, or clinical formats.

    Uses a temp file + atomic rename so the output path is never left
    in a half-written state. JSON and CSV stream from iterators without
    buffering; PDF materializes all variants for pagination. Clinical
    formats (clinical-pdf, clinical-html, clinical-docx) delegate to
    the ClinicalReportGenerator.

    Parameters
    ----------
    config : ReportConfig
        Settings including the desired output format.
    clinical_config : ClinicalReportConfig, optional
        Configuration for clinical report generation. Required when
        the output format starts with "clinical-".
    reference_checksums : dict[str, str], optional
        Mapping of reference file paths to SHA-256 checksums.
        Passed through to ClinicalReportGenerator for the audit
        trail.
    """

    def __init__(
        self,
        config: ReportConfig,
        clinical_config: ClinicalReportConfig | None = None,
        reference_checksums: dict[str, str] | None = None,
    ) -> None:
        self._config = config
        self._clinical_config = clinical_config
        self._reference_checksums: dict[str, str] = (
            reference_checksums if reference_checksums is not None else {}
        )

    def generate(
        self,
        variants: Iterator[ClassifiedVariant] | Sequence[ClassifiedVariant],
        output_path: Path,
        source_vcf_path: Path | None = None,
        mito_results: list[Any] | None = None,
        qc_report: Any = None,
    ) -> Path:
        """Write classified variants to the configured format.

        Writes to a temp file alongside the target, then atomically
        replaces it on success. On failure, cleans up the temp file.

        For VCF format, delegates entirely to ``write_vcf`` which
        handles its own atomic write and tabix indexing internally.

        Parameters
        ----------
        variants : Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]]
            Classified variants in priority order. May be empty.
            JSON/CSV consume iterators incrementally; PDF materializes.
        output_path : Path
            Where the final report lands.
        source_vcf_path : Path, optional
            Path to the original input VCF. Required when format
            is ``"vcf"``; ignored for other formats.

        Returns
        -------
        Path
            The written report path.

        Raises
        ------
        IOError
            On write or encoding failure, or if format is "vcf"
            and ``source_vcf_path`` is None.
        """
        output_path = resolve_path(Path(output_path))
        output_path.parent.mkdir(parents=True, exist_ok=True)

        fmt = self._config.output_format

        if fmt.startswith("clinical-"):
            return self._generate_clinical(
                variants,
                output_path,
                mito_results=mito_results,
                qc_report=qc_report,
            )

        if fmt == "vcf":
            if source_vcf_path is None:
                raise OSError("VCF output format requires source_vcf_path")
            from vartriage.reporting.vcf_writer import write_vcf

            # VCF output materializes all variants into memory for the
            # lookup dict, unlike JSON/CSV which stream incrementally.
            # For whole-genome inputs this can reach hundreds of MB.
            materialized = list(variants)
            write_vcf(materialized, source_vcf_path, output_path)
            return output_path

        tmp_fd = None
        tmp_path: Path | None = None
        try:
            tmp_fd, tmp_name = tempfile.mkstemp(
                dir=output_path.parent,
                prefix=".report_",
                suffix=f".{fmt}.tmp",
            )
            os.close(tmp_fd)
            tmp_fd = None
            tmp_path = Path(tmp_name)

            if fmt == "json":
                write_json(variants, tmp_path)
            elif fmt == "csv":
                write_csv(variants, tmp_path)
            elif fmt == "pdf":
                materialized = list(variants)
                self._write_pdf(materialized, tmp_path)
            else:
                raise OSError(f"Unsupported output format: {fmt}")

            os.replace(str(tmp_path), str(output_path))
            tmp_path = None

            return output_path

        except OSError:
            raise
        except Exception as exc:
            raise OSError(f"Failed to generate {fmt.upper()} report: {exc}") from exc
        finally:
            if tmp_path is not None:
                with contextlib.suppress(OSError):
                    tmp_path.unlink(missing_ok=True)

    def _write_pdf(
        self,
        variants: Sequence[ClassifiedVariant],
        output_path: Path,
    ) -> Path:
        """Render PDF, trying reportlab first then the text fallback.

        Parameters
        ----------
        variants : Sequence[ClassifiedVariant]
            Materialized variant list.
        output_path : Path
            Temp path for the PDF.

        Returns
        -------
        Path
            Path to the rendered PDF.

        Raises
        ------
        IOError
            If no PDF backend is available or rendering fails.
        """
        try:
            from vartriage.reporting.pdf_writer import (
                HAS_REPORTLAB,
                ReportlabPDFRenderer,
            )
        except ImportError:
            pass
        else:
            if HAS_REPORTLAB:
                renderer = ReportlabPDFRenderer()
                return renderer.render(list(variants), output_path)

        from vartriage.reporting.pdf_fallback import PDFFallbackRenderer

        renderer_fallback = PDFFallbackRenderer()
        return renderer_fallback.render(list(variants), output_path)

    def _generate_clinical(
        self,
        variants: Iterator[ClassifiedVariant] | Sequence[ClassifiedVariant],
        output_path: Path,
        mito_results: list[Any] | None = None,
        qc_report: Any = None,
    ) -> Path:
        """Delegate to ClinicalReportGenerator for clinical formats.

        Parameters
        ----------
        variants : Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]]
            Classified variants to include in the report.
        output_path : Path
            Where the final report lands.

        Returns
        -------
        Path
            The written report path.

        Raises
        ------
        IOError
            If clinical_config was not provided.
        """
        if self._clinical_config is None:
            raise OSError(
                "Clinical format requires a ClinicalReportConfig to be provided"
            )

        from vartriage import __version__
        from vartriage.reporting.clinical.generator import ClinicalReportGenerator

        clinical_gen = ClinicalReportGenerator(
            config=self._clinical_config,
            pipeline_version=__version__,
            reference_checksums=self._reference_checksums,
        )
        return clinical_gen.generate(
            variants,
            output_path,
            mito_results=mito_results,
            qc_report=qc_report,
        )
