"""Report generation orchestrator.

Routes classified variant output to the appropriate format writer (JSON, CSV,
or PDF) based on configuration. Writes to a temporary file first and performs
an atomic rename on success, ensuring no partial or corrupted output reaches
the target path on failure.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Sequence

from vartriage.models.config import ReportConfig
from vartriage.models.variant import ClassifiedVariant
from vartriage.reporting.csv_writer import write_csv
from vartriage.reporting.json_writer import write_json


class ReportGenerator:
    """Generate clinical reports in JSON, CSV, or PDF format.

    Routes to the appropriate format writer based on ``ReportConfig``.
    All writes go through a temporary file with an atomic rename on
    success, guaranteeing that the target path never holds partial output.

    Parameters
    ----------
    config : ReportConfig
        Report generation settings including output format.

    Raises
    ------
    IOError
        If write fails. Partial or corrupted files are never produced
        at the target path.
    """

    def __init__(self, config: ReportConfig) -> None:
        self._config = config

    def generate(
        self,
        variants: Sequence[ClassifiedVariant],
        output_path: Path,
    ) -> Path:
        """Serialize classified variants to the configured output format.

        Writes to a temporary file in the same directory as the target,
        then atomically moves the temp file to the final path on success.
        On failure, the temp file is cleaned up and no partial output
        exists at the target path.

        Parameters
        ----------
        variants : Sequence[ClassifiedVariant]
            The classified variants to include in the report, in
            prioritized rank order. May be empty.
        output_path : Path
            Destination file path for the report output.

        Returns
        -------
        Path
            The path to the written report file.

        Raises
        ------
        IOError
            If a write error or encoding error prevents report generation.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        fmt = self._config.output_format

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
                self._write_pdf(variants, tmp_path)
            else:
                raise IOError(
                    f"Unsupported output format: {fmt}"
                )

            os.replace(str(tmp_path), str(output_path))
            tmp_path = None

            return output_path

        except IOError:
            raise
        except Exception as exc:
            raise IOError(
                f"Failed to generate {fmt.upper()} report: {exc}"
            ) from exc
        finally:
            if tmp_path is not None:
                try:
                    tmp_path.unlink(missing_ok=True)
                except OSError:
                    pass

    def _write_pdf(
        self,
        variants: Sequence[ClassifiedVariant],
        output_path: Path,
    ) -> Path:
        """Write PDF output, auto-detecting reportlab availability.

        Parameters
        ----------
        variants : Sequence[ClassifiedVariant]
            Variants to render.
        output_path : Path
            Temporary path for the PDF file.

        Returns
        -------
        Path
            The path to the written PDF.

        Raises
        ------
        IOError
            If reportlab is not installed or rendering fails.
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

        from vartriage.reporting.pdf_fallback import (
            PDFFallbackRenderer,
        )

        renderer_fallback = PDFFallbackRenderer()
        return renderer_fallback.render(list(variants), output_path)
