"""Report generation — routes output to JSON, CSV, or PDF writers.

Writes to a temp file first and does an atomic rename on success, so the
target path never contains partial output if something fails mid-write.

Accepts both iterators and sequences. JSON/CSV stream directly; PDF
materializes everything since page layout needs random access.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Iterator, Sequence, Union

from vartriage.models.config import ReportConfig
from vartriage.models.variant import ClassifiedVariant
from vartriage.reporting.csv_writer import write_csv
from vartriage.reporting.json_writer import write_json


class ReportGenerator:
    """Writes clinical reports in JSON, CSV, or PDF.

    Uses a temp file + atomic rename so the output path is never left
    in a half-written state. JSON and CSV stream from iterators without
    buffering; PDF materializes all variants for pagination.

    Parameters
    ----------
    config : ReportConfig
        Settings including the desired output format.
    """

    def __init__(self, config: ReportConfig) -> None:
        self._config = config

    def generate(
        self,
        variants: Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]],
        output_path: Path,
    ) -> Path:
        """Write classified variants to the configured format.

        Writes to a temp file alongside the target, then atomically
        replaces it on success. On failure, cleans up the temp file.

        Parameters
        ----------
        variants : Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]]
            Classified variants in priority order. May be empty.
            JSON/CSV consume iterators incrementally; PDF materializes.
        output_path : Path
            Where the final report lands.

        Returns
        -------
        Path
            The written report path.

        Raises
        ------
        IOError
            On write or encoding failure.
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
                materialized = list(variants)
                self._write_pdf(materialized, tmp_path)
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

        from vartriage.reporting.pdf_fallback import (
            PDFFallbackRenderer,
        )

        renderer_fallback = PDFFallbackRenderer()
        return renderer_fallback.render(list(variants), output_path)
