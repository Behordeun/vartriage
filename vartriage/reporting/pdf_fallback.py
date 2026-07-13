"""Fallback PDF renderer when reportlab is not installed.

Stub PDFRenderer that raises ImportError with clear installation
instructions when PDF output is requested without the reportlab dependency.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


class PDFFallbackRenderer:
    """Stub PDF renderer that raises ImportError with install instructions.

    This class satisfies the PDFRenderer protocol interface but always raises
    an ImportError directing users to install the [pdf] extra.

    Methods
    -------
    render(variants, output_path)
        Always raises ImportError with installation instructions.
    """

    def render(self, variants: list[Any], output_path: Path) -> Path:
        """Raise ImportError with instructions to install reportlab.

        Parameters
        ----------
        variants : list[Any]
            List of ClassifiedVariant instances (unused).
        output_path : Path
            Target output path (unused).

        Returns
        -------
        Path
            Never returns. Always raises.

        Raises
        ------
        ImportError
            Always raised with installation instructions.
        """
        raise ImportError(
            "PDF output requires reportlab. " "Install with: pip install vartriage[pdf]"
        )
