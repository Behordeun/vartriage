"""Progress display for downloads and transforms.

Writes to stderr so it doesn't interfere with pipeline output.
No external dependencies - uses simple ANSI escape codes.
"""

from __future__ import annotations

import sys
import time
from typing import Any, Optional

from vartriage.bundle._disk import format_bytes


class ProgressBar:
    """Simple stderr progress bar for file downloads.

    Displays: filename, bytes transferred, speed, percentage, ETA.
    """

    def __init__(
        self,
        filename: str,
        total_bytes: Optional[int] = None,
        stream: "Any" = None,
    ) -> None:
        self._filename = filename
        self._total = total_bytes
        self._transferred = 0
        self._start_time = time.monotonic()
        self._last_update = 0.0
        self._stream = stream or sys.stderr

    def update(self, bytes_written: int) -> None:
        """Update progress with new bytes written.

        Limits output to 10 updates per second to avoid flooding.
        """
        self._transferred += bytes_written
        now = time.monotonic()

        # Rate-limit display updates to every 100ms
        if now - self._last_update < 0.1:
            return
        self._last_update = now

        self._render()

    def finish(self) -> None:
        """Mark download complete and render final state."""
        self._render(final=True)
        self._write("\n")

    def _render(self, final: bool = False) -> None:
        elapsed = time.monotonic() - self._start_time
        speed = self._transferred / elapsed if elapsed > 0 else 0

        parts = [f"\r  {self._filename}: {format_bytes(self._transferred)}"]

        if self._total and self._total > 0:
            pct = min(100.0, (self._transferred / self._total) * 100)
            parts.append(f" / {format_bytes(self._total)} ({pct:.0f}%)")

            if speed > 0 and not final:
                remaining = (self._total - self._transferred) / speed
                parts.append(
                    f" [{format_bytes(int(speed))}/s, {self._format_eta(remaining)}]"
                )
        else:
            if speed > 0:
                parts.append(f" [{format_bytes(int(speed))}/s]")

        if final:
            parts.append(" done")

        # Clear rest of line
        line = "".join(parts)
        self._write(f"{line}\033[K")

    def _format_eta(self, seconds: float) -> str:
        if seconds < 60:
            return f"{seconds:.0f}s"
        elif seconds < 3600:
            return f"{seconds / 60:.0f}m"
        else:
            return f"{seconds / 3600:.1f}h"

    def _write(self, text: str) -> None:
        try:
            self._stream.write(text)
            self._stream.flush()
        except (OSError, AttributeError):
            pass
