"""HTTP download engine with resume support.

Uses urllib from stdlib - no external dependencies required.
Falls back gracefully when servers don't support Range requests.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from vartriage.bundle._checksums import compute_sha256, verify_checksum
from vartriage.bundle._disk import check_disk_space
from vartriage.bundle._progress import ProgressBar


@dataclass
class DownloadResult:
    """Result of a single file download.

    Attributes
    ----------
    path : Path
        Final path of the downloaded file.
    bytes_downloaded : int
        Total bytes transferred in this session.
    resumed : bool
        True if download resumed from a partial file.
    duration_seconds : float
        Wall-clock time for the download.
    checksum_verified : bool
        True if checksum was verified after download.
    """

    path: Path
    bytes_downloaded: int
    resumed: bool
    duration_seconds: float
    checksum_verified: bool


class DownloadError(Exception):
    """Raised when a download fails after all retries."""

    def __init__(self, url: str, reason: str) -> None:
        self.url = url
        self.reason = reason
        super().__init__(f"Download failed for {url}: {reason}")


class BundleDownloader:
    """HTTP file downloader with resume, retry, and progress display.

    Features:
    - Resume partial downloads via HTTP Range header
    - Atomic completion (.partial file, rename on success)
    - Exponential backoff retry for transient errors
    - Streaming SHA-256 checksum (no second pass needed)
    - Progress bar on stderr
    """

    # HTTP status codes that warrant a retry
    _RETRYABLE_STATUS = {429, 500, 502, 503, 504}

    def __init__(
        self,
        timeout: tuple[int, int] = (30, 60),
        max_retries: int = 3,
        show_progress: bool = True,
    ) -> None:
        """Initialize the downloader.

        Parameters
        ----------
        timeout : tuple[int, int]
            (connect_timeout, read_timeout) in seconds.
        max_retries : int
            Max retry attempts for transient failures.
        show_progress : bool
            Whether to display a progress bar on stderr.
        """
        self._connect_timeout = timeout[0]
        self._read_timeout = timeout[1]
        self._max_retries = max_retries
        self._show_progress = show_progress

    def download(
        self,
        url: str,
        dest: Path,
        expected_size: Optional[int] = None,
        expected_checksum: str = "",
        resume: bool = True,
    ) -> DownloadResult:
        """Download a single file with resume support.

        Parameters
        ----------
        url : str
            URL to download.
        dest : Path
            Final destination path.
        expected_size : int, optional
            Expected file size for progress display and disk check.
        expected_checksum : str
            SHA-256 checksum to verify after download. Empty to skip.
        resume : bool
            Attempt to resume partial downloads.

        Returns
        -------
        DownloadResult
            Result with download metadata.

        Raises
        ------
        DownloadError
            If download fails after all retries.
        OSError
            If disk space is insufficient.
        """
        dest = dest.resolve()
        dest.parent.mkdir(parents=True, exist_ok=True)

        partial_path = Path(str(dest) + ".partial")
        resumed = False
        start_offset = 0

        # Check disk space if we know the expected size
        if expected_size and expected_size > 0:
            check_disk_space(dest.parent, expected_size)

        # Resume from partial if it exists
        if resume and partial_path.exists():
            start_offset = partial_path.stat().st_size
            resumed = True

        start_time = time.monotonic()
        bytes_downloaded = self._download_with_retry(
            url=url,
            dest=partial_path,
            start_offset=start_offset,
            expected_size=expected_size,
        )
        duration = time.monotonic() - start_time

        # Atomic rename to final destination
        os.replace(str(partial_path), str(dest))

        # Verify checksum
        checksum_ok = True
        if expected_checksum:
            checksum_ok = verify_checksum(dest, expected_checksum)
            if not checksum_ok:
                actual = compute_sha256(dest)
                dest.unlink(missing_ok=True)
                raise DownloadError(
                    url,
                    f"Checksum mismatch: expected {expected_checksum}, got {actual}",
                )

        return DownloadResult(
            path=dest,
            bytes_downloaded=bytes_downloaded,
            resumed=resumed,
            duration_seconds=duration,
            checksum_verified=bool(expected_checksum),
        )

    def _download_with_retry(
        self,
        url: str,
        dest: Path,
        start_offset: int,
        expected_size: Optional[int],
    ) -> int:
        """Download with exponential backoff retry."""
        last_error: Optional[Exception] = None

        for attempt in range(self._max_retries + 1):
            if attempt > 0:
                wait = min(2**attempt, 8)
                time.sleep(wait)

            try:
                return self._stream_download(
                    url=url,
                    dest=dest,
                    start_offset=start_offset,
                    expected_size=expected_size,
                )
            except HTTPError as exc:
                if exc.code not in self._RETRYABLE_STATUS:
                    raise DownloadError(url, f"HTTP {exc.code}: {exc.reason}") from exc
                last_error = exc
            except OSError as exc:
                last_error = exc

            start_offset = self._current_offset(dest, start_offset)

        raise DownloadError(
            url, f"Failed after {self._max_retries + 1} attempts: {last_error}"
        )

    @staticmethod
    def _current_offset(dest: Path, fallback: int) -> int:
        """Return current file size for resume, or fallback if file missing."""
        return dest.stat().st_size if dest.exists() else fallback

    def _stream_download(
        self,
        url: str,
        dest: Path,
        start_offset: int,
        expected_size: Optional[int],
    ) -> int:
        """Stream download to disk with optional resume."""
        request = Request(url)
        request.add_header("User-Agent", "vartriage-bundle-downloader/0.6.0")

        mode = "ab"
        if start_offset > 0:
            request.add_header("Range", f"bytes={start_offset}-")
        else:
            mode = "wb"

        response = urlopen(request, timeout=self._connect_timeout)

        # Check if server supports range (206 Partial Content)
        if start_offset > 0 and response.status != 206:
            # Server doesn't support resume, restart from beginning
            start_offset = 0
            mode = "wb"

        # Set up progress display
        progress: Optional[ProgressBar] = None
        if self._show_progress:
            filename = dest.name.replace(".partial", "")
            progress = ProgressBar(filename, expected_size)
            if start_offset > 0:
                progress.update(start_offset)

        bytes_downloaded = 0
        chunk_size = 65536

        with open(dest, mode) as f:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                bytes_downloaded += len(chunk)
                if progress:
                    progress.update(len(chunk))

        if progress:
            progress.finish()

        return bytes_downloaded
