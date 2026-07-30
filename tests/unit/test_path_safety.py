"""Unit tests for path_safety traversal protection and storage-path constraints."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from vartriage._internal.path_safety import safe_read_path, safe_write_path


class TestSafeReadPath:
    def test_rejects_direct_traversal(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            safe_read_path(tmp_path / ".." / "secret.txt")

    def test_rejects_embedded_traversal(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            safe_read_path(tmp_path / "nested" / ".." / ".." / "secret.txt")

    def test_raises_file_not_found_for_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            safe_read_path(tmp_path / "nonexistent.txt")

    def test_returns_resolved_path_for_existing_file(self, tmp_path: Path) -> None:
        f = tmp_path / "file.txt"
        f.write_text("data")
        result = safe_read_path(f)
        assert result.is_absolute()
        assert result.exists()


class TestSafeWritePath:
    def test_rejects_direct_traversal(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            safe_write_path(tmp_path / ".." / "outside.txt")

    def test_rejects_embedded_traversal(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            safe_write_path(tmp_path / "a" / ".." / ".." / "outside.txt")

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        target = safe_write_path(tmp_path / "subdir" / "nested" / "file.txt")
        assert target.is_absolute()
        assert target.parent.is_dir()

    def test_returned_path_is_writable(self, tmp_path: Path) -> None:
        target = safe_write_path(tmp_path / "out" / "file.txt")
        target.write_text("content")
        assert target.read_text() == "content"
