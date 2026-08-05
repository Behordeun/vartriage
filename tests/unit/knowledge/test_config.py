"""Unit tests for KnowledgeBaseConfig validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.knowledge.config import KnowledgeBaseConfig


def test_valid_hpo_terms_accepted() -> None:
    config = KnowledgeBaseConfig(hpo_terms=frozenset({"HP:0001250", "HP:0001249"}))
    assert len(config.hpo_terms) == 2


def test_invalid_hpo_term_format_raises() -> None:
    with pytest.raises(ValueError, match="HP:NNNNNNN"):
        KnowledgeBaseConfig(hpo_terms=frozenset({"INVALID"}))


def test_short_hpo_term_raises() -> None:
    with pytest.raises(ValueError, match="HP:NNNNNNN"):
        KnowledgeBaseConfig(hpo_terms=frozenset({"HP:123"}))


def test_empty_hpo_terms_valid() -> None:
    config = KnowledgeBaseConfig(hpo_terms=frozenset())
    assert config.hpo_terms == frozenset()


def test_default_data_dir_resolves(tmp_path: Path) -> None:
    config = KnowledgeBaseConfig(data_dir=tmp_path)
    assert config.resolved_data_dir == tmp_path


def test_none_data_dir_uses_bundled_default() -> None:
    config = KnowledgeBaseConfig()
    resolved = config.resolved_data_dir
    # Should point to vartriage/data/knowledge/
    assert resolved.name == "knowledge"
    assert resolved.parent.name == "data"
