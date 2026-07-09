"""Shared pytest fixtures and configuration for the test suite."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Generator

import pytest
from hypothesis import settings

from vartriage.models.config import (
    QualityFilterConfig,
    PrioritizationConfig,
)


settings.register_profile("ci", max_examples=500, deadline=None)
settings.register_profile("dev", max_examples=50, deadline=None)
settings.register_profile("debug", max_examples=10, deadline=None)
settings.load_profile("dev")


@pytest.fixture
def tmp_dir() -> Generator[Path, None, None]:
    """Provide a temporary directory cleaned up after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def default_quality_config() -> QualityFilterConfig:
    """Standard quality filter configuration with default threshold."""
    return QualityFilterConfig(min_qual=20.0)


@pytest.fixture
def default_prioritization_config() -> PrioritizationConfig:
    """Standard prioritization configuration with default thresholds."""
    return PrioritizationConfig(max_allele_frequency=0.01)
