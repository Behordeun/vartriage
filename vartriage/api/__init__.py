"""API-based annotation backend for zero-configuration variant triage.

Provides HTTP clients for Ensembl VEP, ClinVar E-utilities, CADD, and
SpliceAI that implement the same Protocol interfaces as the local file
backends. Requires the `api` optional dependency group (httpx).

Install with: pip install vartriage[api]
"""

from __future__ import annotations

from vartriage.api.config import APIConfig

__all__ = [
    "APIConfig",
    "APIAnnotationEngine",
    "APIScoreProvider",
]


def __getattr__(name: str) -> object:
    """Lazy imports to avoid pulling httpx on module access."""
    if name == "APIAnnotationEngine":
        from vartriage.api.annotation_engine import APIAnnotationEngine

        return APIAnnotationEngine
    if name == "APIScoreProvider":
        from vartriage.api.score_provider import APIScoreProvider

        return APIScoreProvider
    raise AttributeError(f"module 'vartriage.api' has no attribute {name!r}")
