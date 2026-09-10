"""Variant Prioritization Library.

A streaming pipeline for identifying and classifying pathogenic genetic
variants from VCF data. Designed for whole-genome scale datasets with
memory-bounded processing (<2GB RSS for 4M+ variants).

Example
-------
>>> from vartriage import Pipeline, PipelineConfig
>>> config = PipelineConfig(vcf_path=Path("input.vcf.gz"), output_path=Path("report.json"))
>>> pipeline = Pipeline(config)
>>> output = pipeline.run(vcf_path=config.vcf_path)
"""

# PS1/PM5 protein index (v0.14.0)
from vartriage.annotation.clinvar_protein_index import ClinVarProteinIndex
from vartriage.annotation.engine import AnnotationEngine
from vartriage.classification.acmg import ACMGClassifier

# Cohort analysis (v0.11.0)
from vartriage.cohort.aggregator import CohortAggregator
from vartriage.cohort.pipeline import CohortPipeline
from vartriage.cohort.report import CohortReportGenerator
from vartriage.cohort.statistics import CohortStatistics
from vartriage.exceptions import VarTriageWarning
from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.exceptions import (
    ConfigurationError,
    ParseError,
    ReferenceFileError,
    VariantPrioritizationError,
)
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.cohort import (
    CohortConfig,
    CohortSummary,
    CohortVariant,
    GeneBurden,
)
from vartriage.models.config import (
    AnnotationConfig,
    MissingDataConfig,
    PipelineConfig,
    PrioritizationConfig,
    QualityFilterConfig,
    ReportConfig,
)
from vartriage.models.variant import (
    EVIDENCE_STRENGTH_MAP,
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    ClinVarAssertion,
    EvidenceStrength,
    EvidenceTag,
    FunctionalConsequence,
    ProteinChange,
    ScoredVariant,
    Variant,
)
from vartriage.models.warnings import MissingDataWarning
from vartriage.pipeline import Pipeline
from vartriage.prioritization.engine import PrioritizationEngine
from vartriage.reporting.generator import ReportGenerator

# Structural variant triage (v0.13.0)
from vartriage.structural import (
    AnnotatedSV,
    ClassifiedSV,
    ScoredSV,
    StructuralVariant,
    SVAnnotator,
    SVClassification,
    SVClassifier,
    SVConsequence,
    SVEvidenceCategory,
    SVParser,
    SVScorer,
    SVTriageConfig,
    SVTriagePipeline,
    SVType,
)

try:
    from importlib.metadata import version as _get_version

    __version__ = _get_version("vartriage")
except Exception:
    __version__ = "0.18.0"

__all__ = [
    # Pipeline orchestrator
    "Pipeline",
    # Cohort analysis
    "CohortPipeline",
    "CohortAggregator",
    "CohortReportGenerator",
    "CohortStatistics",
    "CohortConfig",
    "CohortVariant",
    "CohortSummary",
    "GeneBurden",
    # Structural variant triage
    "SVTriagePipeline",
    "SVParser",
    "SVAnnotator",
    "SVScorer",
    "SVClassifier",
    "SVTriageConfig",
    "StructuralVariant",
    "AnnotatedSV",
    "ScoredSV",
    "ClassifiedSV",
    "SVType",
    "SVConsequence",
    "SVClassification",
    "SVEvidenceCategory",
    # Processing stages
    "VCFParser",
    "QualityFilter",
    "AnnotationEngine",
    "PrioritizationEngine",
    "ACMGClassifier",
    "ReportGenerator",
    # Core data models
    "Variant",
    "AnnotatedVariant",
    "ScoredVariant",
    "ClassifiedVariant",
    # Enums
    "FunctionalConsequence",
    "ClinVarAssertion",
    "ACMGClassification",
    "EvidenceTag",
    "EvidenceStrength",
    # Constants
    "EVIDENCE_STRENGTH_MAP",
    # Configuration classes
    "QualityFilterConfig",
    "AnnotationConfig",
    "PrioritizationConfig",
    "ReportConfig",
    "MissingDataConfig",
    "PipelineConfig",
    # Exceptions
    "VariantPrioritizationError",
    "ParseError",
    "ConfigurationError",
    "ReferenceFileError",
    # Warnings
    "VarTriageWarning",
    "MissingDataWarning",
    # v0.14.0: PS1/PM5 and gnomAD API
    "ProteinChange",
    "ClinVarProteinIndex",
    "GnomADClient",
]


def __getattr__(name: str) -> object:
    """Lazy import for GnomADClient to avoid httpx import at module level."""
    if name == "GnomADClient":
        from vartriage.api.gnomad_client import GnomADClient

        globals()["GnomADClient"] = GnomADClient
        return GnomADClient
    raise AttributeError(f"module 'vartriage' has no attribute {name!r}")
