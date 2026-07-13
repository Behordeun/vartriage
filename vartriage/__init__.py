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

from vartriage.annotation.engine import AnnotationEngine
from vartriage.classification.acmg import ACMGClassifier
from vartriage.exceptions import VarTriageWarning
from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.exceptions import (ConfigurationError, ParseError,
                                     ReferenceFileError,
                                     VariantPrioritizationError)
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.config import (AnnotationConfig, MissingDataConfig,
                                     PipelineConfig, PrioritizationConfig,
                                     QualityFilterConfig, ReportConfig)
from vartriage.models.variant import (EVIDENCE_STRENGTH_MAP,
                                      ACMGClassification, AnnotatedVariant,
                                      ClassifiedVariant, ClinVarAssertion,
                                      EvidenceStrength, EvidenceTag,
                                      FunctionalConsequence, ScoredVariant,
                                      Variant)
from vartriage.models.warnings import MissingDataWarning
from vartriage.pipeline import Pipeline
from vartriage.prioritization.engine import PrioritizationEngine
from vartriage.reporting.generator import ReportGenerator

try:
    from importlib.metadata import version as _get_version

    __version__ = _get_version("vartriage")
except Exception:
    __version__ = "0.1.0"

__all__ = [
    # Pipeline orchestrator
    "Pipeline",
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
]
