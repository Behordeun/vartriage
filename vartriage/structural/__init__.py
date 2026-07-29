"""Structural variant triage module.

Provides a complete pipeline for parsing, annotating, scoring, and
classifying structural variants (deletions, duplications, inversions,
insertions, translocations, CNVs) using the ClinGen 2020 technical
standards for CNV interpretation.

Example
-------
>>> from vartriage.structural import SVTriagePipeline, SVTriageConfig
>>> config = SVTriageConfig(
...     vcf_path=Path("sv_calls.vcf.gz"),
...     output_path=Path("sv_report.json"),
...     gene_annotation_path=Path("gencode.v44.gtf"),
...     dosage_sensitivity_path=Path("clingen_dosage.tsv"),
... )
>>> pipeline = SVTriagePipeline(config)
>>> output = pipeline.run()
"""

from vartriage.structural.annotator import SVAnnotator
from vartriage.structural.classifier import SVClassifier
from vartriage.structural.config import SVTriageConfig
from vartriage.structural.models import (
    AnnotatedSV,
    Breakpoint,
    ClassifiedSV,
    GeneOverlap,
    ScoredSV,
    StructuralVariant,
    SVClassification,
    SVConsequence,
    SV_CONSEQUENCE_SEVERITY,
    SVEvidenceCategory,
    SV_EVIDENCE_POINTS,
    SVType,
)
from vartriage.structural.parser import SVParser
from vartriage.structural.pipeline import SVTriagePipeline
from vartriage.structural.scoring import SVScorer

__all__ = [
    # Pipeline
    "SVTriagePipeline",
    # Processing stages
    "SVParser",
    "SVAnnotator",
    "SVScorer",
    "SVClassifier",
    # Configuration
    "SVTriageConfig",
    # Data models
    "StructuralVariant",
    "Breakpoint",
    "GeneOverlap",
    "AnnotatedSV",
    "ScoredSV",
    "ClassifiedSV",
    # Enums
    "SVType",
    "SVConsequence",
    "SVClassification",
    "SVEvidenceCategory",
    # Constants
    "SV_CONSEQUENCE_SEVERITY",
    "SV_EVIDENCE_POINTS",
]
