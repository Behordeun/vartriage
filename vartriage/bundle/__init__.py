"""Score bundle downloader for vartriage reference files.

Provides automated downloading, transformation, verification, and
storage of genomics reference data (gnomAD, ClinVar, CADD, REVEL,
SpliceAI, GENCODE) required by the vartriage pipeline.

Public API
----------
BundleRegistry : Registry of available score bundles
BundleStorage : On-disk layout and path resolution
BundleDownloader : HTTP download engine with resume support
BundleConfig : User configuration (TOML-based)
"""

from vartriage.bundle.config import BundleConfig
from vartriage.bundle.manifest import BundleManifest
from vartriage.bundle.registry import BundleEntry, BundleRegistry
from vartriage.bundle.storage import BundleStorage

__all__ = [
    "BundleConfig",
    "BundleEntry",
    "BundleManifest",
    "BundleRegistry",
    "BundleStorage",
]
