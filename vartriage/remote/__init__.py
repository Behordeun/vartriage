"""Remote tabix score backend for querying large reference databases via HTTP.

Provides CADD and gnomAD frequency lookups using pysam's native support for
remote bgzipped/tabix-indexed files. Queries use HTTP byte-range requests so
only the compressed blocks covering each query region are transferred.

Modules
-------
presets : Named URL presets for common databases (CADD, gnomAD).
cache : SQLite-backed score cache with configurable TTL.
cadd : Remote tabix CADD Phred score backend.
gnomad : Remote tabix gnomAD frequency backend.
circuit_breaker : Failure tracking with automatic open/close.
config : RemoteTabixConfig dataclass.
"""

from vartriage.remote.config import RemoteTabixConfig

__all__ = ["RemoteTabixConfig"]
