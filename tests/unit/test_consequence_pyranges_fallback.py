"""Regression tests for the pyranges consequence backend's runtime safety.

Guarantees that a failure in the vectorized overlap path degrades to the
pure-Python annotator for the batch rather than aborting the classification
run, and that the two backends do not collide on a shared reference-cache
file (they serialize incompatible payloads).
"""

from pathlib import Path

import pytest

from vartriage.models.variant import FunctionalConsequence, Variant

pytest.importorskip("pyranges")

from vartriage._internal.cache import cache_path_for  # noqa: E402
from vartriage.annotation.consequence_pyranges import (  # noqa: E402
    PyRangesConsequenceAnnotator,
)

SAMPLE_GTF = (
    'chr1\ttest\texon\t1000\t1200\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";\n'
    'chr1\ttest\tCDS\t1050\t1190\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";\n'
)


def _variant(pos: int) -> Variant:
    return Variant(
        chrom="chr1",
        pos=pos,
        id=None,
        ref="A",
        alt="T",
        qual=30.0,
        filter_status="PASS",
    )


def _gtf(tmp_path: Path) -> Path:
    p = tmp_path / "genes.gtf"
    p.write_text(SAMPLE_GTF, encoding="utf-8")
    return p


def test_vectorized_failure_falls_back_to_pure_python(tmp_path: Path) -> None:
    """When the vectorized join raises at run time, the batch is served by the
    pure-Python backend instead of propagating the exception.
    """
    ann = PyRangesConsequenceAnnotator(_gtf(tmp_path))
    variants = [_variant(1100)]

    expected = ann.assign_batch(variants)

    def _boom(_variants):
        raise RuntimeError("simulated backend failure")

    ann._assign_batch_vectorized = _boom  # type: ignore[method-assign]
    result = ann.assign_batch(variants)

    assert result == expected
    assert len(result) == len(variants)
    assert all(isinstance(c, FunctionalConsequence) for c in result)


def test_pyranges_and_intervaltree_use_distinct_cache_files(tmp_path: Path) -> None:
    """The two backends serialize incompatible payloads, so they must not share
    a cache path. The tag namespaces them.
    """
    gtf = _gtf(tmp_path)
    pyr = cache_path_for(gtf, tag="pyranges")
    tree = cache_path_for(gtf, tag="intervaltree")
    untagged = cache_path_for(gtf)
    assert pyr != tree
    assert pyr != untagged
    assert tree != untagged
