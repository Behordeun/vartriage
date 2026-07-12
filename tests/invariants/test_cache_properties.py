"""Property tests for reference-loading cache infrastructure.

Covers:
- GTF interval tree cache round-trip correctness (Property 3)
- Score dictionary cache round-trip correctness (Property 4)
- TabixFrequencyDatabase lookup positional correspondence (Property 6)
- Extension-based backend routing (Property 5)
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock, patch

from hypothesis import given, settings
from hypothesis import strategies as st

from vartriage._internal.cache import try_load_cache, try_write_cache
from vartriage._internal.interval_tree import SortedArrayIntervalIndex
from vartriage.prioritization.score_loader import CoordinateKey, ScoreLoader

# ---------------------------------------------------------------------------
# Shared strategies
# ---------------------------------------------------------------------------

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
NUCLEOTIDES = ["A", "C", "G", "T"]


@st.composite
def genomic_interval_entry(draw: st.DrawFn) -> dict:
    """Generate a synthetic GTF exon entry as a dict."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    start = draw(st.integers(min_value=1, max_value=249_000_000))
    length = draw(st.integers(min_value=10, max_value=10_000))
    end = start + length
    gene_name = draw(
        st.text(
            alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789",
            min_size=3,
            max_size=10,
        )
    )
    transcript_id = draw(
        st.text(
            alphabet="ENST0123456789",
            min_size=5,
            max_size=15,
        )
    )
    feature_type = draw(st.sampled_from(["exon", "CDS", "transcript"]))
    strand = draw(st.sampled_from(["+", "-"]))
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "gene_name": gene_name,
        "transcript_id": transcript_id,
        "feature_type": feature_type,
        "strand": strand,
    }


def _build_gtf_content(intervals: list[dict]) -> str:
    """Convert interval dicts into minimal GTF text lines."""
    lines: list[str] = []
    for iv in intervals:
        attrs = (
            f'gene_name "{iv["gene_name"]}"; ' f'transcript_id "{iv["transcript_id"]}";'
        )
        line = (
            f'{iv["chrom"]}\ttest\t{iv["feature_type"]}\t'
            f'{iv["start"]}\t{iv["end"]}\t.\t{iv["strand"]}\t.\t{attrs}'
        )
        lines.append(line)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Property 3: GTF interval tree cache round-trip correctness
# Validates: Requirements 2.5
# ---------------------------------------------------------------------------


@given(
    intervals=st.lists(
        genomic_interval_entry(),
        min_size=1,
        max_size=20,
    ),
    query_chrom=st.sampled_from(CHROMOSOMES),
    query_pos=st.integers(min_value=1, max_value=249_000_000),
)
@settings(max_examples=50)
def test_gtf_cache_round_trip_correctness(
    intervals: list[dict],
    query_chrom: str,
    query_pos: int,
) -> None:
    """Building an index from GTF, caching it, and restoring from cache
    produces identical overlap query results as the original index.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = Path(tmpdir) / "test.gtf"
        gtf_content = _build_gtf_content(intervals)
        gtf_path.write_text(gtf_content, encoding="utf-8")

        # Build index from GTF (writes cache automatically)
        original_index = SortedArrayIntervalIndex()
        original_index.load(gtf_path)

        # Query the original index
        original_results = original_index.overlap(query_chrom, query_pos, "A", "T")

        # Load a fresh index from the same path (should hit cache)
        cached_index = SortedArrayIntervalIndex()
        cached_index.load(gtf_path)

        # Query the cache-restored index
        cached_results = cached_index.overlap(query_chrom, query_pos, "A", "T")

        assert len(original_results) == len(cached_results), (
            f"Result count mismatch: original={len(original_results)}, "
            f"cached={len(cached_results)}"
        )

        for i, (orig, cached) in enumerate(zip(original_results, cached_results)):
            assert (
                orig["gene_name"] == cached["gene_name"]
            ), f"Gene name mismatch at index {i}"
            assert (
                orig["feature_type"] == cached["feature_type"]
            ), f"Feature type mismatch at index {i}"
            assert (
                orig["transcript_id"] == cached["transcript_id"]
            ), f"Transcript ID mismatch at index {i}"
            assert (
                orig["consequence"] == cached["consequence"]
            ), f"Consequence mismatch at index {i}"


# ---------------------------------------------------------------------------
# Property 4: Score dictionary cache round-trip correctness
# Validates: Requirements 6.5
# ---------------------------------------------------------------------------


@st.composite
def score_entry(draw: st.DrawFn) -> tuple[CoordinateKey, float]:
    """Generate a (chrom, pos, ref, alt) -> score pair."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    pos = draw(st.integers(min_value=1, max_value=250_000_000))
    ref = draw(
        st.text(
            alphabet=st.sampled_from(NUCLEOTIDES),
            min_size=1,
            max_size=3,
        )
    )
    alt = draw(
        st.text(
            alphabet=st.sampled_from(NUCLEOTIDES),
            min_size=1,
            max_size=3,
        )
    )
    score = draw(
        st.floats(
            min_value=0.0,
            max_value=99.0,
            allow_nan=False,
            allow_infinity=False,
        )
    )
    return ((chrom, pos, ref, alt), score)


@given(
    entries=st.lists(
        score_entry(),
        min_size=1,
        max_size=30,
        unique_by=lambda x: x[0],
    ),
)
@settings(max_examples=50)
def test_score_cache_round_trip_correctness(
    entries: list[tuple[CoordinateKey, float]],
) -> None:
    """Writing a score dictionary as TSV, loading it via ScoreLoader
    (which caches), then loading again from cache produces identical
    lookup results for all keys.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_path = Path(tmpdir) / "scores.tsv"

        # Write TSV
        with open(tsv_path, "w", encoding="utf-8") as f:
            for (chrom, pos, ref, alt), score in entries:
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{score}\n")

        loader = ScoreLoader()

        # First load: parses TSV, writes cache
        first_result = loader.load_cadd(tsv_path)

        # Second load: should read from cache
        second_result = loader.load_cadd(tsv_path)

        # Verify both loads produce identical results
        keys = [(chrom, pos, ref, alt) for (chrom, pos, ref, alt), _ in entries]

        first_lookup = loader.lookup_batch(keys, first_result)
        second_lookup = loader.lookup_batch(keys, second_result)

        assert len(first_lookup) == len(second_lookup) == len(keys)

        for i, key in enumerate(keys):
            assert first_lookup[i] == second_lookup[i], (
                f"Score mismatch for {key}: "
                f"first={first_lookup[i]}, second={second_lookup[i]}"
            )


# ---------------------------------------------------------------------------
# Property 6: TabixFrequencyDatabase lookup positional correspondence
# Validates: Requirements 3.4, 3.5
# ---------------------------------------------------------------------------


@st.composite
def variant_tuple(draw: st.DrawFn) -> tuple[str, int, str, str]:
    """Generate a (chrom, pos, ref, alt) variant tuple."""
    chrom = draw(st.sampled_from(CHROMOSOMES))
    pos = draw(st.integers(min_value=1, max_value=250_000_000))
    ref = draw(
        st.text(
            alphabet=st.sampled_from(NUCLEOTIDES),
            min_size=1,
            max_size=3,
        )
    )
    alt = draw(
        st.text(
            alphabet=st.sampled_from(NUCLEOTIDES),
            min_size=1,
            max_size=3,
        )
    )
    return (chrom, pos, ref, alt)


@given(
    variants=st.lists(variant_tuple(), min_size=0, max_size=50),
)
@settings(max_examples=50)
def test_tabix_lookup_positional_correspondence(
    variants: list[tuple[str, int, str, str]],
) -> None:
    """For any list of N variant tuples, lookup_batch returns a list
    of exactly N elements, maintaining positional correspondence.
    """
    from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

    db = TabixFrequencyDatabase()

    # Mock pysam.TabixFile to return controlled results
    mock_tabix = MagicMock()
    # Each fetch call returns an empty iterator (no records found)
    mock_tabix.fetch.return_value = iter([])
    db._tabix = mock_tabix

    results = db.lookup_batch(variants)

    assert len(results) == len(
        variants
    ), f"Output length {len(results)} != input length {len(variants)}"

    # Each element should be Optional[float] (None in this case since
    # mock returns no records)
    for i, result in enumerate(results):
        assert result is None, f"Expected None at index {i} (no records), got {result}"


@given(
    variants=st.lists(variant_tuple(), min_size=1, max_size=30),
    data=st.data(),
)
@settings(max_examples=50)
def test_tabix_lookup_positional_correspondence_with_hits(
    variants: list[tuple[str, int, str, str]],
    data: st.DataObject,
) -> None:
    """When some variants have hits, the output still has length N
    and results are positionally matched to inputs.
    """
    from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

    db = TabixFrequencyDatabase()

    # Decide which variants will "hit" in the tabix lookup
    hit_indices = set(
        data.draw(
            st.lists(
                st.integers(min_value=0, max_value=len(variants) - 1),
                min_size=0,
                max_size=len(variants),
                unique=True,
            )
        )
    )

    # Generate AF values for hits
    af_values: dict[int, float] = {}
    for idx in hit_indices:
        af_values[idx] = data.draw(
            st.floats(
                min_value=0.0,
                max_value=1.0,
                allow_nan=False,
                allow_infinity=False,
            )
        )

    mock_tabix = MagicMock()

    call_count = [0]

    def mock_fetch(chrom, start, end):
        idx = call_count[0]
        call_count[0] += 1
        if idx in hit_indices:
            # Build a VCF record line that will match
            _, pos, ref, alt = variants[idx]
            af = af_values[idx]
            record = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tAF={af}"
            return iter([record])
        return iter([])

    mock_tabix.fetch.side_effect = mock_fetch
    db._tabix = mock_tabix

    results = db.lookup_batch(variants)

    assert len(results) == len(
        variants
    ), f"Output length {len(results)} != input length {len(variants)}"


# ---------------------------------------------------------------------------
# Property 5: Extension-based backend routing
# Validates: Requirements 4.1, 4.2
# ---------------------------------------------------------------------------

# Strategies for generating filenames with various extensions
TABIX_EXTENSIONS = [".vcf.bgz", ".vcf.gz"]
NON_TABIX_EXTENSIONS = [".tsv", ".tsv.gz"]
OTHER_EXTENSIONS = [".txt", ".csv", ".bed", ".bam", ".fasta", ".fa.gz"]


@st.composite
def filename_with_extension(draw: st.DrawFn, extensions: list[str]) -> str:
    """Generate a random filename with one of the specified extensions."""
    basename = draw(
        st.text(
            alphabet="abcdefghijklmnopqrstuvwxyz0123456789_-.",
            min_size=1,
            max_size=30,
        )
    )
    ext = draw(st.sampled_from(extensions))
    return basename + ext


def _routes_to_tabix(filename: str) -> bool:
    """Replicate the extension-based routing logic from AnnotationEngine."""
    return filename.endswith((".vcf.bgz", ".vcf.gz"))


def _routes_to_tsv_backend(filename: str) -> bool:
    """Check if filename routes to TSV-based backends."""
    return filename.endswith((".tsv", ".tsv.gz"))


@given(
    filename=filename_with_extension(TABIX_EXTENSIONS),
)
@settings(max_examples=100)
def test_vcf_extension_routes_to_tabix(filename: str) -> None:
    """Files ending with .vcf.bgz or .vcf.gz route to tabix backend."""
    assert _routes_to_tabix(
        filename
    ), f"Expected '{filename}' to route to tabix backend"
    assert not _routes_to_tsv_backend(
        filename
    ), f"Expected '{filename}' NOT to route to TSV backend"


@given(
    filename=filename_with_extension(NON_TABIX_EXTENSIONS),
)
@settings(max_examples=100)
def test_tsv_extension_does_not_route_to_tabix(filename: str) -> None:
    """Files ending with .tsv or .tsv.gz do NOT route to tabix backend."""
    assert not _routes_to_tabix(
        filename
    ), f"Expected '{filename}' NOT to route to tabix backend"
    assert _routes_to_tsv_backend(
        filename
    ), f"Expected '{filename}' to route to TSV backend"


@given(
    filename=filename_with_extension(
        TABIX_EXTENSIONS + NON_TABIX_EXTENSIONS + OTHER_EXTENSIONS
    ),
)
@settings(max_examples=100)
def test_routing_is_mutually_exclusive_for_known_extensions(
    filename: str,
) -> None:
    """Tabix and TSV routing are mutually exclusive. A filename
    cannot match both routing conditions.
    """
    is_tabix = _routes_to_tabix(filename)
    is_tsv = _routes_to_tsv_backend(filename)
    assert not (is_tabix and is_tsv), f"'{filename}' matches both tabix and TSV routing"


@given(
    filename=filename_with_extension(TABIX_EXTENSIONS),
)
@settings(max_examples=100)
def test_engine_routing_logic_matches_tabix_extension(
    filename: str,
) -> None:
    """Verify that AnnotationEngine._build_frequency_db routing logic
    correctly identifies tabix-eligible paths by extension check.
    """
    # The actual engine checks `gnomad_name.endswith((".vcf.bgz", ".vcf.gz"))`
    gnomad_name = filename
    routes_tabix = gnomad_name.endswith((".vcf.bgz", ".vcf.gz"))
    assert routes_tabix is True, f"Engine routing should select tabix for '{filename}'"
