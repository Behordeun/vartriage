"""Numpy and polars vectorized helper utilities.

Provides numpy-based batch score normalization helpers (always available)
and optional polars-based batch join helpers for coordinate overlaps and
frequency lookups (guarded by import check).

All operations are designed to stay within 2GB memory overhead by
processing in bounded chunks when needed.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray

# Memory limit for vectorized operations (2 GB)
_MAX_MEMORY_BYTES: int = 2 * 1024 * 1024 * 1024

# Approximate bytes per variant row in numpy arrays
# (4 float64 columns = 32 bytes, plus overhead)
_BYTES_PER_VARIANT_ROW: int = 64

# Max variants per chunk to stay within memory budget
_MAX_CHUNK_SIZE: int = _MAX_MEMORY_BYTES // _BYTES_PER_VARIANT_ROW

# Fallback chunk size on MemoryError
_FALLBACK_CHUNK_SIZE: int = 500_000


# --------------------------------------------------------------------------
# Numpy-based score normalization helpers (hard dependency, always available)
# --------------------------------------------------------------------------


def normalize_cadd_phred_vectorized(
    scores: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Normalize CADD Phred scores to 0.0-1.0 scale via vectorized division.

    Applies the formula min(score / 99.0, 1.0) element-wise using numpy
    broadcasting. NaN values in the input remain NaN in the output.

    Parameters
    ----------
    scores : NDArray[np.float64]
        Array of raw CADD Phred scores. NaN represents missing data,
        negative values are treated as invalid and set to NaN.

    Returns
    -------
    NDArray[np.float64]
        Normalized scores in [0.0, 1.0] range. Invalid (negative) entries
        become NaN.
    """
    result: NDArray[np.float64] = np.array(scores, dtype=np.float64, copy=True)

    # Mark negative scores as invalid
    invalid_mask = result < 0.0
    result[invalid_mask] = np.nan

    # Normalize valid scores: min(score / 99.0, 1.0)
    valid_mask = ~np.isnan(result)
    result[valid_mask] = np.minimum(result[valid_mask] / 99.0, 1.0)

    return result


def validate_revel_vectorized(
    scores: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Validate REVEL scores, setting out-of-range values to NaN.

    REVEL scores must be in [0.0, 1.0]. Values outside this range are
    marked as NaN (invalid).

    Parameters
    ----------
    scores : NDArray[np.float64]
        Array of REVEL scores. NaN represents missing data.

    Returns
    -------
    NDArray[np.float64]
        Validated scores with out-of-range entries set to NaN.
    """
    result: NDArray[np.float64] = np.array(scores, dtype=np.float64, copy=True)

    # Mask out-of-range values (ignoring existing NaN entries)
    valid_mask = ~np.isnan(result)
    out_of_range = valid_mask & ((result < 0.0) | (result > 1.0))
    result[out_of_range] = np.nan

    return result


def compute_composite_vectorized(
    cadd_normalized: NDArray[np.float64],
    revel_scores: NDArray[np.float64],
    revel_weight: float = 0.6,
    cadd_weight: float = 0.4,
) -> NDArray[np.float64]:
    """Compute composite pathogenicity rank from normalized scores.

    Formula when both available: (REVEL * revel_weight) + (CADD * cadd_weight).
    When only one score is available, uses that score directly.
    When neither is available, the result is NaN.

    Parameters
    ----------
    cadd_normalized : NDArray[np.float64]
        Normalized CADD scores (0.0-1.0). NaN for missing.
    revel_scores : NDArray[np.float64]
        Validated REVEL scores (0.0-1.0). NaN for missing.
    revel_weight : float
        Weight for REVEL in the composite formula. Default 0.6.
    cadd_weight : float
        Weight for CADD in the composite formula. Default 0.4.

    Returns
    -------
    NDArray[np.float64]
        Composite scores. NaN where neither score is available.

    Raises
    ------
    ValueError
        If input arrays have different shapes.
    """
    if cadd_normalized.shape != revel_scores.shape:
        raise ValueError(
            f"Array shapes must match: CADD {cadd_normalized.shape} "
            f"vs REVEL {revel_scores.shape}"
        )

    has_cadd = ~np.isnan(cadd_normalized)
    has_revel = ~np.isnan(revel_scores)

    both = has_cadd & has_revel
    cadd_only = has_cadd & ~has_revel
    revel_only = ~has_cadd & has_revel

    composite = np.full_like(cadd_normalized, np.nan)
    composite[both] = (
        revel_scores[both] * revel_weight + cadd_normalized[both] * cadd_weight
    )
    composite[cadd_only] = cadd_normalized[cadd_only]
    composite[revel_only] = revel_scores[revel_only]

    return composite


def batch_normalize_scores(
    cadd_scores: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
    revel_weight: float = 0.6,
    cadd_weight: float = 0.4,
) -> tuple[list[Optional[float]], list[Optional[float]], list[Optional[float]]]:
    """Batch normalize and compute composite scores using numpy vectorization.

    Converts Python sequences to numpy arrays, applies vectorized
    normalization, and returns results as Python lists. Processes in
    chunks if the dataset would exceed memory budget.

    Parameters
    ----------
    cadd_scores : Sequence[Optional[float]]
        Raw CADD Phred scores. None for missing data.
    revel_scores : Sequence[Optional[float]]
        Raw REVEL scores. None for missing data.
    revel_weight : float
        Weight for REVEL in composite formula. Default 0.6.
    cadd_weight : float
        Weight for CADD in composite formula. Default 0.4.

    Returns
    -------
    tuple[list[Optional[float]], list[Optional[float]], list[Optional[float]]]
        Three lists of equal length:
        - Normalized CADD scores (0.0-1.0 or None)
        - Validated REVEL scores (0.0-1.0 or None)
        - Composite ranks (or None when neither score available)

    Raises
    ------
    ValueError
        If input sequences have different lengths.
    """
    n = len(cadd_scores)
    if n != len(revel_scores):
        raise ValueError(
            f"Score sequences must have equal length: "
            f"CADD has {n}, REVEL has {len(revel_scores)}"
        )

    if n == 0:
        return [], [], []

    # Determine chunk size based on memory budget
    chunk_size = min(n, _MAX_CHUNK_SIZE)

    cadd_out: list[Optional[float]] = []
    revel_out: list[Optional[float]] = []
    composite_out: list[Optional[float]] = []

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)

        try:
            cadd_chunk, revel_chunk, comp_chunk = _process_score_chunk(
                cadd_scores[start:end],
                revel_scores[start:end],
                revel_weight,
                cadd_weight,
            )
        except MemoryError:
            # Fall back to smaller chunks on memory pressure
            cadd_chunk, revel_chunk, comp_chunk = _process_with_fallback(
                cadd_scores[start:end],
                revel_scores[start:end],
                revel_weight,
                cadd_weight,
            )

        cadd_out.extend(cadd_chunk)
        revel_out.extend(revel_chunk)
        composite_out.extend(comp_chunk)

    return cadd_out, revel_out, composite_out


def _process_score_chunk(
    cadd_scores: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
    revel_weight: float,
    cadd_weight: float,
) -> tuple[list[Optional[float]], list[Optional[float]], list[Optional[float]]]:
    """Process a single chunk of scores through vectorized normalization."""
    cadd_arr = np.array(
        [s if s is not None else np.nan for s in cadd_scores],
        dtype=np.float64,
    )
    revel_arr = np.array(
        [s if s is not None else np.nan for s in revel_scores],
        dtype=np.float64,
    )

    cadd_norm = normalize_cadd_phred_vectorized(cadd_arr)
    revel_valid = validate_revel_vectorized(revel_arr)
    composite = compute_composite_vectorized(
        cadd_norm, revel_valid, revel_weight, cadd_weight
    )

    return (
        _ndarray_to_optional_list(cadd_norm),
        _ndarray_to_optional_list(revel_valid),
        _ndarray_to_optional_list(composite),
    )


def _process_with_fallback(
    cadd_scores: Sequence[Optional[float]],
    revel_scores: Sequence[Optional[float]],
    revel_weight: float,
    cadd_weight: float,
) -> tuple[list[Optional[float]], list[Optional[float]], list[Optional[float]]]:
    """Process scores in smaller fallback chunks after MemoryError."""
    n = len(cadd_scores)
    chunk_size = min(n, _FALLBACK_CHUNK_SIZE)

    cadd_out: list[Optional[float]] = []
    revel_out: list[Optional[float]] = []
    composite_out: list[Optional[float]] = []

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        cadd_chunk, revel_chunk, comp_chunk = _process_score_chunk(
            cadd_scores[start:end],
            revel_scores[start:end],
            revel_weight,
            cadd_weight,
        )
        cadd_out.extend(cadd_chunk)
        revel_out.extend(revel_chunk)
        composite_out.extend(comp_chunk)

    return cadd_out, revel_out, composite_out


def _ndarray_to_optional_list(
    arr: NDArray[np.float64],
) -> list[Optional[float]]:
    """Convert a numpy array to a list of Optional[float], NaN -> None."""
    mask = np.isnan(arr)
    result: list[Optional[float]] = []
    for i in range(len(arr)):
        if mask[i]:
            result.append(None)
        else:
            result.append(float(arr[i]))
    return result


# --------------------------------------------------------------------------
# Polars-based batch join helpers (optional, guarded by import check)
# --------------------------------------------------------------------------

try:
    import polars as pl

    _POLARS_AVAILABLE = True
except ImportError:
    _POLARS_AVAILABLE = False


def polars_available() -> bool:
    """Check if polars is available for accelerated batch operations.

    Returns
    -------
    bool
        True if polars is importable.
    """
    return _POLARS_AVAILABLE


def batch_frequency_join(
    variants: list[tuple[str, int, str, str]],
    reference_path: str,
) -> list[Optional[float]]:
    """Perform a batch left join for frequency lookups using polars.

    Reads the reference TSV into a LazyFrame and joins against the
    input variant coordinates. Returns allele frequencies in the same
    order as the input, with None for variants not found.

    Parameters
    ----------
    variants : list[tuple[str, int, str, str]]
        List of (chrom, pos, ref, alt) tuples to look up.
    reference_path : str
        Path to the tab-separated reference file with columns:
        chrom, pos, ref, alt, af.

    Returns
    -------
    list[Optional[float]]
        Allele frequencies aligned with input order. None for
        variants absent from the reference.

    Raises
    ------
    ImportError
        If polars is not installed.
    RuntimeError
        If the join operation fails.
    """
    if not _POLARS_AVAILABLE:
        raise ImportError(
            "polars is required for batch_frequency_join. "
            "Install with: pip install vartriage[accelerated]"
        )

    if not variants:
        return []

    try:
        ref_df = pl.read_csv(
            reference_path,
            separator="\t",
            has_header=True,
            schema_overrides={
                "chrom": pl.Utf8,
                "pos": pl.Int64,
                "ref": pl.Utf8,
                "alt": pl.Utf8,
                "af": pl.Float64,
            },
        )
    except Exception as exc:
        raise RuntimeError(f"Failed to read reference file: {exc}") from exc

    # Normalize column names
    ref_df = ref_df.rename({col: col.lower() for col in ref_df.columns})

    query_df = pl.DataFrame(
        {
            "chrom": [v[0] for v in variants],
            "pos": [v[1] for v in variants],
            "ref": [v[2] for v in variants],
            "alt": [v[3] for v in variants],
            "_row_idx": list(range(len(variants))),
        }
    )

    result = (
        query_df.lazy()
        .join(
            ref_df.lazy(),
            on=["chrom", "pos", "ref", "alt"],
            how="left",
        )
        .sort("_row_idx")
        .collect()
    )

    af_series = result["af"]
    return [float(val) if val is not None else None for val in af_series.to_list()]


def batch_coordinate_overlap_join(
    variants: list[tuple[str, int, int]],
    regions: list[tuple[str, int, int, str]],
) -> list[list[str]]:
    """Perform a batch interval overlap join using polars.

    For each variant (chrom, start, end), finds all regions
    (chrom, start, end, label) that overlap it.

    Parameters
    ----------
    variants : list[tuple[str, int, int]]
        List of (chrom, start, end) tuples representing variant spans.
    regions : list[tuple[str, int, int, str]]
        List of (chrom, start, end, label) tuples representing
        genomic regions (e.g., exons, CDS features).

    Returns
    -------
    list[list[str]]
        For each variant, a list of labels from overlapping regions.
        Empty list if no overlaps found for that variant.

    Raises
    ------
    ImportError
        If polars is not installed.
    """
    if not _POLARS_AVAILABLE:
        raise ImportError(
            "polars is required for batch_coordinate_overlap_join. "
            "Install with: pip install vartriage[accelerated]"
        )

    if not variants or not regions:
        return [[] for _ in variants]

    # Build DataFrames
    var_df = pl.DataFrame(
        {
            "var_chrom": [v[0] for v in variants],
            "var_start": [v[1] for v in variants],
            "var_end": [v[2] for v in variants],
            "_var_idx": list(range(len(variants))),
        }
    )

    reg_df = pl.DataFrame(
        {
            "reg_chrom": [r[0] for r in regions],
            "reg_start": [r[1] for r in regions],
            "reg_end": [r[2] for r in regions],
            "label": [r[3] for r in regions],
        }
    )

    # Cross join on chromosome, then filter for overlapping intervals
    # Overlap condition: var_start < reg_end AND var_end > reg_start
    joined = (
        var_df.lazy()
        .join(reg_df.lazy(), left_on="var_chrom", right_on="reg_chrom", how="inner")
        .filter(
            (pl.col("var_start") < pl.col("reg_end"))
            & (pl.col("var_end") > pl.col("reg_start"))
        )
        .select(["_var_idx", "label"])
        .collect()
    )

    # Group by variant index and collect labels
    results: list[list[str]] = [[] for _ in variants]

    if joined.height > 0:
        grouped = joined.group_by("_var_idx").agg(pl.col("label"))
        for row in grouped.iter_rows():
            idx = row[0]
            labels = row[1]
            results[idx] = labels

    return results
