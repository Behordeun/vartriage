"""Hypothesis strategies for generating genomic variant data.

Strategies produce realistic genomic data: valid chromosome names,
biologically plausible positions, nucleotide alleles, and QUAL scores
consistent with real sequencing pipelines.
"""

from __future__ import annotations

from hypothesis import strategies as st
from hypothesis.strategies import SearchStrategy

from vartriage.models.variant import (
    AnnotatedVariant,
    ClassifiedVariant,
    EvidenceTag,
    FunctionalConsequence,
    ClinVarAssertion,
    ACMGClassification,
    ScoredVariant,
    Variant,
)


CHROMOSOMES: list[str] = [
    f"chr{i}" for i in range(1, 23)
] + ["chrX", "chrY", "chrM"]

NUCLEOTIDES: list[str] = ["A", "C", "G", "T"]

FILTER_PASS_VALUES: list[str] = ["PASS", "."]
FILTER_FAIL_VALUES: list[str] = ["LowQual", "q10", "FAIL", "IndelQual"]


@st.composite
def chromosome(draw: st.DrawFn) -> str:
    """Generate a valid human chromosome name (chr1-chr22, chrX, chrY, chrM)."""
    return draw(st.sampled_from(CHROMOSOMES))


@st.composite
def genomic_position(draw: st.DrawFn) -> int:
    """Generate a 1-based genomic position within a realistic range.

    Human chromosomes range from ~50 Mb (chr21) to ~249 Mb (chr1).
    We use a conservative upper bound of 250,000,000.
    """
    return draw(st.integers(min_value=1, max_value=250_000_000))


@st.composite
def allele(draw: st.DrawFn, min_length: int = 1, max_length: int = 50) -> str:
    """Generate a nucleotide allele string (A, C, G, T sequences).

    Parameters
    ----------
    min_length : int
        Minimum allele length (1 for SNVs).
    max_length : int
        Maximum allele length (indels can be longer).
    """
    length = draw(st.integers(min_value=min_length, max_value=max_length))
    bases = draw(
        st.lists(
            st.sampled_from(NUCLEOTIDES),
            min_size=length,
            max_size=length,
        )
    )
    return "".join(bases)


@st.composite
def snv_allele(draw: st.DrawFn) -> str:
    """Generate a single nucleotide variant allele."""
    return draw(st.sampled_from(NUCLEOTIDES))


@st.composite
def qual_score(draw: st.DrawFn) -> float:
    """Generate a Phred-scaled quality score.

    Realistic range is 0 to ~10000 for most variant callers,
    though the spec allows up to 1_000_000.
    """
    return draw(st.floats(min_value=0.0, max_value=10000.0, allow_nan=False))


@st.composite
def valid_variant(draw: st.DrawFn) -> Variant:
    """Generate a valid Variant record with realistic genomic data.

    Produces variants that would pass basic VCF parsing: valid chromosome,
    position, nucleotide alleles, a non-None QUAL score, and a passing
    FILTER status.
    """
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(allele(min_length=1, max_length=10))
    alt = draw(allele(min_length=1, max_length=10))
    qual = draw(qual_score())
    filter_status = draw(st.sampled_from(FILTER_PASS_VALUES))
    variant_id = draw(
        st.one_of(st.none(), st.text(min_size=1, max_size=20, alphabet="rs0123456789"))
    )

    return Variant(
        chrom=chrom,
        pos=pos,
        id=variant_id,
        ref=ref,
        alt=alt,
        qual=qual,
        filter_status=filter_status,
        info={},
    )


@st.composite
def variant_with_filter(
    draw: st.DrawFn,
    filter_values: list[str] | None = None,
    qual_range: tuple[float | None, float | None] = (0.0, 10000.0),
    allow_missing_qual: bool = False,
) -> Variant:
    """Generate a Variant with configurable FILTER and QUAL properties.

    Useful for testing quality filtering logic with specific filter/qual
    combinations.
    """
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())

    if filter_values is None:
        filter_values = FILTER_PASS_VALUES + FILTER_FAIL_VALUES
    filter_status = draw(st.sampled_from(filter_values))

    if allow_missing_qual:
        qual: float | None = draw(
            st.one_of(
                st.none(),
                st.floats(
                    min_value=qual_range[0] or 0.0,
                    max_value=qual_range[1] or 10000.0,
                    allow_nan=False,
                ),
            )
        )
    else:
        qual = draw(
            st.floats(
                min_value=qual_range[0] or 0.0,
                max_value=qual_range[1] or 10000.0,
                allow_nan=False,
            )
        )

    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=qual,
        filter_status=filter_status,
        info={},
    )


@st.composite
def scored_variant(
    draw: st.DrawFn,
    with_cadd: bool | None = None,
    with_revel: bool | None = None,
) -> ScoredVariant:
    """Generate a ScoredVariant with realistic pathogenicity scores.

    Parameters
    ----------
    with_cadd : bool, optional
        If True, always include CADD score. If False, never include.
        If None, randomly decide.
    with_revel : bool, optional
        If True, always include REVEL score. If False, never include.
        If None, randomly decide.
    """
    variant = draw(valid_variant())
    consequence = draw(st.sampled_from(list(FunctionalConsequence)))
    allele_frequency = draw(
        st.one_of(
            st.none(),
            st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
        )
    )
    clinvar_assertion = draw(
        st.one_of(st.none(), st.sampled_from(list(ClinVarAssertion)))
    )
    frequency_unknown = allele_frequency is None

    annotated = AnnotatedVariant(
        variant=variant,
        consequence=consequence,
        allele_frequency=allele_frequency,
        clinvar_assertion=clinvar_assertion,
        frequency_unknown=frequency_unknown,
        clinvar_unknown=clinvar_assertion is None,
    )

    include_cadd = draw(st.booleans()) if with_cadd is None else with_cadd
    include_revel = draw(st.booleans()) if with_revel is None else with_revel

    cadd_phred: float | None = None
    cadd_normalized: float | None = None
    revel_score: float | None = None
    composite_rank: float | None = None

    if include_cadd:
        cadd_phred = draw(st.floats(min_value=0.0, max_value=60.0, allow_nan=False))
        cadd_normalized = min(cadd_phred / 99.0, 1.0)

    if include_revel:
        revel_score = draw(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))

    if revel_score is not None and cadd_normalized is not None:
        composite_rank = (revel_score * 0.6) + (cadd_normalized * 0.4)
    elif revel_score is not None:
        composite_rank = revel_score
    elif cadd_normalized is not None:
        composite_rank = cadd_normalized

    return ScoredVariant(
        annotated=annotated,
        cadd_phred=cadd_phred,
        cadd_normalized=cadd_normalized,
        revel_score=revel_score,
        composite_rank=composite_rank,
    )


@st.composite
def evidence_tag_set(draw: st.DrawFn) -> frozenset[EvidenceTag]:
    """Generate a valid subset of ACMG evidence tags.

    Returns a frozenset containing zero or more evidence tags,
    representing a realistic combination that could be assigned
    to a single variant.
    """
    tags = draw(
        st.frozensets(st.sampled_from(list(EvidenceTag)), min_size=0, max_size=4)
    )
    return tags


def evidence_tag_strategy() -> SearchStrategy[frozenset[EvidenceTag]]:
    """Return a strategy for evidence tag frozensets (non-composite version)."""
    return st.frozensets(st.sampled_from(list(EvidenceTag)), min_size=0, max_size=4)
