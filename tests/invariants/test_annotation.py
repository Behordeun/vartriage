"""Hypothesis tests for annotation engine (consequence and lookups).

Verifies functional consequence severity assignment and reference database lookup consistency.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

from hypothesis import assume, given, settings
from hypothesis import strategies as st

from tests.generators.variants import chromosome, genomic_position, snv_allele
from vartriage.annotation.clinvar import DictClinVarDatabase
from vartriage.annotation.consequence import _most_severe_consequence
from vartriage.annotation.frequency import DictFrequencyDatabase
from vartriage.models.variant import (CONSEQUENCE_SEVERITY_ORDER,
                                      ClinVarAssertion, FunctionalConsequence,
                                      Variant)

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

CONSEQUENCE_VALUES = list(FunctionalConsequence)


@st.composite
def consequence_subset(draw: st.DrawFn) -> list[FunctionalConsequence]:
    """Generate a non-empty subset of FunctionalConsequence values.

    Simulates multiple transcripts returning different consequences
    for a single variant.
    """
    consequences = draw(
        st.lists(
            st.sampled_from(CONSEQUENCE_VALUES),
            min_size=1,
            max_size=len(CONSEQUENCE_VALUES),
        )
    )
    return consequences


@st.composite
def variant_key(draw: st.DrawFn) -> tuple[str, int, str, str]:
    """Generate a (chrom, pos, ref, alt) lookup key."""
    chrom = draw(chromosome())
    pos = draw(genomic_position())
    ref = draw(snv_allele())
    alt = draw(snv_allele())
    return (chrom, pos, ref, alt)


@st.composite
def frequency_database_scenario(
    draw: st.DrawFn,
) -> tuple[dict[tuple[str, int, str, str], float], list[tuple[str, int, str, str]]]:
    """Generate a frequency database (as a dict) and query keys.

    Some query keys are in the database, some are not, testing both
    the match and no-match paths.
    """
    # Generate entries that exist in the database
    num_entries = draw(st.integers(min_value=1, max_value=20))
    db_entries: dict[tuple[str, int, str, str], float] = {}

    for _ in range(num_entries):
        key = draw(variant_key())
        freq = draw(
            st.floats(
                min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False
            )
        )
        db_entries[key] = freq

    # Generate query keys: mix of present and absent
    present_keys = draw(
        st.lists(
            st.sampled_from(list(db_entries.keys())),
            min_size=1,
            max_size=min(5, num_entries),
        )
    )
    absent_keys = draw(st.lists(variant_key(), min_size=1, max_size=5))

    # Filter out absent keys that happen to be present
    truly_absent = [k for k in absent_keys if k not in db_entries]

    query_keys = present_keys + truly_absent
    draw(st.randoms()).shuffle(query_keys)

    return db_entries, query_keys


@st.composite
def clinvar_database_scenario(
    draw: st.DrawFn,
) -> tuple[
    dict[tuple[str, int, str, str], ClinVarAssertion], list[tuple[str, int, str, str]]
]:
    """Generate a ClinVar database (as a dict) and query keys.

    Some keys match, some don't. Verifying both paths.
    """
    num_entries = draw(st.integers(min_value=1, max_value=20))
    db_entries: dict[tuple[str, int, str, str], ClinVarAssertion] = {}

    for _ in range(num_entries):
        key = draw(variant_key())
        assertion = draw(st.sampled_from(list(ClinVarAssertion)))
        db_entries[key] = assertion

    present_keys = draw(
        st.lists(
            st.sampled_from(list(db_entries.keys())),
            min_size=1,
            max_size=min(5, num_entries),
        )
    )
    absent_keys = draw(st.lists(variant_key(), min_size=1, max_size=5))

    truly_absent = [k for k in absent_keys if k not in db_entries]

    query_keys = present_keys + truly_absent
    draw(st.randoms()).shuffle(query_keys)

    return db_entries, query_keys


# ---------------------------------------------------------------------------


@given(consequences=consequence_subset())
@settings(max_examples=200)
def test_most_severe_consequence_selected(
    consequences: list[FunctionalConsequence],
) -> None:
    """Assigns single most severe consequence from overlapping transcripts.

    Given a list of consequences (simulating multiple transcript overlaps),
    the result is always the one with the lowest index in the severity
    ranking. Intergenic is assigned when no overlaps exist (empty list
    case tested separately).
    """
    # Build overlap dicts as the consequence annotator would see them
    overlaps = [{"consequence": c.value} for c in consequences]

    result = _most_severe_consequence(overlaps)

    # Determine expected: the consequence with the lowest index in severity order
    severity_rank = {c: idx for idx, c in enumerate(CONSEQUENCE_SEVERITY_ORDER)}
    expected = min(consequences, key=lambda c: severity_rank[c])

    assert result == expected, (
        f"Expected most severe={expected.value} from {[c.value for c in consequences]}, "
        f"got {result.value}"
    )


@given(data=st.data())
@settings(max_examples=200)
def test_intergenic_assigned_when_no_overlaps(data: st.DataObject) -> None:
    """Intergenic consequence is assigned when no transcript overlaps exist.

    When the overlap list is empty, the consequence must be INTERGENIC.
    """
    # Empty overlap list means no transcripts overlap this variant
    overlaps: list[dict] = []
    result = _most_severe_consequence(overlaps)

    # With empty overlaps, the function returns the initial best which is INTERGENIC
    assert (
        result == FunctionalConsequence.INTERGENIC
    ), f"Expected INTERGENIC for empty overlap list, got {result.value}"


@given(consequence=st.sampled_from(CONSEQUENCE_VALUES))
@settings(max_examples=200)
def test_single_consequence_returns_itself(
    consequence: FunctionalConsequence,
) -> None:
    """A single overlapping transcript returns its own consequence unchanged."""
    overlaps = [{"consequence": consequence.value}]
    result = _most_severe_consequence(overlaps)

    assert (
        result == consequence
    ), f"Single consequence {consequence.value} should return itself, got {result.value}"


@given(
    less_severe=st.sampled_from(CONSEQUENCE_VALUES[1:]),
    more_severe_idx=st.data(),
)
@settings(max_examples=200)
def test_severity_ranking_is_transitive(
    less_severe: FunctionalConsequence,
    more_severe_idx: st.DataObject,
) -> None:
    """Any consequence ranked higher always wins over a lower-ranked one."""
    severity_rank = {c: idx for idx, c in enumerate(CONSEQUENCE_SEVERITY_ORDER)}
    less_severe_rank = severity_rank[less_severe]

    # Pick any consequence that is more severe (lower rank index)
    assume(less_severe_rank > 0)
    more_severe = more_severe_idx.draw(
        st.sampled_from(CONSEQUENCE_SEVERITY_ORDER[:less_severe_rank])
    )

    overlaps = [
        {"consequence": less_severe.value},
        {"consequence": more_severe.value},
    ]
    result = _most_severe_consequence(overlaps)

    assert result == more_severe, (
        f"Expected {more_severe.value} to win over {less_severe.value}, "
        f"got {result.value}"
    )


# ---------------------------------------------------------------------------


@given(scenario=frequency_database_scenario())
@settings(max_examples=100)
def test_frequency_lookup_match_returns_value(
    scenario: tuple[
        dict[tuple[str, int, str, str], float],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """When a variant matches in gnomAD, the allele frequency value is attached.
    When no match exists, null is returned and a warning is emitted.
    """
    db_entries, query_keys = scenario

    # Write the database to a temp TSV file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("chrom\tpos\tref\talt\taf\n")
        for (chrom, pos, ref, alt), freq in db_entries.items():
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq}\n")
        tmp_path = Path(f.name)

    try:
        db = DictFrequencyDatabase()
        db.load(tmp_path)

        results = db.lookup_batch(query_keys)

        assert len(results) == len(query_keys)

        for key, result in zip(query_keys, results):
            if key in db_entries:
                # Match: value should be attached
                assert result is not None, f"Expected frequency for {key}, got None"
                assert abs(result - db_entries[key]) < 1e-7, (
                    f"Expected frequency {db_entries[key]} for {key}, " f"got {result}"
                )
            else:
                # No match: null returned
                assert (
                    result is None
                ), f"Expected None for absent key {key}, got {result}"
    finally:
        tmp_path.unlink(missing_ok=True)


@given(scenario=frequency_database_scenario())
@settings(max_examples=100)
def test_frequency_lookup_missing_emits_warning(
    scenario: tuple[
        dict[tuple[str, int, str, str], float],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """When no match exists in gnomAD, a MissingDataWarning is emitted.

    The warning includes the variant coordinates and the source name.
    """
    db_entries, query_keys = scenario

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("chrom\tpos\tref\talt\taf\n")
        for (chrom, pos, ref, alt), freq in db_entries.items():
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq}\n")
        tmp_path = Path(f.name)

    try:
        db = DictFrequencyDatabase()
        db.load(tmp_path)
        db.warnings.clear()

        db.lookup_batch(query_keys)

        absent_keys = [k for k in query_keys if k not in db_entries]
        assert len(db.warnings) == len(
            absent_keys
        ), f"Expected {len(absent_keys)} warnings, got {len(db.warnings)}"

        for warning, key in zip(db.warnings, absent_keys):
            assert warning.chrom == key[0]
            assert warning.pos == key[1]
            assert warning.ref == key[2]
            assert warning.alt == key[3]
            assert warning.source == "gnomAD"
    finally:
        tmp_path.unlink(missing_ok=True)


@given(scenario=clinvar_database_scenario())
@settings(max_examples=100)
def test_clinvar_lookup_match_returns_assertion(
    scenario: tuple[
        dict[tuple[str, int, str, str], ClinVarAssertion],
        list[tuple[str, int, str, str]],
    ],
) -> None:
    """When a variant matches in ClinVar, the assertion value is attached.
    When no match exists, null is returned.
    """
    db_entries, query_keys = scenario

    # Map ClinVarAssertion enum values to the TSV significance strings
    assertion_to_str = {
        ClinVarAssertion.PATHOGENIC: "Pathogenic",
        ClinVarAssertion.LIKELY_PATHOGENIC: "Likely pathogenic",
        ClinVarAssertion.VUS: "Uncertain significance",
        ClinVarAssertion.LIKELY_BENIGN: "Likely benign",
        ClinVarAssertion.BENIGN: "Benign",
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("chrom\tpos\tref\talt\tclinical_significance\n")
        for (chrom, pos, ref, alt), assertion in db_entries.items():
            sig_str = assertion_to_str[assertion]
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{sig_str}\n")
        tmp_path = Path(f.name)

    try:
        db = DictClinVarDatabase()
        db.load(tmp_path)

        results = db.lookup_batch(query_keys)

        assert len(results) == len(query_keys)

        for key, result in zip(query_keys, results):
            if key in db_entries:
                assert (
                    result is not None
                ), f"Expected ClinVar assertion for {key}, got None"
                assert (
                    result == db_entries[key]
                ), f"Expected {db_entries[key]} for {key}, got {result}"
            else:
                assert (
                    result is None
                ), f"Expected None for absent key {key}, got {result}"
    finally:
        tmp_path.unlink(missing_ok=True)
