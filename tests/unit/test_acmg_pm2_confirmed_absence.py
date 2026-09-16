"""PM2 fires only on a consulted, confirmed-absent frequency result.

The criterion is "absent from population controls". That requires the
population database to have been consulted and to not contain the variant.
A variant with no frequency data at all (never queried, or annotation not
run) is not evidence of rarity and must not manufacture a PM2 tag.
"""

from __future__ import annotations

from vartriage.classification.acmg import ACMGClassifier
from vartriage.models.variant import (
    AnnotatedVariant,
    EvidenceTag,
    FunctionalConsequence,
    PopulationFrequencies,
    ScoredVariant,
    Variant,
)


def _scored(
    *,
    allele_frequency: float | None,
    frequency_unknown: bool,
    population_frequencies: PopulationFrequencies | None = None,
) -> ScoredVariant:
    variant = Variant(
        chrom="chr1",
        pos=1000,
        id=None,
        ref="A",
        alt="G",
        qual=None,
        filter_status="PASS",
        info={},
    )
    annotated = AnnotatedVariant(
        variant=variant,
        consequence=FunctionalConsequence.MISSENSE,
        gene_name="GENEX",
        allele_frequency=allele_frequency,
        frequency_unknown=frequency_unknown,
        population_frequencies=population_frequencies,
    )
    return ScoredVariant(annotated=annotated)


def _pm2_tags(variant: ScoredVariant) -> set[EvidenceTag]:
    clf = ACMGClassifier()
    tags: set[EvidenceTag] = set()
    missing: set[str] = set()
    clf._evaluate_pm2(variant, tags, missing)
    return tags


class TestPM2ConfirmedAbsence:
    def test_no_frequency_data_does_not_fire_pm2(self) -> None:
        # allele_frequency None AND frequency_unknown False = never queried.
        variant = _scored(allele_frequency=None, frequency_unknown=False)
        assert EvidenceTag.PM2 not in _pm2_tags(variant)

    def test_all_none_population_data_does_not_fire_pm2(self) -> None:
        pop = PopulationFrequencies()  # every subfield None
        variant = _scored(
            allele_frequency=None,
            frequency_unknown=False,
            population_frequencies=pop,
        )
        assert EvidenceTag.PM2 not in _pm2_tags(variant)

    def test_confirmed_absent_fires_pm2(self) -> None:
        # gnomAD consulted and returned nothing: frequency_unknown True.
        variant = _scored(allele_frequency=None, frequency_unknown=True)
        assert EvidenceTag.PM2 in _pm2_tags(variant)

    def test_rare_observed_frequency_fires_pm2(self) -> None:
        variant = _scored(allele_frequency=0.00001, frequency_unknown=False)
        assert EvidenceTag.PM2 in _pm2_tags(variant)

    def test_common_frequency_does_not_fire_pm2(self) -> None:
        variant = _scored(allele_frequency=0.02, frequency_unknown=False)
        assert EvidenceTag.PM2 not in _pm2_tags(variant)

    def test_no_data_records_gnomad_as_missing_source(self) -> None:
        variant = _scored(allele_frequency=None, frequency_unknown=False)
        clf = ACMGClassifier()
        tags: set[EvidenceTag] = set()
        missing: set[str] = set()
        clf._evaluate_pm2(variant, tags, missing)
        assert "gnomAD" in missing
