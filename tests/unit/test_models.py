"""Unit tests for core data models (variant, cohort, config)."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.models.cohort import (
    CohortConfig,
    CohortSummary,
    CohortVariant,
    GeneBurden,
    SampleOccurrence,
)
from vartriage.models.config import (
    AnnotationConfig,
    ClinicalReportConfig,
    InheritanceConfig,
    PipelineConfig,
    PrioritizationConfig,
    QualityFilterConfig,
    ReportConfig,
    SampleConfig,
)
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    ClinVarAssertion,
    EvidenceStrength,
    EvidenceTag,
    FunctionalConsequence,
    PopulationFrequencies,
    ScoredVariant,
    Variant,
    VariantQualityMetrics,
    Zygosity,
)


# --------------------------------------------------------------------------
# models/variant.py
# --------------------------------------------------------------------------


class TestFunctionalConsequence:
    def test_all_values_are_strings(self) -> None:
        for member in FunctionalConsequence:
            assert isinstance(member.value, str)


class TestPopulationFrequencies:
    def test_max_population_af_returns_highest(self) -> None:
        pf = PopulationFrequencies(
            global_af=0.01,
            afr=0.02,
            amr=0.005,
            eas=None,
            nfe=0.008,
            sas=None,
        )
        assert pf.max_population_af == 0.02

    def test_max_population_af_all_none(self) -> None:
        pf = PopulationFrequencies(
            global_af=None, afr=None, amr=None,
            eas=None, nfe=None, sas=None,
        )
        assert pf.max_population_af is None

    def test_any_exceeds_true(self) -> None:
        pf = PopulationFrequencies(afr=0.06)
        assert pf.any_exceeds(0.05) is True

    def test_any_exceeds_false(self) -> None:
        pf = PopulationFrequencies(afr=0.01)
        assert pf.any_exceeds(0.05) is False

    def test_all_below_true(self) -> None:
        pf = PopulationFrequencies(afr=0.001, amr=0.002)
        assert pf.all_below(0.01) is True

    def test_all_below_false(self) -> None:
        pf = PopulationFrequencies(afr=0.001, nfe=0.05)
        assert pf.all_below(0.01) is False


class TestVariant:
    def test_creation(self) -> None:
        v = Variant(chrom="chr17", pos=7577120, id=None, ref="G", alt="A", qual=99.0, filter_status="PASS")
        assert v.chrom == "chr17"
        assert v.pos == 7577120


class TestScoredVariant:
    def test_composite_rank_syncs_to_prioritization_score(self) -> None:
        v = Variant(chrom="chr1", pos=100, id=None, ref="A", alt="T", qual=99.0, filter_status="PASS")
        av = AnnotatedVariant(
            variant=v,
            gene_name="TEST",
            consequence=FunctionalConsequence.MISSENSE,
        )
        scored = ScoredVariant(annotated=av, composite_rank=0.75)
        assert scored.prioritization_score == 0.75

    def test_prioritization_score_syncs_to_composite_rank(self) -> None:
        v = Variant(chrom="chr1", pos=100, id=None, ref="A", alt="T", qual=99.0, filter_status="PASS")
        av = AnnotatedVariant(
            variant=v,
            gene_name="TEST",
            consequence=FunctionalConsequence.MISSENSE,
        )
        scored = ScoredVariant(annotated=av, prioritization_score=0.6)
        assert scored.composite_rank == 0.6


# --------------------------------------------------------------------------
# models/cohort.py
# --------------------------------------------------------------------------


class TestCohortConfig:
    def test_valid_config(self, tmp_path: Path) -> None:
        config = CohortConfig(
            sample_vcfs=[tmp_path / "a.vcf", tmp_path / "b.vcf"],
            output_path=tmp_path / "out",
        )
        assert config.sample_count == 2

    def test_fewer_than_two_samples_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="at least 2 samples"):
            CohortConfig(
                sample_vcfs=[tmp_path / "a.vcf"],
                output_path=tmp_path / "out",
            )

    def test_min_recurrence_zero_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="min_recurrence"):
            CohortConfig(
                sample_vcfs=[tmp_path / "a.vcf", tmp_path / "b.vcf"],
                output_path=tmp_path / "out",
                min_recurrence=0,
            )

    def test_max_af_threshold_out_of_range(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="max_af_threshold"):
            CohortConfig(
                sample_vcfs=[tmp_path / "a.vcf", tmp_path / "b.vcf"],
                output_path=tmp_path / "out",
                max_af_threshold=2.0,
            )

    def test_max_workers_zero_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="max_workers"):
            CohortConfig(
                sample_vcfs=[tmp_path / "a.vcf", tmp_path / "b.vcf"],
                output_path=tmp_path / "out",
                max_workers=0,
            )

    def test_label_for_uses_mapping(self, tmp_path: Path) -> None:
        config = CohortConfig(
            sample_vcfs=[tmp_path / "a.vcf", tmp_path / "b.vcf"],
            output_path=tmp_path / "out",
            sample_labels={"a": "Patient A"},
        )
        assert config.label_for(tmp_path / "a.vcf") == "Patient A"

    def test_label_for_falls_back_to_stem(self, tmp_path: Path) -> None:
        config = CohortConfig(
            sample_vcfs=[tmp_path / "sample1.vcf", tmp_path / "b.vcf"],
            output_path=tmp_path / "out",
        )
        assert config.label_for(tmp_path / "sample1.vcf") == "sample1"

    def test_label_for_strips_vcf_extension(self, tmp_path: Path) -> None:
        config = CohortConfig(
            sample_vcfs=[tmp_path / "s.vcf.gz", tmp_path / "b.vcf"],
            output_path=tmp_path / "out",
        )
        # stem of "s.vcf.gz" is "s.vcf", which should strip .vcf → "s"
        assert config.label_for(tmp_path / "s.vcf.gz") == "s"


class TestCohortVariant:
    def test_cohort_frequency(self) -> None:
        cv = CohortVariant(
            chrom="chr1", pos=100, ref="A", alt="T",
            gene_name="TP53",
            consequence=FunctionalConsequence.MISSENSE,
            sample_count=3, total_samples=10,
            occurrences=(), max_classification=ACMGClassification.VUS,
            all_evidence_tags=frozenset(),
        )
        assert cv.cohort_frequency == pytest.approx(0.3)

    def test_is_singleton(self) -> None:
        cv = CohortVariant(
            chrom="chr1", pos=100, ref="A", alt="T",
            gene_name=None,
            consequence=FunctionalConsequence.MISSENSE,
            sample_count=1, total_samples=5,
            occurrences=(), max_classification=ACMGClassification.VUS,
            all_evidence_tags=frozenset(),
        )
        assert cv.is_singleton is True

    def test_is_universal(self) -> None:
        cv = CohortVariant(
            chrom="chr1", pos=100, ref="A", alt="T",
            gene_name=None,
            consequence=FunctionalConsequence.SYNONYMOUS,
            sample_count=5, total_samples=5,
            occurrences=(), max_classification=ACMGClassification.BENIGN,
            all_evidence_tags=frozenset(),
        )
        assert cv.is_universal is True

    def test_key_property(self) -> None:
        cv = CohortVariant(
            chrom="chrX", pos=999, ref="C", alt="G",
            gene_name=None,
            consequence=FunctionalConsequence.MISSENSE,
            sample_count=2, total_samples=4,
            occurrences=(), max_classification=ACMGClassification.VUS,
            all_evidence_tags=frozenset(),
        )
        assert cv.key == ("chrX", 999, "C", "G")


class TestGeneBurden:
    def test_penetrance(self) -> None:
        gb = GeneBurden(
            gene_name="BRCA1",
            total_variants=5,
            pathogenic_count=2,
            samples_affected=3,
            total_samples=10,
            most_severe=ACMGClassification.PATHOGENIC,
        )
        assert gb.penetrance == pytest.approx(0.3)

    def test_penetrance_zero_samples(self) -> None:
        gb = GeneBurden(
            gene_name="X",
            total_variants=0,
            pathogenic_count=0,
            samples_affected=0,
            total_samples=0,
            most_severe=ACMGClassification.BENIGN,
        )
        assert gb.penetrance == 0.0


# --------------------------------------------------------------------------
# models/config.py
# --------------------------------------------------------------------------


class TestQualityFilterConfig:
    def test_valid_defaults(self) -> None:
        config = QualityFilterConfig()
        assert config.min_qual >= 0

    def test_negative_qual_raises(self) -> None:
        with pytest.raises(ValueError):
            QualityFilterConfig(min_qual=-1)


class TestAnnotationConfig:
    def test_valid_defaults(self, tmp_path: Path) -> None:
        config = AnnotationConfig(
            gene_annotation_path=tmp_path / "genes.gtf",
            gnomad_path=tmp_path / "gnomad.vcf.gz",
        )
        assert config.batch_size == 10_000

    def test_batch_size_out_of_range(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="batch_size"):
            AnnotationConfig(
                gene_annotation_path=tmp_path / "genes.gtf",
                gnomad_path=tmp_path / "gnomad.vcf.gz",
                batch_size=500,
            )


class TestPrioritizationConfig:
    def test_valid_defaults(self) -> None:
        config = PrioritizationConfig()
        assert config.max_allele_frequency == 0.01

    def test_max_af_out_of_range(self) -> None:
        with pytest.raises(ValueError, match="max_allele_frequency"):
            PrioritizationConfig(max_allele_frequency=-0.1)


class TestReportConfig:
    def test_valid_json_format(self) -> None:
        config = ReportConfig(output_format="json")
        assert config.output_format == "json"


class TestPipelineConfig:
    def test_valid_config(self, tmp_path: Path) -> None:
        config = PipelineConfig(
            vcf_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.json",
        )
        assert config.vcf_path == tmp_path / "input.vcf"
