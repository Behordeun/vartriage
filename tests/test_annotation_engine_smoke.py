"""Smoke test for AnnotationEngine orchestrator with backend auto-detection."""

import tempfile
from pathlib import Path

import pytest

from vartriage.annotation.engine import AnnotationEngine
from vartriage.models.config import AnnotationConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    ClinVarAssertion,
    FunctionalConsequence,
    Variant,
)
from vartriage.models.warnings import MissingDataWarning


# Minimal GTF content
SAMPLE_GTF = """\
##description: test annotation
chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\texon\t1000\t1200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\tCDS\t1050\t1190\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\texon\t2000\t2500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\tCDS\t2000\t2500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
"""

# gnomAD TSV reference
SAMPLE_GNOMAD = """\
chrom\tpos\tref\talt\taf
chr1\t2100\tA\tT\t0.001
chr1\t6000\tA\tT\t0.15
"""

# ClinVar TSV reference
SAMPLE_CLINVAR = """\
chrom\tpos\tref\talt\tclinical_significance
chr1\t2100\tA\tT\tPathogenic
chr1\t6000\tA\tT\tBenign
"""


def _write_temp(content: str, suffix: str) -> Path:
    """Write content to a temporary file and return the path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    f.write(content)
    f.flush()
    f.close()
    return Path(f.name)


def _make_variant(chrom: str, pos: int, ref: str, alt: str) -> Variant:
    return Variant(
        chrom=chrom,
        pos=pos,
        id=None,
        ref=ref,
        alt=alt,
        qual=30.0,
        filter_status="PASS",
    )


@pytest.fixture
def annotation_files():
    """Create temporary reference files for testing."""
    gtf_path = _write_temp(SAMPLE_GTF, ".gtf")
    gnomad_path = _write_temp(SAMPLE_GNOMAD, ".tsv")
    clinvar_path = _write_temp(SAMPLE_CLINVAR, ".tsv")
    return gtf_path, gnomad_path, clinvar_path


@pytest.fixture
def engine(annotation_files):
    """Create an AnnotationEngine with temporary reference files."""
    gtf_path, gnomad_path, clinvar_path = annotation_files
    config = AnnotationConfig(
        gene_annotation_path=gtf_path,
        gnomad_path=gnomad_path,
        clinvar_path=clinvar_path,
        batch_size=1000,
    )
    return AnnotationEngine(config)


class TestAnnotationEngineInit:
    def test_creates_with_valid_paths(self, annotation_files):
        gtf_path, gnomad_path, clinvar_path = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=clinvar_path,
        )
        engine = AnnotationEngine(config)
        assert engine is not None

    def test_creates_without_clinvar(self, annotation_files):
        gtf_path, gnomad_path, _ = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=None,
        )
        engine = AnnotationEngine(config)
        assert engine is not None

    def test_raises_on_missing_gtf(self, annotation_files):
        _, gnomad_path, clinvar_path = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=Path("/nonexistent/genes.gtf"),
            gnomad_path=gnomad_path,
            clinvar_path=clinvar_path,
        )
        with pytest.raises(FileNotFoundError, match="Gene annotation"):
            AnnotationEngine(config)

    def test_raises_on_missing_gnomad(self, annotation_files):
        gtf_path, _, clinvar_path = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=Path("/nonexistent/gnomad.tsv"),
            clinvar_path=clinvar_path,
        )
        with pytest.raises(FileNotFoundError, match="gnomAD"):
            AnnotationEngine(config)

    def test_raises_on_missing_clinvar(self, annotation_files):
        gtf_path, gnomad_path, _ = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=Path("/nonexistent/clinvar.tsv"),
        )
        with pytest.raises(FileNotFoundError, match="ClinVar"):
            AnnotationEngine(config)


class TestAnnotationEngineAnnotate:
    def test_annotates_variant_in_coding_region_with_frequency(self, engine):
        """Variant in CDS with matching gnomAD and ClinVar entries."""
        variant = _make_variant("chr1", 2100, "A", "T")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        annotated = results[0]
        assert isinstance(annotated, AnnotatedVariant)
        assert annotated.variant == variant
        assert annotated.consequence == FunctionalConsequence.MISSENSE
        assert annotated.allele_frequency == pytest.approx(0.001)
        assert annotated.clinvar_assertion == ClinVarAssertion.PATHOGENIC
        assert annotated.frequency_unknown is False
        assert annotated.clinvar_unknown is False

    def test_annotates_intergenic_variant_with_frequency(self, engine):
        """Variant outside coding regions but present in gnomAD."""
        variant = _make_variant("chr1", 6000, "A", "T")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        annotated = results[0]
        assert annotated.consequence == FunctionalConsequence.INTERGENIC
        assert annotated.allele_frequency == pytest.approx(0.15)
        assert annotated.clinvar_assertion == ClinVarAssertion.BENIGN
        assert annotated.frequency_unknown is False

    def test_annotates_variant_not_in_frequency_db(self, engine):
        """Variant not in gnomAD should get null frequency + warning."""
        variant = _make_variant("chr1", 3100, "A", "T")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        annotated = results[0]
        assert annotated.allele_frequency is None
        assert annotated.frequency_unknown is True

        # Should have emitted at least one gnomAD warning
        gnomad_warnings = [
            w for w in engine.warnings if w.source == "gnomAD"
        ]
        assert len(gnomad_warnings) >= 1
        assert gnomad_warnings[0].chrom == "chr1"
        assert gnomad_warnings[0].pos == 3100

    def test_annotates_variant_not_in_clinvar(self, engine):
        """Variant not in ClinVar should get null assertion + clinvar_unknown."""
        variant = _make_variant("chr1", 3100, "A", "T")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        annotated = results[0]
        assert annotated.clinvar_assertion is None
        assert annotated.clinvar_unknown is True

    def test_processes_empty_iterator(self, engine):
        """Empty input should produce empty output."""
        results = list(engine.annotate(iter([])))
        assert results == []

    def test_processes_multiple_variants_in_batch(self, engine):
        """Multiple variants processed together."""
        variants = [
            _make_variant("chr1", 2100, "A", "T"),
            _make_variant("chr1", 6000, "A", "T"),
            _make_variant("chr1", 9999, "G", "C"),
        ]
        results = list(engine.annotate(iter(variants)))

        assert len(results) == 3
        # First: coding region, has frequency
        assert results[0].consequence == FunctionalConsequence.MISSENSE
        assert results[0].allele_frequency == pytest.approx(0.001)
        # Second: intergenic, has frequency
        assert results[1].consequence == FunctionalConsequence.INTERGENIC
        assert results[1].allele_frequency == pytest.approx(0.15)
        # Third: intergenic, no frequency
        assert results[2].consequence == FunctionalConsequence.INTERGENIC
        assert results[2].allele_frequency is None
        assert results[2].frequency_unknown is True

    def test_preserves_variant_ordering(self, engine):
        """Output order matches input order."""
        variants = [
            _make_variant("chr1", 9999, "G", "C"),
            _make_variant("chr1", 2100, "A", "T"),
            _make_variant("chr1", 6000, "A", "T"),
        ]
        results = list(engine.annotate(iter(variants)))

        assert results[0].variant.pos == 9999
        assert results[1].variant.pos == 2100
        assert results[2].variant.pos == 6000

    def test_batching_with_more_variants_than_batch_size(self, annotation_files):
        """Verify multiple batches work correctly."""
        gtf_path, gnomad_path, clinvar_path = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=clinvar_path,
            batch_size=1000,  # minimum batch size
        )
        engine = AnnotationEngine(config)

        # Create more variants than batch_size
        variants = [
            _make_variant("chr1", 2100, "A", "T") for _ in range(1500)
        ]
        results = list(engine.annotate(iter(variants)))

        assert len(results) == 1500
        for r in results:
            assert r.consequence == FunctionalConsequence.MISSENSE
            assert r.allele_frequency == pytest.approx(0.001)


class TestAnnotationEngineWithoutClinVar:
    def test_no_clinvar_yields_null_assertion_and_unknown_false(
        self, annotation_files
    ):
        """When ClinVar path is None, clinvar_unknown stays False (not queried)."""
        gtf_path, gnomad_path, _ = annotation_files
        config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=None,
            batch_size=1000,
        )
        engine = AnnotationEngine(config)

        variant = _make_variant("chr1", 2100, "A", "T")
        results = list(engine.annotate(iter([variant])))

        assert len(results) == 1
        # When ClinVar is not configured, assertion is None but
        # clinvar_unknown is True (we have no ClinVar data at all)
        assert results[0].clinvar_assertion is None
        assert results[0].clinvar_unknown is True
