"""Memory and performance integration tests for the variant prioritization pipeline.

Tests verify:
- Pipeline RSS stays below 2GB on 4M+ variant synthetic datasets
- Annotation completes within 30 seconds for 10K variants
- Report generation completes within 10 seconds for 10K variants

These tests are marked @pytest.mark.slow so they can be skipped during fast CI runs.
"""

from __future__ import annotations

import platform
import resource
import tempfile
import time
from pathlib import Path

import pytest

from vartriage.annotation.engine import AnnotationEngine
from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.config import (AnnotationConfig, QualityFilterConfig,
                                     ReportConfig)
from vartriage.models.variant import (ACMGClassification, AnnotatedVariant,
                                      ClassifiedVariant, ClinVarAssertion,
                                      EvidenceTag, FunctionalConsequence,
                                      ScoredVariant, Variant)
from vartriage.reporting.generator import ReportGenerator


def _get_rss_bytes() -> int:
    """Get current maximum RSS in bytes, cross-platform."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    if platform.system() == "Darwin":
        return usage.ru_maxrss  # macOS returns bytes
    return usage.ru_maxrss * 1024  # Linux returns KB


def _write_synthetic_vcf(path: Path, num_variants: int) -> None:
    """Write a syntactically valid VCF with the given number of data lines.

    The file is minimal but parseable by pysam. Each line uses simple
    single-nucleotide variants on chr1 with incrementing positions.
    """
    chromosomes = [f"chr{i}" for i in range(1, 23)]
    with open(path, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i in range(num_variants):
            chrom = chromosomes[i % len(chromosomes)]
            pos = (i // 22) + 1
            f.write(f"{chrom}\t{pos}\tvar_{i}\tA\tG\t50\tPASS\tDP=30\n")


def _write_minimal_gtf(path: Path) -> None:
    """Write a minimal GTF with a single gene for annotation testing."""
    lines = [
        "##format: gtf",
        "chr1\thavana\tgene\t1\t100000\t.\t+\t.\t"
        'gene_id "TEST1"; gene_name "TEST1";',
        "chr1\thavana\ttranscript\t1\t100000\t.\t+\t.\t"
        'gene_id "TEST1"; transcript_id "TEST1.1"; gene_name "TEST1";',
        "chr1\thavana\texon\t1\t100000\t.\t+\t.\t"
        'gene_id "TEST1"; transcript_id "TEST1.1"; gene_name "TEST1";',
        "chr1\thavana\tCDS\t100\t99900\t.\t+\t0\t"
        'gene_id "TEST1"; transcript_id "TEST1.1"; gene_name "TEST1";',
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_minimal_gnomad(path: Path) -> None:
    """Write a minimal gnomAD TSV with a few entries."""
    lines = ["chrom\tpos\tref\talt\taf\n"]
    for i in range(100):
        lines.append(f"chr1\t{i + 1}\tA\tG\t0.001\n")
    path.write_text("".join(lines), encoding="utf-8")


def _write_minimal_clinvar(path: Path) -> None:
    """Write a minimal ClinVar TSV with a few entries."""
    lines = ["chrom\tpos\tref\talt\tclinical_significance\n"]
    for i in range(50):
        lines.append(f"chr1\t{i + 1}\tA\tG\tPathogenic\n")
    path.write_text("".join(lines), encoding="utf-8")


def _build_classified_variants(count: int) -> list[ClassifiedVariant]:
    """Build a list of synthetic ClassifiedVariant objects for report testing."""
    variants: list[ClassifiedVariant] = []
    consequences = list(FunctionalConsequence)

    for i in range(count):
        raw = Variant(
            chrom=f"chr{(i % 22) + 1}",
            pos=i + 1,
            id=f"var_{i}",
            ref="A",
            alt="G",
            qual=50.0,
            filter_status="PASS",
            info={"DP": 30},
        )
        annotated = AnnotatedVariant(
            variant=raw,
            consequence=consequences[i % len(consequences)],
            allele_frequency=0.001 if i % 3 != 0 else None,
            clinvar_assertion=ClinVarAssertion.PATHOGENIC if i % 5 == 0 else None,
            frequency_unknown=(i % 3 == 0),
            clinvar_unknown=(i % 5 != 0),
        )
        scored = ScoredVariant(
            annotated=annotated,
            cadd_phred=25.0 if i % 2 == 0 else None,
            cadd_normalized=min(25.0 / 99.0, 1.0) if i % 2 == 0 else None,
            revel_score=0.75 if i % 3 != 2 else None,
            composite_rank=0.65 if i % 2 == 0 else 0.5,
        )
        classified = ClassifiedVariant(
            scored=scored,
            evidence_tags=frozenset({EvidenceTag.PP3}),
            classification=ACMGClassification.VUS,
            missing_data_sources=frozenset(),
        )
        variants.append(classified)

    return variants


# ---------------------------------------------------------------------------
# Memory Bound Tests
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.skip(reason="requires 4M+ variant file generation. Run manually")
class TestMemoryBounds:
    """Verify pipeline RSS stays below 2GB on large datasets."""

    TWO_GB = 2 * 1024 * 1024 * 1024  # 2 GB in bytes
    VARIANT_COUNT = 4_000_001  # 4M+ variants

    def test_pipeline_rss_below_2gb_on_4m_variants(self, tmp_path: Path) -> None:
        """Streaming 4M+ variants through parse + filter keeps RSS < 2GB.

        This test generates a synthetic VCF with 4M+ simple data lines and
        streams them through the VCF parser and quality filter, verifying
        that peak RSS never exceeds the 2GB memory budget.
        """
        vcf_path = tmp_path / "large_input.vcf"
        _write_synthetic_vcf(vcf_path, self.VARIANT_COUNT)

        rss_before = _get_rss_bytes()

        parser = VCFParser(vcf_path)
        qf = QualityFilter(QualityFilterConfig(min_qual=20.0))

        count = 0
        for variant in qf.apply(iter(parser)):
            count += 1
            if count % 500_000 == 0:
                current_rss = _get_rss_bytes()
                assert current_rss < self.TWO_GB, (
                    f"RSS exceeded 2GB at {count} variants: "
                    f"{current_rss / (1024**3):.2f} GB"
                )

        parser.close()

        final_rss = _get_rss_bytes()
        assert final_rss < self.TWO_GB, (
            f"Final RSS exceeded 2GB after processing {count} variants: "
            f"{final_rss / (1024**3):.2f} GB"
        )
        assert count == self.VARIANT_COUNT


# ---------------------------------------------------------------------------
# Performance Tests
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestAnnotationPerformance:
    """Verify annotation completes within 30 seconds for 10K variants."""

    VARIANT_COUNT = 10_000
    TIME_LIMIT_SECONDS = 30

    def test_annotation_completes_within_30_seconds(self, tmp_path: Path) -> None:
        """Annotating 10K variants against a small reference finishes in < 30s.

        Generates a synthetic VCF with 10K variants, builds minimal
        reference files, and measures annotation wall-clock time.
        """
        vcf_path = tmp_path / "input_10k.vcf"
        _write_synthetic_vcf(vcf_path, self.VARIANT_COUNT)

        gtf_path = tmp_path / "genes.gtf"
        _write_minimal_gtf(gtf_path)

        gnomad_path = tmp_path / "gnomad.tsv"
        _write_minimal_gnomad(gnomad_path)

        clinvar_path = tmp_path / "clinvar.tsv"
        _write_minimal_clinvar(clinvar_path)

        parser = VCFParser(vcf_path)
        variants = list(parser)
        parser.close()
        assert len(variants) == self.VARIANT_COUNT

        ann_config = AnnotationConfig(
            gene_annotation_path=gtf_path,
            gnomad_path=gnomad_path,
            clinvar_path=clinvar_path,
            batch_size=5_000,
        )
        engine = AnnotationEngine(ann_config)

        start_time = time.time()
        annotated = list(engine.annotate(iter(variants)))
        elapsed = time.time() - start_time

        assert len(annotated) == self.VARIANT_COUNT
        assert elapsed < self.TIME_LIMIT_SECONDS, (
            f"Annotation took {elapsed:.2f}s for {self.VARIANT_COUNT} variants, "
            f"exceeding the {self.TIME_LIMIT_SECONDS}s limit"
        )


@pytest.mark.slow
class TestReportGenerationPerformance:
    """Verify report generation completes within 10 seconds for 10K variants."""

    VARIANT_COUNT = 10_000
    TIME_LIMIT_SECONDS = 10

    def test_json_report_within_10_seconds(self, tmp_path: Path) -> None:
        """JSON serialization of 10K classified variants finishes in < 10s."""
        variants = _build_classified_variants(self.VARIANT_COUNT)

        config = ReportConfig(output_format="json")
        generator = ReportGenerator(config)
        output_path = tmp_path / "report.json"

        start_time = time.time()
        generator.generate(variants, output_path)
        elapsed = time.time() - start_time

        assert output_path.exists()
        assert output_path.stat().st_size > 0
        assert elapsed < self.TIME_LIMIT_SECONDS, (
            f"JSON report generation took {elapsed:.2f}s for "
            f"{self.VARIANT_COUNT} variants, "
            f"exceeding the {self.TIME_LIMIT_SECONDS}s limit"
        )

    def test_csv_report_within_10_seconds(self, tmp_path: Path) -> None:
        """CSV serialization of 10K classified variants finishes in < 10s."""
        variants = _build_classified_variants(self.VARIANT_COUNT)

        config = ReportConfig(output_format="csv")
        generator = ReportGenerator(config)
        output_path = tmp_path / "report.csv"

        start_time = time.time()
        generator.generate(variants, output_path)
        elapsed = time.time() - start_time

        assert output_path.exists()
        assert output_path.stat().st_size > 0
        assert elapsed < self.TIME_LIMIT_SECONDS, (
            f"CSV report generation took {elapsed:.2f}s for "
            f"{self.VARIANT_COUNT} variants, "
            f"exceeding the {self.TIME_LIMIT_SECONDS}s limit"
        )
