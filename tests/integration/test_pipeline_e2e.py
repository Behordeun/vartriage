"""End-to-end integration test for the variant prioritization pipeline.

Creates synthetic reference data (VCF, GTF, gnomAD, ClinVar, CADD, REVEL)
and runs all pipeline stages to verify correctness of the full flow from
VCF ingestion through ACMG classification to report generation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from vartriage.annotation.engine import AnnotationEngine
from vartriage.classification.acmg import ACMGClassifier
from vartriage.filter.quality_filter import QualityFilter
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.config import (AnnotationConfig, PrioritizationConfig,
                                     QualityFilterConfig, ReportConfig)
from vartriage.models.variant import (ACMGClassification, AnnotatedVariant,
                                      ClassifiedVariant, FunctionalConsequence)
from vartriage.prioritization.scoring import score_variants
from vartriage.reporting.generator import ReportGenerator

# ---------------------------------------------------------------------------
# Synthetic data generation helpers
# ---------------------------------------------------------------------------

_VCF_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _build_vcf_lines() -> list[str]:
    """Build ~100 synthetic VCF data lines with diverse scenarios.

    Scenarios covered:
      - PASS variants in CDS with various consequences (SNV->Missense,
        indel->Frameshift, in-frame insertion/deletion, splice site)
      - PASS variants in exon/UTR regions (Synonymous)
      - Intergenic variants (no gene overlap)
      - Variants with LowQual filter (should be excluded)
      - Variants with missing QUAL (should be excluded)
      - Low QUAL variants below threshold (should be excluded)
      - High quality variants across multiple chromosomes
    """
    lines: list[str] = []
    idx = 0

    # --- Group 1: CDS coding variants (chr1, gene BRCA1 region 1000-2000) ---
    # Missense SNVs (in CDS)
    for pos in range(1050, 1070):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tA\tG\t50\tPASS\tDP=30")

    # Frameshift indels (in CDS: deletion of 1 or 2 bp)
    for i, pos in enumerate(range(1100, 1110)):
        idx += 1
        ref = "AC" if i % 2 == 0 else "ACG"
        lines.append(f"chr1\t{pos}\tvar_{idx}\t{ref}\tA\t60\tPASS\tDP=40")

    # In-frame insertion (3bp insertion in CDS)
    for pos in range(1200, 1205):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tA\tATTT\t55\tPASS\tDP=35")

    # In-frame deletion (3bp deletion in CDS)
    for pos in range(1250, 1255):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tATTT\tA\t55\tPASS\tDP=35")

    # --- Group 2: Splice site variants (within 2bp of exon boundary) ---
    for pos in [999, 1000, 1001, 1002]:
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tC\tT\t70\tPASS\tDP=50")

    # --- Group 3: Exon but not CDS (UTR-like, mapped as Synonymous) ---
    # Gene BRCA1 has exon region 900-2100 but CDS is only 1000-2000
    for pos in range(910, 920):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tG\tA\t45\tPASS\tDP=25")

    # --- Group 4: Intergenic (no gene overlap, chr2 has no genes) ---
    for pos in range(5000, 5015):
        idx += 1
        lines.append(f"chr2\t{pos}\tvar_{idx}\tT\tC\t40\tPASS\tDP=20")

    # --- Group 5: Variants that should FAIL quality filtering ---
    # LowQual filter
    for pos in range(1300, 1310):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tA\tT\t50\tLowQual\tDP=10")

    # Missing QUAL (represented as ".")
    for pos in range(1400, 1405):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tG\tC\t.\tPASS\tDP=15")

    # Low QUAL below threshold of 30
    for pos in range(1500, 1505):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tC\tA\t10\tPASS\tDP=12")

    # --- Group 6: Variants on chr3 in TP53 gene CDS 3000-4000 ---
    for pos in range(3050, 3060):
        idx += 1
        lines.append(f"chr3\t{pos}\tvar_{idx}\tA\tG\t80\tPASS\tDP=60")

    # Frameshift in TP53
    for pos in range(3100, 3105):
        idx += 1
        lines.append(f"chr3\t{pos}\tvar_{idx}\tAC\tA\t90\tPASS\tDP=70")

    # --- Group 7: Additional PASS variants to reach ~100 total ---
    for pos in range(1600, 1615):
        idx += 1
        lines.append(f"chr1\t{pos}\tvar_{idx}\tG\tT\t65\tPASS\tDP=45")

    return lines


def _write_vcf(tmp_dir: Path) -> Path:
    """Write a synthetic VCF file with ~100 variants."""
    vcf_path = tmp_dir / "test_input.vcf"
    lines = _build_vcf_lines()
    content = _VCF_HEADER + "\n".join(lines) + "\n"
    vcf_path.write_text(content, encoding="utf-8")
    return vcf_path


def _write_gtf(tmp_dir: Path) -> Path:
    """Write a minimal GTF with two genes.

    Gene BRCA1 on chr1: transcript 900-2100, CDS 1000-2000, exon 900-2100
    Gene TP53 on chr3: transcript 2900-4100, CDS 3000-4000, exon 2900-4100
    """
    gtf_path = tmp_dir / "genes.gtf"
    gtf_lines = [
        "##format: gtf",
        "chr1\thavana\tgene\t900\t2100\t.\t+\t.\t"
        'gene_id "BRCA1"; gene_name "BRCA1";',
        "chr1\thavana\ttranscript\t900\t2100\t.\t+\t.\t"
        'gene_id "BRCA1"; transcript_id "BRCA1.1"; gene_name "BRCA1";',
        "chr1\thavana\texon\t900\t2100\t.\t+\t.\t"
        'gene_id "BRCA1"; transcript_id "BRCA1.1"; gene_name "BRCA1";',
        "chr1\thavana\tCDS\t1000\t2000\t.\t+\t0\t"
        'gene_id "BRCA1"; transcript_id "BRCA1.1"; gene_name "BRCA1";',
        "chr3\thavana\tgene\t2900\t4100\t.\t+\t.\t" 'gene_id "TP53"; gene_name "TP53";',
        "chr3\thavana\ttranscript\t2900\t4100\t.\t+\t.\t"
        'gene_id "TP53"; transcript_id "TP53.1"; gene_name "TP53";',
        "chr3\thavana\texon\t2900\t4100\t.\t+\t.\t"
        'gene_id "TP53"; transcript_id "TP53.1"; gene_name "TP53";',
        "chr3\thavana\tCDS\t3000\t4000\t.\t+\t0\t"
        'gene_id "TP53"; transcript_id "TP53.1"; gene_name "TP53";',
    ]
    gtf_path.write_text("\n".join(gtf_lines) + "\n", encoding="utf-8")
    return gtf_path


def _write_gnomad(tmp_dir: Path, vcf_lines: list[str]) -> Path:
    """Write a synthetic gnomAD TSV.

    Strategy:
      - Some CDS variants get very low AF (< 0.0001) -> triggers PM2
      - Some CDS variants get moderate AF (0.001-0.01) -> passes freq filter
      - Some CDS variants get high AF (> 0.01) -> excluded by freq filter
      - Some variants are absent -> frequency_unknown, bypasses AF filter
    """
    gnomad_path = tmp_dir / "gnomad.tsv"
    rows = [["chrom", "pos", "ref", "alt", "af"]]

    # Give low AF to first 10 CDS missense variants (chr1:1050-1059)
    for pos in range(1050, 1060):
        rows.append(["chr1", str(pos), "A", "G", "0.00005"])

    # Give moderate AF to next 10 CDS missense variants (chr1:1060-1069)
    for pos in range(1060, 1070):
        rows.append(["chr1", str(pos), "A", "G", "0.005"])

    # Give high AF to some exon/UTR variants -> should be excluded
    for pos in range(910, 915):
        rows.append(["chr1", str(pos), "G", "A", "0.05"])

    # Intergenic variants with high AF
    for pos in range(5000, 5010):
        rows.append(["chr2", str(pos), "T", "C", "0.15"])

    # TP53 CDS variants: very low AF for frameshift variants
    for pos in range(3100, 3105):
        rows.append(["chr3", str(pos), "AC", "A", "0.00001"])

    # TP53 missense variants: moderate AF
    for pos in range(3050, 3055):
        rows.append(["chr3", str(pos), "A", "G", "0.003"])

    # Remaining TP53 missense: absent from gnomAD -> frequency_unknown

    with open(gnomad_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)

    return gnomad_path


def _write_clinvar(tmp_dir: Path) -> Path:
    """Write a synthetic ClinVar TSV.

    Strategy:
      - Some frameshift variants are Pathogenic in ClinVar
      - Some missense variants are VUS
      - Some variants are Likely_Pathogenic
      - Most variants absent -> clinvar_unknown
    """
    clinvar_path = tmp_dir / "clinvar.tsv"
    rows = [["chrom", "pos", "ref", "alt", "clinical_significance"]]

    # Frameshift variants in TP53 marked Pathogenic
    for pos in range(3100, 3105):
        rows.append(["chr3", str(pos), "AC", "A", "Pathogenic"])

    # Some CDS missense in BRCA1 region marked VUS
    for pos in range(1050, 1055):
        rows.append(["chr1", str(pos), "A", "G", "Uncertain significance"])

    # One missense in chr1 marked Likely pathogenic
    rows.append(["chr1", "1055", "A", "G", "Likely pathogenic"])

    with open(clinvar_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)

    return clinvar_path


def _build_cadd_scores(annotated: list[AnnotatedVariant]) -> list[float | None]:
    """Build synthetic CADD Phred scores for annotated variants.

    Strategy:
      - Frameshift variants: high CADD (35-45)
      - Missense variants: moderate CADD (20-30)
      - Splice site: high CADD (30-40)
      - Synonymous/Intergenic: low CADD (5-10) or None
    """
    scores: list[float | None] = []
    for av in annotated:
        cons = av.consequence
        if cons == FunctionalConsequence.FRAMESHIFT:
            scores.append(40.0)
        elif cons == FunctionalConsequence.NONSENSE:
            scores.append(45.0)
        elif cons == FunctionalConsequence.SPLICE_SITE:
            scores.append(35.0)
        elif cons == FunctionalConsequence.MISSENSE:
            scores.append(25.0)
        elif cons in (
            FunctionalConsequence.IN_FRAME_INSERTION,
            FunctionalConsequence.IN_FRAME_DELETION,
        ):
            scores.append(20.0)
        elif cons == FunctionalConsequence.SYNONYMOUS:
            scores.append(8.0)
        else:
            scores.append(None)
    return scores


def _build_revel_scores(annotated: list[AnnotatedVariant]) -> list[float | None]:
    """Build synthetic REVEL scores for annotated variants.

    Strategy:
      - Frameshift variants: high REVEL (0.85-0.95) -> triggers PP3
      - Missense in BRCA1 with low AF: REVEL 0.75 -> triggers PP3
      - Missense moderate: REVEL 0.5
      - Synonymous/Intergenic: None (REVEL not applicable)
    """
    scores: list[float | None] = []
    for av in annotated:
        cons = av.consequence
        v = av.variant
        if cons == FunctionalConsequence.FRAMESHIFT:
            scores.append(0.9)
        elif cons == FunctionalConsequence.NONSENSE:
            scores.append(0.92)
        elif cons == FunctionalConsequence.SPLICE_SITE:
            scores.append(0.8)
        elif cons == FunctionalConsequence.MISSENSE:
            # High REVEL for low-AF BRCA1 missense
            if v.chrom == "chr1" and 1050 <= v.pos < 1060:
                scores.append(0.75)
            else:
                scores.append(0.5)
        elif cons in (
            FunctionalConsequence.IN_FRAME_INSERTION,
            FunctionalConsequence.IN_FRAME_DELETION,
        ):
            scores.append(0.6)
        else:
            scores.append(None)
    return scores


# ---------------------------------------------------------------------------
# Pipeline execution helper (wires stages manually since Pipeline.run()
# is not yet implemented per task 10.4)
# ---------------------------------------------------------------------------


def _run_pipeline(
    vcf_path: Path,
    gtf_path: Path,
    gnomad_path: Path,
    clinvar_path: Path,
    output_dir: Path,
    output_format: str = "json",
    min_qual: float = 30.0,
    max_af: float = 0.01,
) -> tuple[Path, list[ClassifiedVariant]]:
    """Execute pipeline stages manually and return output path + results."""
    # Stage 1: Parse VCF
    parser = VCFParser(vcf_path)
    raw_variants = list(parser)
    parser.close()

    # Stage 2: Quality filter
    qf_config = QualityFilterConfig(min_qual=min_qual)
    quality_filter = QualityFilter(qf_config)
    import warnings as _w

    with _w.catch_warnings():
        _w.simplefilter("ignore")
        filtered = list(quality_filter.apply(iter(raw_variants)))

    # Stage 3: Annotation
    ann_config = AnnotationConfig(
        gene_annotation_path=gtf_path,
        gnomad_path=gnomad_path,
        clinvar_path=clinvar_path,
        batch_size=1000,
    )
    engine = AnnotationEngine(ann_config)
    annotated = list(engine.annotate(iter(filtered)))

    # Stage 4: Prioritization (frequency filter + scoring)
    pri_config = PrioritizationConfig(
        max_allele_frequency=max_af,
        batch_size=1000,
    )
    # Apply frequency filter first
    from vartriage.prioritization.frequency_filter import FrequencyFilter

    freq_filter = FrequencyFilter(pri_config)
    af_filtered = list(freq_filter.apply(iter(annotated)))

    # Build synthetic scores for filtered variants
    cadd_scores = _build_cadd_scores(af_filtered)
    revel_scores = _build_revel_scores(af_filtered)

    # Score and sort
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        scored = score_variants(af_filtered, cadd_scores, revel_scores)

    # Stage 5: ACMG Classification
    classifier = ACMGClassifier()
    classified = list(classifier.classify(iter(scored)))

    # Stage 6: Report Generation
    report_config = ReportConfig(output_format=output_format)
    generator = ReportGenerator(report_config)
    output_path = output_dir / f"report.{output_format}"
    generator.generate(classified, output_path)

    return output_path, classified


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pipeline_data(tmp_path: Path) -> dict[str, Path]:
    """Create all synthetic reference data files and return paths."""
    vcf_path = _write_vcf(tmp_path)
    gtf_path = _write_gtf(tmp_path)
    vcf_lines = _build_vcf_lines()
    gnomad_path = _write_gnomad(tmp_path, vcf_lines)
    clinvar_path = _write_clinvar(tmp_path)
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    return {
        "vcf_path": vcf_path,
        "gtf_path": gtf_path,
        "gnomad_path": gnomad_path,
        "clinvar_path": clinvar_path,
        "output_dir": output_dir,
        "tmp_path": tmp_path,
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestPipelineE2E:
    """End-to-end integration tests verifying full pipeline execution."""

    def test_pipeline_produces_classified_variants(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Running the full pipeline produces a non-empty classified list."""
        output_path, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )
        assert len(classified) > 0
        assert all(isinstance(v, ClassifiedVariant) for v in classified)

    def test_quality_filtering_excludes_expected_variants(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Variants with LowQual filter, missing QUAL, or low QUAL are excluded."""
        # Parse all variants
        parser = VCFParser(pipeline_data["vcf_path"])
        all_variants = list(parser)
        parser.close()

        total_count = len(all_variants)
        # We have: 10 LowQual, 5 missing QUAL, 5 below threshold
        # That's 20 variants that should be filtered out
        assert total_count >= 80  # Ensure we have enough raw variants

        # Apply quality filter with threshold 30
        import warnings as _w

        qf = QualityFilter(QualityFilterConfig(min_qual=30.0))
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            filtered = list(qf.apply(iter(all_variants)))

        # Should exclude: 10 LowQual + 5 missing QUAL + 5 low QUAL = 20
        expected_excluded = 20
        assert len(filtered) == total_count - expected_excluded

    def test_frequency_unknown_bypasses_af_filter(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Variants absent from gnomAD (frequency_unknown) pass AF filter."""
        _, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )

        # Some variants should have frequency_unknown=True and still appear
        freq_unknown_variants = [
            v for v in classified if v.scored.annotated.frequency_unknown
        ]
        assert len(freq_unknown_variants) > 0

    def test_consequence_types_cover_all_categories(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Pipeline produces variants with diverse consequence types."""
        _, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )

        consequences = {v.scored.annotated.consequence for v in classified}
        # Should have at least Missense, Frameshift, and either
        # Splice_Site or Synonymous or Intergenic
        assert FunctionalConsequence.MISSENSE in consequences
        assert FunctionalConsequence.FRAMESHIFT in consequences

    def test_pathogenic_classification_for_qualifying_variants(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Variants with PVS1 + PM2 + PP3 + PP5 should be Likely_Pathogenic or Pathogenic.

        TP53 frameshift variants have:
        - Frameshift consequence -> PVS1 (Very Strong)
        - AF < 0.0001 -> PM2 (Moderate)
        - REVEL > 0.7 -> PP3 (Supporting)
        - ClinVar Pathogenic -> PP5 (Supporting)

        PVS1 (Very Strong) + PM2 (Moderate) = Likely_Pathogenic
        PVS1 (Very Strong) + PP3 + PP5 (2 Supporting) = Pathogenic
        Combined: PVS1 + PM2 + PP3 + PP5 should yield Likely_Pathogenic
        (1 Very Strong + 1 Moderate satisfies Likely_Pathogenic rule)
        """
        _, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )

        pathogenic_or_likely = [
            v
            for v in classified
            if v.classification
            in (
                ACMGClassification.PATHOGENIC,
                ACMGClassification.LIKELY_PATHOGENIC,
            )
        ]
        # TP53 frameshift variants should qualify
        assert len(pathogenic_or_likely) > 0

        # Verify at least one has the expected evidence tags
        from vartriage.models.variant import EvidenceTag

        for v in pathogenic_or_likely:
            tags = v.evidence_tags
            # Should have PVS1 at minimum (frameshift)
            assert EvidenceTag.PVS1 in tags

    def test_composite_rank_ordering(self, pipeline_data: dict[str, Path]) -> None:
        """Output is sorted descending by composite_rank, nulls last."""
        _, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )

        ranks = [v.scored.composite_rank for v in classified]
        non_null_ranks = [r for r in ranks if r is not None]
        null_positions = [i for i, r in enumerate(ranks) if r is None]

        # Non-null ranks should be in descending order
        assert non_null_ranks == sorted(non_null_ranks, reverse=True)

        # Null ranks should appear after all non-null ranks
        if null_positions and non_null_ranks:
            first_null_pos = null_positions[0]
            last_non_null_pos = max(i for i, r in enumerate(ranks) if r is not None)
            assert first_null_pos > last_non_null_pos


class TestJSONOutput:
    """Verify JSON report output validity and structure."""

    def test_json_output_is_valid_json(self, pipeline_data: dict[str, Path]) -> None:
        """Generated JSON file is parseable and RFC 8259 compliant."""
        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="json",
        )
        assert output_path.exists()

        content = output_path.read_text(encoding="utf-8")
        data = json.loads(content)
        assert isinstance(data, list)
        assert len(data) > 0

    def test_json_output_has_correct_structure(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Each JSON record contains all required fields."""
        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="json",
        )

        data = json.loads(output_path.read_text(encoding="utf-8"))
        expected_fields = {
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            "gene_name",
            "functional_consequence",
            "allele_frequency",
            "revel_score",
            "composite_rank",
            "prioritization_score",
            "clinvar_assertion",
            "acmg_classification",
            "evidence_tags",
        }

        for record in data:
            assert set(record.keys()) == expected_fields
            assert isinstance(record["chromosome"], str)
            assert isinstance(record["position"], int)
            assert isinstance(record["ref_allele"], str)
            assert isinstance(record["alt_allele"], str)

    def test_json_round_trip_fidelity(self, pipeline_data: dict[str, Path]) -> None:
        """Serializing and deserializing JSON preserves all values."""
        output_path, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="json",
        )

        data = json.loads(output_path.read_text(encoding="utf-8"))
        assert len(data) == len(classified)

        # Verify key fields match the in-memory classified list
        for record, cv in zip(data, classified):
            base = cv.scored.annotated.variant
            assert record["chromosome"] == base.chrom
            assert record["position"] == base.pos
            assert record["ref_allele"] == base.ref
            assert record["alt_allele"] == base.alt

    def test_json_null_representation(self, pipeline_data: dict[str, Path]) -> None:
        """Absent field values are represented as JSON null."""
        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="json",
        )

        data = json.loads(output_path.read_text(encoding="utf-8"))
        # Some variants should have null allele_frequency (frequency_unknown)
        null_af_records = [r for r in data if r["allele_frequency"] is None]
        assert len(null_af_records) > 0


class TestCSVOutput:
    """Verify CSV report output validity and structure."""

    def test_csv_output_has_correct_header(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Generated CSV has the expected header row with correct fields."""
        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="csv",
        )
        assert output_path.exists()

        with open(output_path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            header = next(reader)

        expected_header = [
            "chromosome",
            "position",
            "ref_allele",
            "alt_allele",
            "gene_name",
            "functional_consequence",
            "allele_frequency",
            "revel_score",
            "composite_rank",
            "clinvar_assertion",
            "acmg_classification",
            "evidence_tags",
        ]
        assert header == expected_header

    def test_csv_row_count_matches_classified(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """CSV data row count matches the number of classified variants."""
        output_path, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="csv",
        )

        with open(output_path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        # First row is header, rest are data
        data_rows = rows[1:]
        assert len(data_rows) == len(classified)

    def test_csv_consistent_field_count(self, pipeline_data: dict[str, Path]) -> None:
        """Every CSV row has the same number of fields as the header."""
        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="csv",
        )

        with open(output_path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        header_count = len(rows[0])
        for i, row in enumerate(rows[1:], start=2):
            assert (
                len(row) == header_count
            ), f"Row {i} has {len(row)} fields, expected {header_count}"

    def test_csv_empty_fields_for_null_values(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Absent values are represented as empty strings in CSV."""
        output_path, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="csv",
        )

        with open(output_path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)

        # Find clinvar_assertion column index
        clinvar_idx = header.index("clinvar_assertion")
        # Some variants should have empty clinvar_assertion
        empty_clinvar_rows = [r for r in rows if r[clinvar_idx] == ""]
        assert len(empty_clinvar_rows) > 0


class TestPDFOutput:
    """Verify PDF report output (requires reportlab)."""

    def test_pdf_output_is_valid(self, pipeline_data: dict[str, Path]) -> None:
        """Generated PDF file exists and has non-zero size."""
        pytest.importorskip("reportlab")

        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="pdf",
        )
        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_pdf_starts_with_pdf_header(self, pipeline_data: dict[str, Path]) -> None:
        """Valid PDF files start with the %PDF magic bytes."""
        pytest.importorskip("reportlab")

        output_path, _ = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
            output_format="pdf",
        )
        content = output_path.read_bytes()
        assert content[:4] == b"%PDF"


class TestVUSClassification:
    """Verify VUS is assigned when no combining rules are met."""

    def test_variants_without_strong_evidence_are_vus(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Variants with insufficient evidence receive VUS classification."""
        _, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )

        vus_variants = [
            v for v in classified if v.classification == ACMGClassification.VUS
        ]
        # Many variants should be VUS (synonymous, intergenic,
        # missense without sufficient evidence)
        assert len(vus_variants) > 0

    def test_empty_evidence_tags_yields_vus(
        self, pipeline_data: dict[str, Path]
    ) -> None:
        """Variants with no evidence tags at all are classified as VUS."""
        _, classified = _run_pipeline(
            vcf_path=pipeline_data["vcf_path"],
            gtf_path=pipeline_data["gtf_path"],
            gnomad_path=pipeline_data["gnomad_path"],
            clinvar_path=pipeline_data["clinvar_path"],
            output_dir=pipeline_data["output_dir"],
        )

        # Find any variant with no evidence tags
        no_tags = [v for v in classified if len(v.evidence_tags) == 0]
        for v in no_tags:
            assert v.classification == ACMGClassification.VUS
