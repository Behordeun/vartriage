"""Tests for SVAnnotator reference data loading (GTF, dosage, gnomAD-SV)."""

from __future__ import annotations

from pathlib import Path

import pytest

from vartriage.structural.annotator import (
    SVAnnotator,
    _extract_attribute,
    _parse_score,
)
from vartriage.structural.models import SVConsequence, SVType, StructuralVariant


class TestExtractAttribute:
    def test_quoted_value(self) -> None:
        attrs = 'gene_id "ENSG00000141510"; gene_name "TP53"; gene_type "protein_coding";'
        assert _extract_attribute(attrs, "gene_name") == "TP53"
        assert _extract_attribute(attrs, "gene_type") == "protein_coding"

    def test_missing_key_returns_none(self) -> None:
        attrs = 'gene_id "ENSG00000141510";'
        assert _extract_attribute(attrs, "gene_name") is None

    def test_unquoted_value(self) -> None:
        attrs = "gene_id ENSG00000141510"
        result = _extract_attribute(attrs, "gene_id")
        assert result == "ENSG00000141510"


class TestParseScore:
    def test_valid_float(self) -> None:
        assert _parse_score("3.0") == 3.0

    def test_none_input(self) -> None:
        assert _parse_score(None) is None

    def test_empty_string(self) -> None:
        assert _parse_score("") is None

    def test_na_values(self) -> None:
        assert _parse_score("N/A") is None
        assert _parse_score("NA") is None
        assert _parse_score("-") is None
        assert _parse_score("Not yet evaluated") is None

    def test_invalid_string(self) -> None:
        assert _parse_score("unknown") is None

    def test_whitespace_stripped(self) -> None:
        assert _parse_score("  2.5  ") == 2.5


class TestGTFLoading:
    def test_loads_protein_coding_genes(self, tmp_path: Path) -> None:
        gtf = tmp_path / "genes.gtf"
        gtf.write_text(
            '# GTF file\n'
            'chr17\tENSEMBL\tgene\t7565097\t7590856\t.\t-\t.\t'
            'gene_id "ENSG00000141510"; gene_name "TP53"; gene_type "protein_coding";\n'
            'chr17\tENSEMBL\texon\t7565097\t7565332\t.\t-\t.\t'
            'gene_id "ENSG00000141510"; gene_name "TP53"; gene_type "protein_coding";\n'
            'chr17\tENSEMBL\texon\t7572927\t7573008\t.\t-\t.\t'
            'gene_id "ENSG00000141510"; gene_name "TP53"; gene_type "protein_coding";\n'
        )
        annotator = SVAnnotator(gene_annotation_path=gtf)
        assert "chr17" in annotator._genes
        assert len(annotator._genes["chr17"]) == 1
        gene = annotator._genes["chr17"][0]
        assert gene.symbol == "TP53"
        assert gene.exon_count == 2

    def test_skips_non_protein_coding(self, tmp_path: Path) -> None:
        gtf = tmp_path / "genes.gtf"
        gtf.write_text(
            'chr1\tENSEMBL\tgene\t100\t500\t.\t+\t.\t'
            'gene_id "ENSG0001"; gene_name "lncRNA1"; gene_type "lncRNA";\n'
        )
        annotator = SVAnnotator(gene_annotation_path=gtf)
        assert annotator._genes.get("chr1", []) == []

    def test_missing_gtf_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            SVAnnotator(gene_annotation_path=tmp_path / "missing.gtf")


class TestDosageLoading:
    def test_loads_dosage_scores(self, tmp_path: Path) -> None:
        dosage_file = tmp_path / "dosage.tsv"
        dosage_file.write_text(
            "gene_symbol\thi_score\tts_score\n"
            "BRCA1\t3.0\t0.0\n"
            "PMP22\t3.0\t3.0\n"
            "UNKNOWN\tN/A\tN/A\n"
        )
        annotator = SVAnnotator(dosage_sensitivity_path=dosage_file)
        assert annotator._dosage["BRCA1"].hi_score == 3.0
        assert annotator._dosage["PMP22"].ts_score == 3.0
        assert annotator._dosage["UNKNOWN"].hi_score is None

    def test_missing_dosage_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            SVAnnotator(dosage_sensitivity_path=tmp_path / "missing.tsv")


class TestFrequencyDBLoading:
    def test_loads_sv_frequency_records(self, tmp_path: Path) -> None:
        freq_file = tmp_path / "gnomad_sv.tsv"
        freq_file.write_text(
            "# gnomAD-SV\n"
            "chr1\t1000\t5000\tDEL\t0.005\n"
            "chr1\t10000\t20000\tDUP\t0.01\n"
        )
        annotator = SVAnnotator(gnomad_sv_path=freq_file)
        assert len(annotator._sv_database["chr1"]) == 2
        # Sorted by start position
        assert annotator._sv_database["chr1"][0].start == 1000
        assert annotator._sv_database["chr1"][1].start == 10000

    def test_missing_frequency_db_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            SVAnnotator(gnomad_sv_path=tmp_path / "missing.bed")

    def test_skips_malformed_lines(self, tmp_path: Path) -> None:
        freq_file = tmp_path / "gnomad_sv.tsv"
        freq_file.write_text(
            "chr1\t1000\t5000\tDEL\t0.005\n"
            "chr1\tbad\t5000\tDEL\t0.005\n"
            "chr1\t2000\t3000\tINV\tnot_a_float\n"
            "short_line\n"
        )
        annotator = SVAnnotator(gnomad_sv_path=freq_file)
        # Only the first valid line should load
        assert len(annotator._sv_database.get("chr1", [])) == 1


class TestAnnotatorEndToEnd:
    def test_full_annotation_with_all_references(self, tmp_path: Path) -> None:
        gtf = tmp_path / "genes.gtf"
        gtf.write_text(
            'chr22\tENSEMBL\tgene\t19744226\t19771115\t.\t+\t.\t'
            'gene_id "ENSG0001"; gene_name "TBX1"; gene_type "protein_coding";\n'
            'chr22\tENSEMBL\texon\t19744226\t19750000\t.\t+\t.\t'
            'gene_id "ENSG0001"; gene_name "TBX1"; gene_type "protein_coding";\n'
        )
        dosage = tmp_path / "dosage.tsv"
        dosage.write_text("gene_symbol\thi_score\tts_score\nTBX1\t3.0\t0\n")
        freq = tmp_path / "gnomad.tsv"
        freq.write_text("chr22\t18900000\t21500000\tDEL\t0.0001\n")

        annotator = SVAnnotator(
            gene_annotation_path=gtf,
            dosage_sensitivity_path=dosage,
            gnomad_sv_path=freq,
            reciprocal_overlap=0.5,
            whole_gene_threshold=0.8,
        )

        sv = StructuralVariant(
            chrom="chr22", start=18916842, end=21465659,
            sv_type=SVType.DEL,
        )
        result = annotator._annotate_single(sv)

        assert result.consequence == SVConsequence.WHOLE_GENE_DELETION
        assert result.genes_affected == 1
        assert result.gene_overlaps[0].gene_symbol == "TBX1"
        assert result.gene_overlaps[0].is_haploinsufficient is True
        assert result.population_frequency == 0.0001
        assert result.frequency_unknown is False
