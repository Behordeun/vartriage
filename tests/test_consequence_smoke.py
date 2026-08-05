"""Smoke test for functional consequence assignment."""

import tempfile
from pathlib import Path

from vartriage._internal.interval_tree import (
    SortedArrayIntervalIndex,
    _parse_attributes,
)
from vartriage.annotation.consequence import ConsequenceAnnotator
from vartriage.models.variant import FunctionalConsequence, Variant

# Minimal GTF content for testing
SAMPLE_GTF = """\
##description: test annotation
chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";
chr1\ttest\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\texon\t1000\t1200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\tCDS\t1050\t1190\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\texon\t2000\t2500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\tCDS\t2000\t2500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\texon\t3000\t3500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
chr1\ttest\tCDS\t3000\t3500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";
"""


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


class TestAttributeParsing:
    def test_gtf_format(self):
        result = _parse_attributes('gene_id "BRCA1"; transcript_id "TX001";')
        assert result["gene_id"] == "BRCA1"
        assert result["transcript_id"] == "TX001"

    def test_gff3_format(self):
        result = _parse_attributes("gene_id=BRCA1;transcript_id=TX001")
        assert result["gene_id"] == "BRCA1"
        assert result["transcript_id"] == "TX001"


class TestSortedArrayIntervalIndex:
    def _create_index(self) -> tuple[SortedArrayIntervalIndex, Path]:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as tmp:
            tmp.write(SAMPLE_GTF)
            tmp.flush()
            tmp_name = tmp.name

        index = SortedArrayIntervalIndex()
        index.load(Path(tmp_name))
        return index, Path(tmp_name)

    def test_load_gtf(self):
        index, tmp_path = self._create_index()
        try:
            assert index._loaded is True
            assert "chr1" in index._chromosomes
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_overlap_coding_region(self):
        index, tmp_path = self._create_index()
        try:
            # Position 1100 is inside the CDS (1050-1190)
            hits = index.overlap("chr1", 1100, "A", "T")
            assert len(hits) > 0
            cds_hits = [h for h in hits if h["feature_type"] == "CDS"]
            assert len(cds_hits) > 0
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_no_overlap_intergenic(self):
        index, tmp_path = self._create_index()
        try:
            # Position 6000 is outside all gene regions
            hits = index.overlap("chr1", 6000, "A", "T")
            assert hits == []
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_splice_site_detection(self):
        index, tmp_path = self._create_index()
        try:
            hits = index.overlap("chr1", 1200, "A", "T")
            splice_hits = [h for h in hits if h.get("is_splice_site")]
            assert len(splice_hits) > 0
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_missing_file_raises(self):
        index = SortedArrayIntervalIndex()
        try:
            index.load(Path("/nonexistent/file.gtf"))
            raise AssertionError("Should have raised FileNotFoundError")
        except FileNotFoundError:
            pass

    def test_unknown_chromosome(self):
        index, tmp_path = self._create_index()
        try:
            hits = index.overlap("chrZ", 100, "A", "T")
            assert hits == []
        finally:
            tmp_path.unlink(missing_ok=True)


class TestConsequenceAnnotator:
    def _create_annotator(self) -> tuple[ConsequenceAnnotator, Path]:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as tmp:
            tmp.write(SAMPLE_GTF)
            tmp.flush()
            tmp_name = tmp.name
        return ConsequenceAnnotator(Path(tmp_name)), Path(tmp_name)

    def test_intergenic_variant(self):
        annotator, tmp_path = self._create_annotator()
        try:
            variant = _make_variant("chr1", 6000, "A", "T")
            result = annotator.assign(variant)
            assert result == FunctionalConsequence.INTERGENIC
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_missense_snv_in_cds(self):
        annotator, tmp_path = self._create_annotator()
        try:
            # Position 2100 is in CDS (2000-2500)
            variant = _make_variant("chr1", 2100, "A", "T")
            result = annotator.assign(variant)
            assert result == FunctionalConsequence.MISSENSE
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_frameshift_in_cds(self):
        annotator, tmp_path = self._create_annotator()
        try:
            # Insertion of 1 base (not divisible by 3) in CDS
            variant = _make_variant("chr1", 2100, "A", "AT")
            result = annotator.assign(variant)
            assert result == FunctionalConsequence.FRAMESHIFT
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_in_frame_insertion_in_cds(self):
        annotator, tmp_path = self._create_annotator()
        try:
            # Insertion of 3 bases (divisible by 3) in CDS
            variant = _make_variant("chr1", 2100, "A", "ATCG")
            result = annotator.assign(variant)
            assert result == FunctionalConsequence.IN_FRAME_INSERTION
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_in_frame_deletion_in_cds(self):
        annotator, tmp_path = self._create_annotator()
        try:
            # Deletion of 3 bases (divisible by 3) in CDS
            variant = _make_variant("chr1", 2100, "ATCG", "A")
            result = annotator.assign(variant)
            assert result == FunctionalConsequence.IN_FRAME_DELETION
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_severity_ranking_most_severe_wins(self):
        annotator, tmp_path = self._create_annotator()
        try:
            # A variant that overlaps both CDS and exon should get the CDS consequence
            variant = _make_variant("chr1", 2100, "A", "T")
            result = annotator.assign(variant)
            severity = {
                c: i
                for i, c in enumerate(
                    [
                        FunctionalConsequence.FRAMESHIFT,
                        FunctionalConsequence.NONSENSE,
                        FunctionalConsequence.SPLICE_SITE,
                        FunctionalConsequence.MISSENSE,
                        FunctionalConsequence.IN_FRAME_INSERTION,
                        FunctionalConsequence.IN_FRAME_DELETION,
                        FunctionalConsequence.SYNONYMOUS,
                        FunctionalConsequence.INTERGENIC,
                    ]
                )
            }
            assert severity[result] <= severity[FunctionalConsequence.MISSENSE]
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_assign_batch(self):
        annotator, tmp_path = self._create_annotator()
        try:
            variants = [
                _make_variant("chr1", 6000, "A", "T"),
                _make_variant("chr1", 2100, "A", "T"),
            ]
            results = annotator.assign_batch(variants)
            assert len(results) == 2
            assert results[0] == FunctionalConsequence.INTERGENIC
            assert results[1] == FunctionalConsequence.MISSENSE
        finally:
            tmp_path.unlink(missing_ok=True)
