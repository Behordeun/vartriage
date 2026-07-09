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
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".gtf", delete=False
        )
        tmp.write(SAMPLE_GTF)
        tmp.flush()
        tmp.close()

        index = SortedArrayIntervalIndex()
        index.load(Path(tmp.name))
        return index, Path(tmp.name)

    def test_load_gtf(self):
        index, _ = self._create_index()
        assert index._loaded is True
        assert "chr1" in index._chromosomes

    def test_overlap_coding_region(self):
        index, _ = self._create_index()
        # Position 1100 is inside the CDS (1050-1190)
        hits = index.overlap("chr1", 1100, "A", "T")
        assert len(hits) > 0
        cds_hits = [h for h in hits if h["feature_type"] == "CDS"]
        assert len(cds_hits) > 0

    def test_no_overlap_intergenic(self):
        index, _ = self._create_index()
        # Position 6000 is outside all gene regions
        hits = index.overlap("chr1", 6000, "A", "T")
        assert hits == []

    def test_splice_site_detection(self):
        index, _ = self._create_index()
        # Position 1199 is within 2 bases of exon end (1200)
        # The exon is chr1:1000-1200 (GTF 1-based), so 0-based: 999-1200
        # Exon end at 1200, donor site is 1198-1202
        # VCF pos 1200 -> 0-based start 1199, is in donor zone
        hits = index.overlap("chr1", 1200, "A", "T")
        splice_hits = [h for h in hits if h.get("is_splice_site")]
        assert len(splice_hits) > 0

    def test_missing_file_raises(self):
        index = SortedArrayIntervalIndex()
        try:
            index.load(Path("/nonexistent/file.gtf"))
            assert False, "Should have raised FileNotFoundError"
        except FileNotFoundError:
            pass

    def test_unknown_chromosome(self):
        index, _ = self._create_index()
        hits = index.overlap("chrZ", 100, "A", "T")
        assert hits == []


class TestConsequenceAnnotator:
    def _create_annotator(self) -> ConsequenceAnnotator:
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".gtf", delete=False
        )
        tmp.write(SAMPLE_GTF)
        tmp.flush()
        tmp.close()
        return ConsequenceAnnotator(Path(tmp.name))

    def test_intergenic_variant(self):
        annotator = self._create_annotator()
        variant = _make_variant("chr1", 6000, "A", "T")
        result = annotator.assign(variant)
        assert result == FunctionalConsequence.INTERGENIC

    def test_missense_snv_in_cds(self):
        annotator = self._create_annotator()
        # Position 2100 is in CDS (2000-2500)
        variant = _make_variant("chr1", 2100, "A", "T")
        result = annotator.assign(variant)
        assert result == FunctionalConsequence.MISSENSE

    def test_frameshift_in_cds(self):
        annotator = self._create_annotator()
        # Insertion of 1 base (not divisible by 3) in CDS
        variant = _make_variant("chr1", 2100, "A", "AT")
        result = annotator.assign(variant)
        assert result == FunctionalConsequence.FRAMESHIFT

    def test_in_frame_insertion_in_cds(self):
        annotator = self._create_annotator()
        # Insertion of 3 bases (divisible by 3) in CDS
        variant = _make_variant("chr1", 2100, "A", "ATCG")
        result = annotator.assign(variant)
        assert result == FunctionalConsequence.IN_FRAME_INSERTION

    def test_in_frame_deletion_in_cds(self):
        annotator = self._create_annotator()
        # Deletion of 3 bases (divisible by 3) in CDS
        variant = _make_variant("chr1", 2100, "ATCG", "A")
        result = annotator.assign(variant)
        assert result == FunctionalConsequence.IN_FRAME_DELETION

    def test_severity_ranking_most_severe_wins(self):
        annotator = self._create_annotator()
        # A variant that overlaps both CDS and exon should get the CDS consequence
        # since CDS gives Missense while exon gives Synonymous for an SNV
        variant = _make_variant("chr1", 2100, "A", "T")
        result = annotator.assign(variant)
        # Should be at least Missense (more severe than Synonymous)
        severity = {c: i for i, c in enumerate(
            [FunctionalConsequence.FRAMESHIFT,
             FunctionalConsequence.NONSENSE,
             FunctionalConsequence.SPLICE_SITE,
             FunctionalConsequence.MISSENSE,
             FunctionalConsequence.IN_FRAME_INSERTION,
             FunctionalConsequence.IN_FRAME_DELETION,
             FunctionalConsequence.SYNONYMOUS,
             FunctionalConsequence.INTERGENIC]
        )}
        assert severity[result] <= severity[FunctionalConsequence.MISSENSE]

    def test_assign_batch(self):
        annotator = self._create_annotator()
        variants = [
            _make_variant("chr1", 6000, "A", "T"),
            _make_variant("chr1", 2100, "A", "T"),
        ]
        results = annotator.assign_batch(variants)
        assert len(results) == 2
        assert results[0] == FunctionalConsequence.INTERGENIC
        assert results[1] == FunctionalConsequence.MISSENSE
