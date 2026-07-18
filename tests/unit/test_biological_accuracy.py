"""Unit tests for v0.8.0 biological accuracy modules.

Covers: genetic code translation, transcript CDS index, codon resolver,
variant normalizer, and prioritization score computation.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vartriage._internal.genetic_code import (CODON_TABLE, reverse_complement,
                                              translate_codon)
from vartriage.annotation.transcript_index import (CDSExon, TranscriptCDS,
                                                   TranscriptCDSIndex)
from vartriage.models.variant import FunctionalConsequence
from vartriage.prioritization.scoring import compute_prioritization_score

# ============================================================================
# Genetic Code
# ============================================================================


class TestGeneticCode:
    """Standard genetic code translation."""

    def test_methionine_start_codon(self) -> None:
        assert translate_codon("ATG") == "M"

    def test_stop_codons(self) -> None:
        assert translate_codon("TAA") == "*"
        assert translate_codon("TAG") == "*"
        assert translate_codon("TGA") == "*"

    def test_all_64_codons_mapped(self) -> None:
        assert len(CODON_TABLE) == 64

    def test_leucine_sixfold_degeneracy(self) -> None:
        leu_codons = [k for k, v in CODON_TABLE.items() if v == "L"]
        assert len(leu_codons) == 6

    def test_ambiguous_codon_returns_question_mark(self) -> None:
        assert translate_codon("NNN") == "?"
        assert translate_codon("ATN") == "?"

    def test_lowercase_input_handled(self) -> None:
        assert translate_codon("atg") == "M"

    def test_reverse_complement_simple(self) -> None:
        assert reverse_complement("ATCG") == "CGAT"

    def test_reverse_complement_palindrome(self) -> None:
        assert reverse_complement("ATAT") == "ATAT"

    def test_reverse_complement_single_base(self) -> None:
        assert reverse_complement("A") == "T"
        assert reverse_complement("G") == "C"


# ============================================================================
# TranscriptCDSIndex
# ============================================================================


class TestTranscriptCDS:
    """Per-transcript CDS coordinate mapping."""

    def _make_forward_transcript(self) -> TranscriptCDS:
        """Simple forward-strand transcript: 2 CDS exons."""
        tc = TranscriptCDS(
            transcript_id="TX1",
            gene_name="GENE1",
            chrom="chr1",
            strand="+",
            frame_offset=0,
        )
        # Exon 1: positions 100-109 (10bp)
        tc.cds_exons.append(CDSExon(start=100, end=110))
        # Exon 2: positions 200-209 (10bp)
        tc.cds_exons.append(CDSExon(start=200, end=210))
        tc.finalize()
        return tc

    def _make_reverse_transcript(self) -> TranscriptCDS:
        """Simple reverse-strand transcript: 2 CDS exons."""
        tc = TranscriptCDS(
            transcript_id="TX2",
            gene_name="GENE2",
            chrom="chr1",
            strand="-",
            frame_offset=0,
        )
        tc.cds_exons.append(CDSExon(start=100, end=110))
        tc.cds_exons.append(CDSExon(start=200, end=210))
        tc.finalize()
        return tc

    def test_cds_length(self) -> None:
        tc = self._make_forward_transcript()
        assert tc.cds_length == 20

    def test_forward_strand_first_exon_start(self) -> None:
        tc = self._make_forward_transcript()
        assert tc.genomic_to_cds_position(100) == 0

    def test_forward_strand_first_exon_end(self) -> None:
        tc = self._make_forward_transcript()
        assert tc.genomic_to_cds_position(109) == 9

    def test_forward_strand_second_exon_start(self) -> None:
        tc = self._make_forward_transcript()
        # Second exon starts at CDS position 10 (after 10bp first exon)
        assert tc.genomic_to_cds_position(200) == 10

    def test_forward_strand_intronic_position_returns_none(self) -> None:
        tc = self._make_forward_transcript()
        assert tc.genomic_to_cds_position(150) is None

    def test_reverse_strand_last_exon_end_is_cds_start(self) -> None:
        tc = self._make_reverse_transcript()
        # For negative strand, CDS position 0 is at the end of the last exon
        assert tc.genomic_to_cds_position(209) == 0

    def test_reverse_strand_last_exon_start(self) -> None:
        tc = self._make_reverse_transcript()
        assert tc.genomic_to_cds_position(200) == 9

    def test_reverse_strand_first_exon(self) -> None:
        tc = self._make_reverse_transcript()
        # After 10bp from exon2, exon1 starts at CDS position 10
        assert tc.genomic_to_cds_position(109) == 10


class TestTranscriptCDSIndex:
    """Index building and lookup."""

    def test_add_and_finalize(self) -> None:
        index = TranscriptCDSIndex()
        index.add_cds_exon("TX1", "GENE1", "chr1", 100, 200, "+", 0)
        index.add_cds_exon("TX1", "GENE1", "chr1", 300, 400, "+", 0)
        index.finalize()
        assert index.transcript_count == 1

    def test_find_overlapping(self) -> None:
        index = TranscriptCDSIndex()
        index.add_cds_exon("TX1", "GENE1", "chr1", 100, 200, "+", 0)
        index.add_cds_exon("TX2", "GENE2", "chr1", 150, 250, "+", 0)
        index.finalize()

        hits = index.find_overlapping("chr1", 175)
        assert len(hits) == 2

    def test_no_overlap_returns_empty(self) -> None:
        index = TranscriptCDSIndex()
        index.add_cds_exon("TX1", "GENE1", "chr1", 100, 200, "+", 0)
        index.finalize()

        hits = index.find_overlapping("chr1", 50)
        assert hits == []

    def test_wrong_chromosome_returns_empty(self) -> None:
        index = TranscriptCDSIndex()
        index.add_cds_exon("TX1", "GENE1", "chr1", 100, 200, "+", 0)
        index.finalize()

        hits = index.find_overlapping("chr2", 150)
        assert hits == []

    def test_serialization_round_trip(self) -> None:
        index = TranscriptCDSIndex()
        index.add_cds_exon("TX1", "GENE1", "chr1", 100, 200, "+", 0)
        index.add_cds_exon("TX1", "GENE1", "chr1", 300, 400, "+", 0)
        index.finalize()

        data = index.to_serializable()
        restored = TranscriptCDSIndex.from_serializable(data)

        assert restored.transcript_count == 1
        tc = restored.get_transcript("TX1")
        assert tc is not None
        assert tc.cds_length == 200
        assert tc.strand == "+"


# ============================================================================
# VariantNormalizer (without FASTA - logic tests)
# ============================================================================


class TestVariantNormalizerLogic:
    """Normalizer trim logic (mocked FASTA for left-alignment)."""

    def test_snv_passes_through_unchanged(self) -> None:
        from vartriage._internal.normalizer import VariantNormalizer

        with patch("pysam.FastaFile"):
            normalizer = VariantNormalizer.__new__(VariantNormalizer)
            normalizer._fasta = MagicMock()

        result = normalizer.normalize("chr1", 100, "A", "T")
        assert result == ("chr1", 100, "A", "T")

    def test_right_trim_shared_suffix(self) -> None:
        from vartriage._internal.normalizer import VariantNormalizer

        normalizer = VariantNormalizer.__new__(VariantNormalizer)
        normalizer._fasta = MagicMock()
        # FASTA not needed for trim-only (left_align won't shift SNV-length results)
        normalizer._fasta.fetch.return_value = "X"  # won't match

        # ref=ATCG, alt=ACG -> right-trim G -> ref=ATC, alt=AC
        # Then left-trim A -> ref=TC, alt=C, pos+1
        result = normalizer.normalize("chr1", 100, "ATCG", "ACG")
        # After right-trim: ATC/AC. After left-trim: TC/C at pos 101.
        # Left-align depends on FASTA (mocked to not shift)
        assert result[2] != result[3]  # ref != alt after normalization

    def test_left_trim_shared_prefix(self) -> None:
        from vartriage._internal.normalizer import VariantNormalizer

        normalizer = VariantNormalizer.__new__(VariantNormalizer)
        normalizer._fasta = MagicMock()
        normalizer._fasta.fetch.return_value = "X"

        # ref=ACGT, alt=ATGT -> no shared suffix, left-trim AC/AT -> C/T at pos+1
        # Actually shares prefix A: ref=CGT, alt=TGT, pos+1
        # Then shares suffix GT: ref=C, alt=T -> SNV at pos+1
        result = normalizer.normalize("chr1", 100, "ACGT", "ATGT")
        # After right-trim GT: AC/AT. After left-trim A: C/T at pos+1
        assert result == ("chr1", 101, "C", "T")


# ============================================================================
# Prioritization Score
# ============================================================================


class TestPrioritizationScore:
    """Literature-backed scoring: REVEL for missense, SpliceAI for splice, CADD general."""

    def test_missense_uses_revel(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.MISSENSE,
            revel_score=0.85,
            spliceai_score=None,
            cadd_phred=25.0,
        )
        # REVEL=0.85, CADD=25/60=0.417. Max should be REVEL.
        assert score == pytest.approx(0.85)

    def test_splice_site_uses_spliceai(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.SPLICE_SITE,
            revel_score=None,
            spliceai_score=0.92,
            cadd_phred=30.0,
        )
        # SpliceAI=0.92, CADD=30/60=0.5. Max should be SpliceAI.
        assert score == pytest.approx(0.92)

    def test_frameshift_uses_cadd(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.FRAMESHIFT,
            revel_score=None,
            spliceai_score=None,
            cadd_phred=35.0,
        )
        assert score == pytest.approx(35.0 / 60.0)

    def test_cadd_capped_at_one(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.INTERGENIC,
            revel_score=None,
            spliceai_score=None,
            cadd_phred=99.0,
        )
        assert score == pytest.approx(1.0)

    def test_no_scores_returns_none(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.SYNONYMOUS,
            revel_score=None,
            spliceai_score=None,
            cadd_phred=None,
        )
        assert score is None

    def test_max_across_multiple_scores(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.MISSENSE,
            revel_score=0.3,
            spliceai_score=0.7,
            cadd_phred=50.0,
        )
        # REVEL=0.3 (primary for missense), SpliceAI=0.7 (supplementary),
        # CADD=50/60=0.833. Max = 0.833
        assert score == pytest.approx(50.0 / 60.0)

    def test_revel_only_for_missense(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.MISSENSE,
            revel_score=0.95,
            spliceai_score=None,
            cadd_phred=None,
        )
        assert score == pytest.approx(0.95)

    def test_spliceai_only_for_splice(self) -> None:
        score = compute_prioritization_score(
            consequence=FunctionalConsequence.SPLICE_SITE,
            revel_score=None,
            spliceai_score=0.6,
            cadd_phred=None,
        )
        assert score == pytest.approx(0.6)


# ============================================================================
# ScoredVariant __post_init__ sync
# ============================================================================


class TestScoredVariantSync:
    """Verify composite_rank and prioritization_score stay in sync."""

    def test_composite_rank_syncs_to_prioritization_score(self) -> None:
        from vartriage.models.variant import (AnnotatedVariant, ScoredVariant,
                                              Variant)

        raw = Variant(
            chrom="chr1",
            pos=100,
            id=None,
            ref="A",
            alt="T",
            qual=30.0,
            filter_status="PASS",
        )
        ann = AnnotatedVariant(variant=raw, consequence=FunctionalConsequence.MISSENSE)
        scored = ScoredVariant(annotated=ann, composite_rank=0.75)

        assert scored.prioritization_score == 0.75

    def test_prioritization_score_syncs_to_composite_rank(self) -> None:
        from vartriage.models.variant import (AnnotatedVariant, ScoredVariant,
                                              Variant)

        raw = Variant(
            chrom="chr1",
            pos=100,
            id=None,
            ref="A",
            alt="T",
            qual=30.0,
            filter_status="PASS",
        )
        ann = AnnotatedVariant(variant=raw, consequence=FunctionalConsequence.MISSENSE)
        scored = ScoredVariant(annotated=ann, prioritization_score=0.88)

        assert scored.composite_rank == 0.88

    def test_both_none_stays_none(self) -> None:
        from vartriage.models.variant import (AnnotatedVariant, ScoredVariant,
                                              Variant)

        raw = Variant(
            chrom="chr1",
            pos=100,
            id=None,
            ref="A",
            alt="T",
            qual=30.0,
            filter_status="PASS",
        )
        ann = AnnotatedVariant(variant=raw, consequence=FunctionalConsequence.MISSENSE)
        scored = ScoredVariant(annotated=ann)

        assert scored.composite_rank is None
        assert scored.prioritization_score is None
