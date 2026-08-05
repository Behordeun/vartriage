"""ClinGen-based structural variant classification.

Implements the ACMG/ClinGen Technical Standards for interpretation of
copy-number losses and gains (Riggs et al. 2020). Evaluates scored SVs
against evidence categories from Sections 1-4 and produces a 5-tier
classification (Pathogenic through Benign).

The point-based scoring system maps accumulated evidence to final
classification using ClinGen-defined thresholds:
  >= 0.99  Pathogenic
  0.90-0.98  Likely Pathogenic
  -0.89 to 0.89  VUS (uncertain)
  -0.98 to -0.90  Likely Benign
  <= -0.99  Benign
"""

from __future__ import annotations

import logging
from collections.abc import Iterator

from vartriage.structural.models import (
    AnnotatedSV,
    ClassifiedSV,
    ScoredSV,
    SVClassification,
    SVEvidenceCategory,
    SVType,
)

logger = logging.getLogger(__name__)

# Classification thresholds from ClinGen scoring framework
_PATHOGENIC_THRESHOLD: float = 0.99
_LIKELY_PATHOGENIC_THRESHOLD: float = 0.90
_LIKELY_BENIGN_THRESHOLD: float = -0.90
_BENIGN_THRESHOLD: float = -0.99


class SVClassifier:
    """Classify scored structural variants using ClinGen evidence framework.

    Evaluates each ScoredSV against applicable evidence categories,
    accumulates points, and maps the total to a 5-tier classification.

    The classifier handles losses (DEL, CNV with CN<2) and gains
    (DUP, CNV with CN>2) separately, applying the appropriate
    evidence sections for each.

    Parameters
    ----------
    pathogenic_regions : list[tuple[str, int, int]]
        Known pathogenic CNV regions as (chrom, start, end) tuples.
        Used for Section 2 overlap evaluation.
    benign_regions : list[tuple[str, int, int]]
        Known benign CNV regions as (chrom, start, end) tuples.
        Used for Section 2 benign overlap evaluation.
    """

    def __init__(
        self,
        pathogenic_regions: list[tuple[str, int, int]] | None = None,
        benign_regions: list[tuple[str, int, int]] | None = None,
        pathogenic_region_names: dict[tuple[str, int, int], str] | None = None,
    ) -> None:
        self._pathogenic_regions = pathogenic_regions or []
        self._benign_regions = benign_regions or []
        self._pathogenic_region_names = pathogenic_region_names or {}

    def classify(self, variants: Iterator[ScoredSV]) -> Iterator[ClassifiedSV]:
        """Classify a stream of scored SVs.

        Parameters
        ----------
        variants : Iterator[ScoredSV]
            Scored SVs from the scoring engine.

        Yields
        ------
        ClassifiedSV
            Each SV with evidence categories, accumulated score,
            and final classification.
        """
        for scored in variants:
            yield self._classify_single(scored)

    def _classify_single(self, scored: ScoredSV) -> ClassifiedSV:
        """Evaluate all applicable evidence categories for one SV."""
        categories: set[SVEvidenceCategory] = set()
        missing_sources: set[str] = set()
        evidence_score = 0.0
        syndrome_name: str | None = None

        annotated = scored.annotated
        is_loss = self._is_loss(annotated)

        # Section 1: Initial assessment of genomic content
        evidence_score += self._evaluate_section1(
            annotated, categories, missing_sources
        )

        # Section 2: Overlap with established regions (also resolves syndrome)
        s2_score, matched_syndrome = self._evaluate_section2_with_syndrome(
            annotated, is_loss, categories, missing_sources
        )
        evidence_score += s2_score
        if matched_syndrome is not None:
            syndrome_name = matched_syndrome

        # Section 3: Gene-level evaluation (losses)
        if is_loss:
            evidence_score += self._evaluate_section3(
                annotated, categories, missing_sources
            )

        # Section 4: Duplication-specific evaluation (gains)
        if not is_loss and annotated.sv.sv_type in (SVType.DUP, SVType.CNV):
            evidence_score += self._evaluate_section4(
                annotated, categories, missing_sources
            )

        # Frequency-based benign evidence
        evidence_score += self._evaluate_frequency_evidence(annotated, categories)

        classification = self._score_to_classification(evidence_score)

        return ClassifiedSV(
            scored=scored,
            classification=classification,
            evidence_categories=frozenset(categories),
            evidence_score=evidence_score,
            missing_data_sources=frozenset(missing_sources),
            syndrome_name=syndrome_name,
        )

    def _evaluate_section1(
        self,
        sv: AnnotatedSV,
        categories: set[SVEvidenceCategory],
        missing_sources: set[str],
    ) -> float:
        """Section 1: Initial assessment of genomic content.

        1A: Contains protein-coding genes (informational, 0 points)
        1B: Contains established HI/TS gene (strong positive evidence)
        """
        score = 0.0

        if sv.genes_affected > 0:
            categories.add(SVEvidenceCategory.CONTAINS_PROTEIN_CODING)

        if sv.hi_genes_affected > 0:
            categories.add(SVEvidenceCategory.CONTAINS_ESTABLISHED_HI_GENE)
            # Strong evidence when established HI gene is fully contained
            score += 0.45

        return score

    def _evaluate_section2(
        self,
        sv: AnnotatedSV,
        is_loss: bool,
        categories: set[SVEvidenceCategory],
        missing_sources: set[str],
    ) -> float:
        """Section 2: Overlap with established pathogenic/benign CNV regions.

        Compares the SV span against known pathogenic and benign regions.
        """
        score = 0.0
        sv_chrom = sv.sv.chrom
        sv_start = sv.sv.start
        sv_end = sv.sv.end

        # Check pathogenic region overlap
        if not self._pathogenic_regions:
            missing_sources.add("pathogenic_regions")
        else:
            path_overlap = self._best_region_overlap(
                sv_chrom, sv_start, sv_end, self._pathogenic_regions
            )
            if path_overlap is not None:
                frac_query, frac_ref = path_overlap
                if frac_ref >= 0.9:
                    # SV completely overlaps a pathogenic region
                    categories.add(SVEvidenceCategory.COMPLETE_OVERLAP_PATHOGENIC)
                    score += 0.45
                elif frac_ref >= 0.5:
                    # Partial overlap with pathogenic region
                    categories.add(SVEvidenceCategory.PARTIAL_OVERLAP_PATHOGENIC)
                    score += 0.25
                elif frac_query >= 0.5:
                    # SV is smaller than pathogenic region but overlaps it
                    categories.add(SVEvidenceCategory.OVERLAP_SMALLER_THAN_PATHOGENIC)

        # Check benign region overlap
        if not self._benign_regions:
            missing_sources.add("benign_regions")
        else:
            benign_overlap = self._best_region_overlap(
                sv_chrom, sv_start, sv_end, self._benign_regions
            )
            if benign_overlap is not None:
                frac_query, frac_ref = benign_overlap
                if frac_query >= 0.9:
                    # SV is fully contained within a benign region
                    categories.add(SVEvidenceCategory.CONTAINED_WITHIN_BENIGN)
                    score -= 0.60
                elif frac_query >= 0.5:
                    # Partial overlap with benign region
                    categories.add(SVEvidenceCategory.PARTIAL_OVERLAP_BENIGN)
                    score -= 0.30
                elif frac_ref >= 0.9:
                    # SV completely contains a benign region (larger than it)
                    categories.add(SVEvidenceCategory.COMPLETELY_CONTAINS_BENIGN)
                    # No score change: containing a benign region doesn't
                    # make the larger SV benign

        return score

    def _evaluate_section2_with_syndrome(
        self,
        sv: AnnotatedSV,
        is_loss: bool,
        categories: set[SVEvidenceCategory],
        missing_sources: set[str],
    ) -> tuple[float, str | None]:
        """Section 2 with syndrome name resolution.

        Same logic as _evaluate_section2 but also returns the matched
        syndrome name from the pathogenic region lookup.
        """
        score = 0.0
        syndrome: str | None = None
        sv_chrom = sv.sv.chrom
        sv_start = sv.sv.start
        sv_end = sv.sv.end

        if not self._pathogenic_regions:
            missing_sources.add("pathogenic_regions")
        else:
            best_region, path_overlap = self._best_region_overlap_with_key(
                sv_chrom, sv_start, sv_end, self._pathogenic_regions
            )
            if path_overlap is not None:
                frac_query, frac_ref = path_overlap
                if frac_ref >= 0.9:
                    categories.add(SVEvidenceCategory.COMPLETE_OVERLAP_PATHOGENIC)
                    score += 0.45
                elif frac_ref >= 0.5:
                    categories.add(SVEvidenceCategory.PARTIAL_OVERLAP_PATHOGENIC)
                    score += 0.25
                elif frac_query >= 0.5:
                    categories.add(SVEvidenceCategory.OVERLAP_SMALLER_THAN_PATHOGENIC)

                # Resolve syndrome name from the matched region
                if best_region is not None and frac_ref >= 0.5:
                    syndrome = self._pathogenic_region_names.get(best_region)

        if not self._benign_regions:
            missing_sources.add("benign_regions")
        else:
            benign_overlap = self._best_region_overlap(
                sv_chrom, sv_start, sv_end, self._benign_regions
            )
            if benign_overlap is not None:
                frac_query, frac_ref = benign_overlap
                if frac_query >= 0.9:
                    categories.add(SVEvidenceCategory.CONTAINED_WITHIN_BENIGN)
                    score -= 0.60
                elif frac_query >= 0.5:
                    categories.add(SVEvidenceCategory.PARTIAL_OVERLAP_BENIGN)
                    score -= 0.30
                elif frac_ref >= 0.9:
                    categories.add(SVEvidenceCategory.COMPLETELY_CONTAINS_BENIGN)

        return score, syndrome

    def _evaluate_section3(
        self,
        sv: AnnotatedSV,
        categories: set[SVEvidenceCategory],
        missing_sources: set[str],
    ) -> float:
        """Section 3: Gene-level evaluation for losses (deletions).

        3A: Established HI gene fully contained in the SV
        3B: Gene partially deleted (one breakpoint within gene)
        3C: Breakpoint within gene but functional impact uncertain
        """
        score = 0.0

        if not sv.gene_overlaps:
            return score

        for overlap in sv.gene_overlaps:
            if overlap.is_whole_gene and overlap.is_haploinsufficient:
                # 3A: HI gene fully deleted - very strong evidence
                categories.add(SVEvidenceCategory.GENE_FULLY_CONTAINED)
                score += 0.45
                break
            elif not overlap.is_whole_gene and overlap.overlap_fraction >= 0.1:
                if overlap.is_haploinsufficient:
                    # 3B: HI gene partially deleted
                    categories.add(SVEvidenceCategory.GENE_PARTIALLY_DELETED)
                    score += 0.25
                    break
                else:
                    # 3C: Non-HI gene has breakpoint within it
                    categories.add(SVEvidenceCategory.BREAKPOINT_WITHIN_GENE)
                    score += 0.10
                    break

        return score

    def _evaluate_section4(
        self,
        sv: AnnotatedSV,
        categories: set[SVEvidenceCategory],
        missing_sources: set[str],
    ) -> float:
        """Section 4: Duplication-specific evaluation.

        4F: Established TS gene fully contained (tandem dup assumed)
        4G: Gene disrupted by duplication breakpoint
        4H: Intragenic duplication without predicted disruption
        """
        score = 0.0

        if not sv.gene_overlaps:
            return score

        # Check pathogenic region overlap for duplications
        if self._pathogenic_regions:
            path_overlap = self._best_region_overlap(
                sv.sv.chrom, sv.sv.start, sv.sv.end, self._pathogenic_regions
            )
            if path_overlap is not None:
                frac_query, frac_ref = path_overlap
                if frac_query >= 0.9 and frac_ref >= 0.9:
                    categories.add(SVEvidenceCategory.DUP_IDENTICAL_TO_PATHOGENIC)
                    score += 0.45
                    return score
                elif frac_ref >= 0.9:
                    categories.add(SVEvidenceCategory.DUP_COMPLETE_OVERLAP_PATHOGENIC)
                    score += 0.35
                    return score

        for overlap in sv.gene_overlaps:
            if overlap.is_whole_gene and overlap.is_triplosensitive:
                # 4F: TS gene fully contained in duplication
                categories.add(SVEvidenceCategory.DUP_TS_GENE_CONTAINED)
                score += 0.35
                break
            elif not overlap.is_whole_gene and overlap.overlap_fraction >= 0.1:
                # Breakpoint within a gene - potential disruption
                if overlap.exons_affected < overlap.total_exons:
                    # 4G: Partial dup likely disrupts the gene
                    categories.add(SVEvidenceCategory.DUP_GENE_DISRUPTED)
                    score += 0.20
                    break
                else:
                    # 4H: Whole gene dup without clear disruption
                    categories.add(SVEvidenceCategory.DUP_INTRAGENIC_NO_DISRUPTION)
                    break

        return score

    def _evaluate_frequency_evidence(
        self,
        sv: AnnotatedSV,
        categories: set[SVEvidenceCategory],
    ) -> float:
        """Apply frequency-based evidence adjustment.

        Common SVs (present in population databases) get negative
        evidence points, pushing toward benign classification.
        """
        if sv.population_frequency is None:
            return 0.0

        af = sv.population_frequency

        # Very common SVs that somehow passed filtering get strong
        # benign evidence
        if af >= 0.01:
            return -0.60
        if af >= 0.005:
            return -0.30
        if af >= 0.001:
            return -0.15

        return 0.0

    def _score_to_classification(self, score: float) -> SVClassification:
        """Map accumulated evidence score to 5-tier classification."""
        if score >= _PATHOGENIC_THRESHOLD:
            return SVClassification.PATHOGENIC
        if score >= _LIKELY_PATHOGENIC_THRESHOLD:
            return SVClassification.LIKELY_PATHOGENIC
        if score <= _BENIGN_THRESHOLD:
            return SVClassification.BENIGN
        if score <= _LIKELY_BENIGN_THRESHOLD:
            return SVClassification.LIKELY_BENIGN
        return SVClassification.VUS

    def _is_loss(self, sv: AnnotatedSV) -> bool:
        """Determine if the SV represents a copy-number loss."""
        if sv.sv.sv_type == SVType.DEL:
            return True
        if sv.sv.sv_type == SVType.CNV:
            return sv.sv.copy_number is not None and sv.sv.copy_number < 2
        return False

    def _best_region_overlap(
        self,
        chrom: str,
        start: int,
        end: int,
        regions: list[tuple[str, int, int]],
    ) -> tuple[float, float] | None:
        """Find best overlapping region and return overlap fractions.

        Returns (fraction_of_query_covered, fraction_of_reference_covered)
        for the best-matching region, or None if no overlap found.
        """
        best: tuple[float, float] | None = None
        best_total = 0.0
        sv_length = end - start + 1

        for r_chrom, r_start, r_end in regions:
            if r_chrom != chrom:
                # Try chr prefix normalization
                alt_chrom = (
                    chrom.replace("chr", "")
                    if chrom.startswith("chr")
                    else f"chr{chrom}"
                )
                if r_chrom != alt_chrom:
                    continue

            overlap_start = max(start, r_start)
            overlap_end = min(end, r_end)

            if overlap_start > overlap_end:
                continue

            overlap_bp = overlap_end - overlap_start + 1
            ref_length = r_end - r_start + 1

            frac_query = overlap_bp / sv_length if sv_length > 0 else 0.0
            frac_ref = overlap_bp / ref_length if ref_length > 0 else 0.0
            total = frac_query + frac_ref

            if total > best_total:
                best_total = total
                best = (frac_query, frac_ref)

        return best

    def _best_region_overlap_with_key(
        self,
        chrom: str,
        start: int,
        end: int,
        regions: list[tuple[str, int, int]],
    ) -> tuple[tuple[str, int, int] | None, tuple[float, float] | None]:
        """Like _best_region_overlap but also returns the matched region tuple.

        Returns (matched_region_key, (frac_query, frac_ref)) or (None, None).
        """
        best_fracs: tuple[float, float] | None = None
        best_total = 0.0
        best_key: tuple[str, int, int] | None = None
        sv_length = end - start + 1

        for r_chrom, r_start, r_end in regions:
            if r_chrom != chrom:
                alt_chrom = (
                    chrom.replace("chr", "")
                    if chrom.startswith("chr")
                    else f"chr{chrom}"
                )
                if r_chrom != alt_chrom:
                    continue

            overlap_start = max(start, r_start)
            overlap_end = min(end, r_end)

            if overlap_start > overlap_end:
                continue

            overlap_bp = overlap_end - overlap_start + 1
            ref_length = r_end - r_start + 1

            frac_query = overlap_bp / sv_length if sv_length > 0 else 0.0
            frac_ref = overlap_bp / ref_length if ref_length > 0 else 0.0
            total = frac_query + frac_ref

            if total > best_total:
                best_total = total
                best_fracs = (frac_query, frac_ref)
                best_key = (r_chrom, r_start, r_end)

        return best_key, best_fracs
