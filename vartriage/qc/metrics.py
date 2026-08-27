"""QC metrics dataclass and streaming computation engine."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pysam

_TRANSITIONS = frozenset({("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")})


def _is_transition(ref: str, alt: str) -> bool:
    return (ref.upper(), alt.upper()) in _TRANSITIONS


def _is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref != alt


@dataclass(frozen=True)
class QCMetrics:
    """Computed QC metrics from a single VCF pass.

    All counts are based on biallelic decomposition. Multi-allelic sites
    contribute one count per ALT allele.
    """

    total_variants: int
    snv_count: int
    indel_count: int
    insertion_count: int
    deletion_count: int
    transition_count: int
    transversion_count: int
    ti_tv_ratio: float
    het_count: int | None
    hom_alt_count: int | None
    het_hom_ratio: float | None
    ins_del_ratio: float
    per_chrom_counts: dict[str, int]


class QCComputer:
    """Single-pass QC metric computation from a VCF stream.

    Iterates pysam records once, counting transitions, transversions,
    het/hom genotypes, insertions/deletions, and per-chromosome variant
    counts. O(N) time, O(chromosomes) memory.

    Parameters
    ----------
    sample_id : str | None
        Sample name for genotype extraction. When None, het/hom metrics
        are skipped (reported as None in QCMetrics).
    """

    def __init__(self, sample_id: str | None = None) -> None:
        self._sample_id = sample_id
        self._total = 0
        self._snv_count = 0
        self._transition_count = 0
        self._transversion_count = 0
        self._insertion_count = 0
        self._deletion_count = 0
        self._het_count = 0
        self._hom_alt_count = 0
        self._chrom_counts: dict[str, int] = {}
        self._has_sample = sample_id is not None

    def process_record(self, record: pysam.VariantRecord) -> None:
        """Process a single VCF record, updating internal counters.

        Multi-allelic sites are decomposed: each valid ALT allele
        contributes its own count to variant totals, Ti/Tv, and
        indel metrics. Genotype het/hom counting stays per-record,
        since a sample has one genotype call regardless of ALT count.
        """
        chrom = record.contig
        ref = record.ref
        alts = record.alts
        if not alts:
            self._chrom_counts[chrom] = self._chrom_counts.get(chrom, 0) + 1
            self._total += 1
            return

        for alt in alts:
            if alt is None or alt == "." or ref is None:
                continue
            self._chrom_counts[chrom] = self._chrom_counts.get(chrom, 0) + 1
            self._total += 1
            self._classify_allele(ref, alt)

        if self._has_sample:
            self._extract_genotype(record)

    def _classify_allele(self, ref: str, alt: str) -> None:
        """Update SNV/indel counters for one decomposed ALT allele."""
        if _is_snv(ref, alt):
            self._snv_count += 1
            if _is_transition(ref, alt):
                self._transition_count += 1
            else:
                self._transversion_count += 1
            return

        # Indel classification by length comparison
        if len(alt) > len(ref):
            self._insertion_count += 1
        elif len(ref) > len(alt):
            self._deletion_count += 1

    def _extract_genotype(self, record: pysam.VariantRecord) -> None:
        """Count het and hom-alt genotypes for the target sample."""
        if self._sample_id is None:
            return
        try:
            sample = record.samples[self._sample_id]
            gt = sample.get("GT")
            if gt is None:
                return
            alleles = gt
            # Filter out missing calls (None alleles)
            non_missing = [a for a in alleles if a is not None]
            if len(non_missing) < 2:
                return
            # het: different alleles, at least one non-ref
            # hom-alt: all alleles non-ref and identical
            if non_missing[0] != non_missing[1]:
                if any(a != 0 for a in non_missing):
                    self._het_count += 1
            elif non_missing[0] != 0:
                self._hom_alt_count += 1
        except (KeyError, TypeError, IndexError):
            pass

    def finalize(self) -> QCMetrics:
        """Compute derived metrics and return frozen QCMetrics."""
        ti_tv = (
            self._transition_count / self._transversion_count
            if self._transversion_count > 0
            else 0.0
        )

        ins_del = (
            self._insertion_count / self._deletion_count
            if self._deletion_count > 0
            else 0.0
        )

        het_count: int | None = None
        hom_alt_count: int | None = None
        het_hom_ratio: float | None = None

        if self._has_sample:
            het_count = self._het_count
            hom_alt_count = self._hom_alt_count
            het_hom_ratio = (
                self._het_count / self._hom_alt_count
                if self._hom_alt_count > 0
                else 0.0
            )

        return QCMetrics(
            total_variants=self._total,
            snv_count=self._snv_count,
            indel_count=self._insertion_count + self._deletion_count,
            insertion_count=self._insertion_count,
            deletion_count=self._deletion_count,
            transition_count=self._transition_count,
            transversion_count=self._transversion_count,
            ti_tv_ratio=ti_tv,
            het_count=het_count,
            hom_alt_count=hom_alt_count,
            het_hom_ratio=het_hom_ratio,
            ins_del_ratio=ins_del,
            per_chrom_counts=dict(self._chrom_counts),
        )


def compute_qc_metrics(
    vcf_path: Any,
    sample_id: str | None = None,
) -> QCMetrics:
    """Run a full QC pass over a VCF file and return metrics.

    Opens the VCF independently for the QC pass (separate from the
    annotation pass). Streaming: O(N) time, O(chromosomes) memory.

    Parameters
    ----------
    vcf_path : Path
        Path to a .vcf or .vcf.gz file.
    sample_id : str | None
        Sample to extract genotypes from. None skips het/hom.

    Returns
    -------
    QCMetrics
        Frozen dataclass with all computed metrics.
    """
    vcf = pysam.VariantFile(str(vcf_path), "r")

    # Auto-detect the sample when none is given and the VCF has exactly one.
    effective_sample = sample_id
    if effective_sample is None:
        samples = list(vcf.header.samples)
        if len(samples) == 1:
            effective_sample = samples[0]

    computer = QCComputer(sample_id=effective_sample)

    try:
        for record in vcf:
            computer.process_record(record)
    finally:
        vcf.close()

    return computer.finalize()
