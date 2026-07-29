"""Structural variant gene annotation and impact classification.

Determines which genes an SV overlaps, classifies the impact type
(whole-gene deletion, partial overlap, gene disruption, etc.), and
integrates dosage sensitivity data from ClinGen to identify
haploinsufficient and triplosensitive genes.

Gene coordinates are loaded from a GTF file. Dosage sensitivity is
loaded from a ClinGen TSV. Population frequency is matched via
reciprocal overlap against a gnomAD-SV reference.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

from vartriage.structural.models import (
    AnnotatedSV,
    GeneOverlap,
    SVConsequence,
    SVType,
    StructuralVariant,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class GeneRecord:
    """Minimal gene annotation from GTF."""

    symbol: str
    chrom: str
    start: int
    end: int
    strand: str
    exon_count: int


@dataclass(frozen=True, slots=True)
class DosageEntry:
    """ClinGen dosage sensitivity for a single gene."""

    gene_symbol: str
    hi_score: Optional[float]
    ts_score: Optional[float]


@dataclass(frozen=True, slots=True)
class SVFrequencyRecord:
    """Reference SV from gnomAD-SV for frequency matching."""

    chrom: str
    start: int
    end: int
    sv_type: str
    allele_frequency: float


class SVAnnotator:
    """Annotate structural variants with gene overlap and dosage data.

    Loads reference data at construction time, then annotates a stream
    of StructuralVariant records. Each SV gets:
    - Gene overlap assessment (which genes, overlap fraction, exon count)
    - Impact classification (whole-gene deletion, partial, disruption)
    - Dosage sensitivity flags (HI/TS scores from ClinGen)
    - Population frequency (reciprocal overlap against gnomAD-SV)

    Parameters
    ----------
    gene_annotation_path : Optional[Path]
        GTF file for gene coordinates. When None, all SVs get
        INTERGENIC consequence with no gene overlaps.
    dosage_sensitivity_path : Optional[Path]
        ClinGen dosage sensitivity TSV. When None, dosage scores
        are not attached.
    gnomad_sv_path : Optional[Path]
        gnomAD-SV BED/TSV for frequency matching. When None, all
        SVs are frequency-unknown.
    reciprocal_overlap : float
        Minimum reciprocal overlap fraction for frequency matching.
    whole_gene_threshold : float
        Fraction of gene covered to classify as whole-gene event.
    """

    def __init__(
        self,
        gene_annotation_path: Optional[Path] = None,
        dosage_sensitivity_path: Optional[Path] = None,
        gnomad_sv_path: Optional[Path] = None,
        reciprocal_overlap: float = 0.5,
        whole_gene_threshold: float = 0.8,
    ) -> None:
        self._reciprocal_overlap = reciprocal_overlap
        self._whole_gene_threshold = whole_gene_threshold

        self._genes: dict[str, list[GeneRecord]] = {}
        self._dosage: dict[str, DosageEntry] = {}
        self._sv_database: dict[str, list[SVFrequencyRecord]] = {}

        if gene_annotation_path is not None:
            self._load_genes(gene_annotation_path)
        if dosage_sensitivity_path is not None:
            self._load_dosage(dosage_sensitivity_path)
        if gnomad_sv_path is not None:
            self._load_sv_frequency(gnomad_sv_path)

    def annotate(
        self, variants: Iterator[StructuralVariant]
    ) -> Iterator[AnnotatedSV]:
        """Annotate a stream of structural variants.

        Parameters
        ----------
        variants : Iterator[StructuralVariant]
            Parsed SV records from the VCF parser.

        Yields
        ------
        AnnotatedSV
            Each SV enriched with gene overlap, consequence, and
            population frequency data.
        """
        for sv in variants:
            yield self._annotate_single(sv)

    def _annotate_single(self, sv: StructuralVariant) -> AnnotatedSV:
        """Annotate one SV with gene overlaps and frequency."""
        gene_overlaps = self._find_gene_hits(sv)
        consequence = self._sv_consequence(sv, gene_overlaps)
        pop_freq, freq_unknown = self._lookup_frequency(sv)

        genes_affected = len(gene_overlaps)
        hi_count = sum(1 for g in gene_overlaps if g.is_haploinsufficient)

        return AnnotatedSV(
            sv=sv,
            consequence=consequence,
            gene_overlaps=tuple(gene_overlaps),
            population_frequency=pop_freq,
            frequency_unknown=freq_unknown,
            genes_affected=genes_affected,
            hi_genes_affected=hi_count,
        )

    def _gene_overlap(self, sv: StructuralVariant, gene: "GeneRecord") -> GeneOverlap | None:
        """Return a GeneOverlap for one gene, or None if no overlap."""
        overlap_start = max(sv.start, gene.start)
        overlap_end = min(sv.end, gene.end)
        if overlap_start > overlap_end:
            return None

        gene_length = gene.end - gene.start + 1
        overlap_bp = overlap_end - overlap_start + 1
        overlap_fraction = overlap_bp / gene_length if gene_length > 0 else 0.0
        is_whole = overlap_fraction >= self._whole_gene_threshold
        exons_affected = gene.exon_count if is_whole else max(1, int(gene.exon_count * overlap_fraction))

        dosage = self._dosage.get(gene.symbol)
        hi_score = dosage.hi_score if dosage else None
        ts_score = dosage.ts_score if dosage else None

        return GeneOverlap(
            gene_symbol=gene.symbol,
            gene_chrom=gene.chrom,
            gene_start=gene.start,
            gene_end=gene.end,
            overlap_fraction=overlap_fraction,
            is_whole_gene=is_whole,
            exons_affected=exons_affected,
            total_exons=gene.exon_count,
            is_haploinsufficient=hi_score is not None and hi_score >= 3.0,
            is_triplosensitive=ts_score is not None and ts_score >= 3.0,
            hi_score=hi_score,
            ts_score=ts_score,
        )

    def _find_gene_hits(self, sv: StructuralVariant) -> list[GeneOverlap]:
        """Find all protein-coding genes overlapping the SV span."""
        chrom_genes = self._genes.get(sv.chrom, [])
        if not chrom_genes:
            alt_chrom = sv.chrom.replace("chr", "") if sv.chrom.startswith("chr") else f"chr{sv.chrom}"
            chrom_genes = self._genes.get(alt_chrom, [])

        hits = [self._gene_overlap(sv, g) for g in chrom_genes]
        result = [h for h in hits if h is not None]
        result.sort(key=lambda g: g.overlap_fraction, reverse=True)
        return result

    @staticmethod
    def _del_dup_consequence(sv_type: SVType, is_whole: bool) -> SVConsequence:
        """Consequence for DEL or DUP based on whole-gene status."""
        if sv_type == SVType.DEL:
            return SVConsequence.WHOLE_GENE_DELETION if is_whole else SVConsequence.PARTIAL_GENE_DELETION
        return SVConsequence.WHOLE_GENE_DUPLICATION if is_whole else SVConsequence.PARTIAL_GENE_DUPLICATION

    @staticmethod
    def _cnv_consequence(sv: StructuralVariant, is_whole: bool) -> SVConsequence:
        """Consequence for CNV based on copy number."""
        if sv.copy_number is not None:
            if sv.copy_number < 2:
                return SVConsequence.WHOLE_GENE_DELETION if is_whole else SVConsequence.PARTIAL_GENE_DELETION
            if sv.copy_number > 2:
                return SVConsequence.WHOLE_GENE_DUPLICATION if is_whole else SVConsequence.PARTIAL_GENE_DUPLICATION
        return SVConsequence.GENE_DISRUPTION

    def _sv_consequence(
        self,
        sv: StructuralVariant,
        gene_hits: list[GeneOverlap],
    ) -> SVConsequence:
        """Determine the most severe gene-level consequence."""
        if not gene_hits:
            return SVConsequence.INTERGENIC

        is_whole = gene_hits[0].is_whole_gene

        if sv.sv_type in (SVType.DEL, SVType.DUP):
            return self._del_dup_consequence(sv.sv_type, is_whole)
        if sv.sv_type in (SVType.INV, SVType.BND, SVType.INS):
            return SVConsequence.GENE_DISRUPTION
        if sv.sv_type == SVType.CNV:
            return self._cnv_consequence(sv, is_whole)
        return SVConsequence.INTERGENIC

    def _reciprocal_freq(
        self, sv_start: int, sv_end: int, sv_type_str: str, ref_sv: object
    ) -> tuple[float, Optional[float]] | None:
        """Return (reciprocal_overlap, allele_frequency) for a matching ref SV, or None."""
        if ref_sv.sv_type != sv_type_str:  # type: ignore[attr-defined]
            return None
        overlap_start = max(sv_start, ref_sv.start)  # type: ignore[attr-defined]
        overlap_end = min(sv_end, ref_sv.end)  # type: ignore[attr-defined]
        if overlap_start > overlap_end:
            return None
        sv_length = sv_end - sv_start + 1
        ref_length = ref_sv.end - ref_sv.start + 1  # type: ignore[attr-defined]
        overlap_bp = overlap_end - overlap_start + 1
        frac_query = overlap_bp / sv_length if sv_length > 0 else 0.0
        frac_ref = overlap_bp / ref_length if ref_length > 0 else 0.0
        return min(frac_query, frac_ref), ref_sv.allele_frequency  # type: ignore[attr-defined]

    def _lookup_frequency(
        self, sv: StructuralVariant
    ) -> tuple[Optional[float], bool]:
        """Match SV against gnomAD-SV using reciprocal overlap.

        Returns (frequency, frequency_unknown). If no reference database
        is loaded or no match found, returns (None, True).
        """
        if not self._sv_database:
            return None, True

        chrom_svs = self._sv_database.get(sv.chrom, [])
        if not chrom_svs:
            alt_chrom = sv.chrom.replace("chr", "") if sv.chrom.startswith("chr") else f"chr{sv.chrom}"
            chrom_svs = self._sv_database.get(alt_chrom, [])

        best_freq: Optional[float] = None
        best_overlap: float = 0.0

        for ref_sv in chrom_svs:
            match = self._reciprocal_freq(sv.start, sv.end, sv.sv_type.value, ref_sv)
            if match is not None and match[0] >= self._reciprocal_overlap and match[0] > best_overlap:
                best_overlap, best_freq = match

        return (best_freq, False) if best_freq is not None else (None, True)

    @staticmethod
    def _parse_gtf_gene_line(
        fields: list[str],
        gene_records: dict[str, tuple[str, int, int, str]],
        gene_exon_counts: dict[str, int],
    ) -> None:
        """Update gene_records and gene_exon_counts from one GTF fields list."""
        feature_type = fields[2]
        if feature_type not in ("gene", "exon"):
            return
        attributes = fields[8]
        gene_name = _extract_attribute(attributes, "gene_name")
        if gene_name is None:
            return
        gene_type = _extract_attribute(attributes, "gene_type") or _extract_attribute(attributes, "gene_biotype")
        if gene_type is not None and gene_type != "protein_coding":
            return
        if feature_type == "gene":
            gene_records[gene_name] = (fields[0], int(fields[3]), int(fields[4]), fields[6])
            gene_exon_counts.setdefault(gene_name, 0)
        else:
            gene_exon_counts[gene_name] = gene_exon_counts.get(gene_name, 0) + 1

    def _load_genes(self, gtf_path: Path) -> None:
        """Load protein-coding gene coordinates from GTF.

        Extracts gene-level features and counts exons per gene
        for overlap assessment.
        """
        gtf_path = gtf_path.resolve()
        if not gtf_path.exists():
            raise FileNotFoundError(f"Gene annotation not found: {gtf_path}")

        gene_exon_counts: dict[str, int] = {}
        gene_records: dict[str, tuple[str, int, int, str]] = {}

        with open(gtf_path, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 9:
                    self._parse_gtf_gene_line(fields, gene_records, gene_exon_counts)

        for symbol, (chrom, start, end, strand) in gene_records.items():
            record = GeneRecord(
                symbol=symbol, chrom=chrom, start=start, end=end, strand=strand,
                exon_count=max(1, gene_exon_counts.get(symbol, 1)),
            )
            self._genes.setdefault(chrom, []).append(record)

        for chrom in self._genes:
            self._genes[chrom].sort(key=lambda g: g.start)

        total_genes = sum(len(gs) for gs in self._genes.values())
        logger.info("Loaded %d protein-coding genes from %s", total_genes, gtf_path)

    def _load_dosage(self, path: Path) -> None:
        """Load ClinGen dosage sensitivity scores.

        Expects a TSV with columns: gene_symbol, hi_score, ts_score.
        Scores are numeric (0-3 scale where 3 = sufficient evidence).
        """
        path = path.resolve()
        if not path.exists():
            raise FileNotFoundError(f"Dosage sensitivity file not found: {path}")

        with open(path, "r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                symbol = row.get("gene_symbol") or row.get("Gene Symbol") or ""
                symbol = symbol.strip()
                if not symbol:
                    continue

                hi_raw = row.get("hi_score") or row.get("Haploinsufficiency Score")
                ts_raw = row.get("ts_score") or row.get("Triplosensitivity Score")

                hi_score = _parse_score(hi_raw)
                ts_score = _parse_score(ts_raw)

                self._dosage[symbol] = DosageEntry(
                    gene_symbol=symbol,
                    hi_score=hi_score,
                    ts_score=ts_score,
                )

        logger.info("Loaded dosage sensitivity for %d genes from %s", len(self._dosage), path)

    def _load_sv_frequency(self, path: Path) -> None:
        """Load gnomAD-SV reference for population frequency matching.

        Expects a TSV/BED with columns: chrom, start, end, sv_type, af.
        """
        path = path.resolve()
        if not path.exists():
            raise FileNotFoundError(f"gnomAD-SV file not found: {path}")

        with open(path, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue

                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue

                sv_type = parts[3].upper()
                try:
                    af = float(parts[4])
                except ValueError:
                    continue

                record = SVFrequencyRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    sv_type=sv_type,
                    allele_frequency=af,
                )
                if chrom not in self._sv_database:
                    self._sv_database[chrom] = []
                self._sv_database[chrom].append(record)

        # Sort for efficient overlap scanning
        for chrom in self._sv_database:
            self._sv_database[chrom].sort(key=lambda r: r.start)

        total = sum(len(rs) for rs in self._sv_database.values())
        logger.info("Loaded %d reference SVs from %s", total, path)


def _extract_attribute(attributes: str, key: str) -> Optional[str]:
    """Extract a value from GTF attribute string (key "value"; format)."""
    for attr in attributes.split(";"):
        attr = attr.strip()
        if attr.startswith(key):
            parts = attr.split('"')
            if len(parts) >= 2:
                return parts[1]
            # Handle key value (no quotes)
            parts = attr.split()
            if len(parts) >= 2:
                return parts[1].strip('"')
    return None


def _parse_score(raw: Optional[str]) -> Optional[float]:
    """Parse a numeric score, returning None for missing/invalid values."""
    if raw is None:
        return None
    raw = raw.strip()
    if not raw or raw in ("", "N/A", "NA", "-", "Not yet evaluated"):
        return None
    try:
        return float(raw)
    except ValueError:
        return None
