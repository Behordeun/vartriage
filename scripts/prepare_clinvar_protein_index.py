#!/usr/bin/env python3
"""Generate a ClinVar protein index TSV for PS1/PM5 evidence criteria.

Reads a ClinVar VCF, filters for Pathogenic/Likely_Pathogenic missense
variants, resolves amino acid changes using a reference FASTA and gene
annotation (GTF), and writes a TSV keyed by (gene, amino_acid_position).

Output format:
    gene    position    ref_aa    alt_aa    chrom    pos    ref    alt    significance

Usage:
    python scripts/prepare_clinvar_protein_index.py \\
        --clinvar-vcf clinvar_20240101.vcf.gz \\
        --reference-fasta GRCh38.fa \\
        --gene-annotation gencode.v44.gtf \\
        --output clinvar_protein_index.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)


def _gtf_attr(attrs: str, key: str) -> str | None:
    """Extract a value from a GTF attribute string: key \"value\"; ..."""
    marker = f'{key} "'
    i = attrs.find(marker)
    if i == -1:
        return None
    i += len(marker)
    j = attrs.find('"', i)
    return attrs[i:j] if j != -1 else None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate ClinVar protein index for PS1/PM5 criteria"
    )
    parser.add_argument(
        "--clinvar-vcf",
        type=Path,
        required=True,
        help="Path to ClinVar VCF (bgzipped with .tbi index)",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        required=True,
        help="Indexed reference FASTA (.fa with .fai)",
    )
    parser.add_argument(
        "--gene-annotation",
        type=Path,
        required=True,
        help="GENCODE GTF file for transcript CDS mapping",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV path",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    return parser.parse_args()


def _is_pathogenic(clnsig: str) -> bool:
    """Check if a ClinVar CLNSIG value indicates pathogenicity."""
    sig_lower = clnsig.lower()
    return "pathogenic" in sig_lower and "conflicting" not in sig_lower


def _validate_inputs(args: argparse.Namespace) -> None:
    for path, label in [
        (args.clinvar_vcf, "ClinVar VCF"),
        (args.reference_fasta, "Reference FASTA"),
        (args.gene_annotation, "Gene annotation"),
    ]:
        if not path.exists():
            logger.error("%s not found: %s", label, path)
            sys.exit(1)


def _resolve_snv_entry(
    record: pysam.VariantRecord,
    alt: str,
    clnsig: str,
    resolver: object,
    skipped: list[int],
) -> tuple | None:
    if len(record.ref) != 1 or len(alt) != 1:
        return None
    # ClinVar VCF uses '1'/'MT'; the FASTA/GTF use 'chr1'/'chrM'. Normalize.
    chrom = record.chrom
    if not chrom.startswith("chr"):
        chrom = "chrM" if chrom in ("MT", "M") else f"chr{chrom}"
    try:
        context = resolver.resolve(chrom, record.pos, record.ref, alt)  # type: ignore[attr-defined]
    except (KeyError, ValueError, AttributeError, IndexError):
        skipped[0] += 1
        return None
    if context is None:
        skipped[0] += 1
        return None
    # Index only true missense substitutions: exclude synonymous, stop-gain
    # (altered "*") and stop-loss (reference "*").
    if (
        context.reference_aa == context.altered_aa
        or context.altered_aa == "*"
        or context.reference_aa == "*"
    ):
        return None
    return (
        context.gene_name or "UNKNOWN",
        context.codon_index + 1,
        context.reference_aa,
        context.altered_aa,
        chrom,
        record.pos,
        record.ref,
        alt,
        clnsig,
    )


def _process_vcf(
    vcf_path: Path,
    resolver: object,
) -> tuple[list[tuple], int, int]:
    vcf = pysam.VariantFile(str(vcf_path))
    entries: list[tuple] = []
    skipped = [0]
    processed = 0
    for record in vcf:
        processed += 1
        if processed % 50000 == 0:
            logger.info(
                "Processed %d records, %d entries collected...", processed, len(entries)
            )
        clnsig = record.info.get("CLNSIG")
        if clnsig is None:
            continue
        if isinstance(clnsig, tuple):
            clnsig = ",".join(str(s) for s in clnsig)
        if not _is_pathogenic(str(clnsig)):
            continue
        for alt in record.alts or []:
            entry = _resolve_snv_entry(record, alt, str(clnsig), resolver, skipped)
            if entry is not None:
                entries.append(entry)
    vcf.close()
    return entries, processed, skipped[0]


def _write_output(output: Path, entries: list[tuple]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as f:
        f.write("gene\tposition\tref_aa\talt_aa\tchrom\tpos\tref\talt\tsignificance\n")
        for entry in sorted(entries):
            f.write("\t".join(str(v) for v in entry) + "\n")


def main() -> None:
    args = _parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    _validate_inputs(args)

    from vartriage.annotation.codon_resolver import CodonResolver
    from vartriage.annotation.transcript_index import TranscriptCDSIndex

    logger.info("Building transcript CDS index from %s...", args.gene_annotation)
    cds_index = TranscriptCDSIndex()
    gtf_path = str(args.gene_annotation)
    gtf_open = gzip.open if gtf_path.endswith(".gz") else open
    with gtf_open(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GTF is 1-based; convert to 0-based
            end = int(parts[4])  # GTF end inclusive -> exclusive
            strand = parts[6]
            try:
                frame = int(parts[7]) if parts[7] != "." else 0
            except ValueError:
                frame = 0
            attrs = parts[8]
            transcript_id = _gtf_attr(attrs, "transcript_id")
            if not transcript_id:
                continue
            cds_index.add_cds_exon(
                transcript_id=transcript_id,
                gene_name=_gtf_attr(attrs, "gene_name") or "unknown",
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                frame=frame,
            )
    cds_index.finalize()
    count = (
        cds_index.transcript_count()
        if callable(cds_index.transcript_count)
        else cds_index.transcript_count
    )
    logger.info("Loaded %d transcripts", count)

    resolver = CodonResolver(
        fasta_path=args.reference_fasta, transcript_index=cds_index
    )

    logger.info("Processing ClinVar VCF: %s", args.clinvar_vcf)
    try:
        entries, processed, skipped = _process_vcf(args.clinvar_vcf, resolver)
    finally:
        resolver.close()

    logger.info(
        "Done. %d pathogenic missense entries from %d records (%d skipped)",
        len(entries),
        processed,
        skipped,
    )
    _write_output(args.output, entries)
    logger.info("Written to %s", args.output)


if __name__ == "__main__":
    main()
