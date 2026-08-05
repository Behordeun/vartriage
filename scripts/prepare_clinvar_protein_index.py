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
import logging
import sys
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)


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
    try:
        context = resolver.resolve(record.chrom, record.pos, record.ref, alt)  # type: ignore[attr-defined]
    except Exception:
        skipped[0] += 1
        return None
    if context is None:
        skipped[0] += 1
        return None
    if context.ref_aa == context.alt_aa or context.alt_aa == "*":
        return None
    return (
        context.gene_name or "UNKNOWN",
        context.aa_position,
        context.ref_aa,
        context.alt_aa,
        record.chrom,
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
    cds_index.load_from_gtf(args.gene_annotation)
    logger.info("Loaded %d transcripts", cds_index.transcript_count)

    fasta = pysam.FastaFile(str(args.reference_fasta))
    resolver = CodonResolver(fasta=fasta, transcript_index=cds_index)

    logger.info("Processing ClinVar VCF: %s", args.clinvar_vcf)
    entries, processed, skipped = _process_vcf(args.clinvar_vcf, resolver)
    fasta.close()

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
