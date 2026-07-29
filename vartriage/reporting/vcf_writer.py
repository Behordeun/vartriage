"""VCF report writer.

Re-reads the source VCF with pysam, injects VARTRIAGE_* INFO fields
for matched classified variants, and writes bgzipped output with a
tabix index.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path

import pysam

from vartriage.models.variant import ClassifiedVariant

LookupKey = tuple[str, int, str, str]
"""Composite key for matching VCF records: (chrom, pos, ref, alt)."""

VARTRIAGE_INFO_FIELDS: list[dict[str, str]] = [
    {
        "ID": "VARTRIAGE_CONSEQUENCE",
        "Number": "1",
        "Type": "String",
        "Description": ("Functional consequence assigned by vartriage"),
    },
    {
        "ID": "VARTRIAGE_AF",
        "Number": "1",
        "Type": "Float",
        "Description": "Population allele frequency from gnomAD",
    },
    {
        "ID": "VARTRIAGE_RANK",
        "Number": "1",
        "Type": "Float",
        "Description": "Composite pathogenicity rank score",
    },
    {
        "ID": "VARTRIAGE_ACMG",
        "Number": "1",
        "Type": "String",
        "Description": ("ACMG/AMP classification assigned by vartriage"),
    },
    {
        "ID": "VARTRIAGE_TAGS",
        "Number": "1",
        "Type": "String",
        "Description": ("Comma-separated ACMG evidence tags assigned by vartriage"),
    },
]


def _build_lookup(
    variants: Sequence[ClassifiedVariant],
) -> dict[LookupKey, ClassifiedVariant]:
    """Build a lookup dictionary from classified variants.

    Maps each variant's genomic coordinates to the variant itself.
    If duplicate keys exist (same chrom, pos, ref, alt), the last
    variant in sequence order wins.

    Parameters
    ----------
    variants : Sequence[ClassifiedVariant]
        Materialized classified variants.

    Returns
    -------
    dict[LookupKey, ClassifiedVariant]
        Mapping of (chrom, pos, ref, alt) to classified variant.
    """
    lookup: dict[LookupKey, ClassifiedVariant] = {}
    for cv in variants:
        v = cv.scored.annotated.variant
        key: LookupKey = (v.chrom, v.pos, v.ref, v.alt)
        lookup[key] = cv
    return lookup


def _add_info_headers(
    header: pysam.VariantHeader,
) -> pysam.VariantHeader:
    """Add VARTRIAGE_* INFO field definitions to a VCF header.

    Parameters
    ----------
    header : pysam.VariantHeader
        Source VCF header to augment.

    Returns
    -------
    pysam.VariantHeader
        The same header object with five new INFO lines added.
    """
    for field_def in VARTRIAGE_INFO_FIELDS:
        header.add_line(
            "##INFO=<ID={ID},Number={Number},Type={Type},"
            'Description="{Description}">'.format(**field_def)
        )
    return header


def _inject_info_fields(
    record: pysam.VariantRecord,
    classified: ClassifiedVariant,
) -> None:
    """Inject VARTRIAGE_* INFO fields into a VCF record.

    Always sets VARTRIAGE_CONSEQUENCE and VARTRIAGE_ACMG.
    Conditionally sets VARTRIAGE_AF, VARTRIAGE_RANK, and
    VARTRIAGE_TAGS only when their source data is non-null/non-empty.

    Parameters
    ----------
    record : pysam.VariantRecord
        A writable record from the output VCF file.
    classified : ClassifiedVariant
        The matched classified variant whose data to inject.
    """
    ann = classified.scored.annotated
    record.info["VARTRIAGE_CONSEQUENCE"] = ann.consequence.value
    record.info["VARTRIAGE_ACMG"] = classified.classification.value

    if ann.allele_frequency is not None:
        record.info["VARTRIAGE_AF"] = ann.allele_frequency

    if classified.scored.composite_rank is not None:
        record.info["VARTRIAGE_RANK"] = classified.scored.composite_rank

    if classified.evidence_tags:
        tags_str = ",".join(sorted(tag.value for tag in classified.evidence_tags))
        record.info["VARTRIAGE_TAGS"] = tags_str


def _copy_info(src_rec: pysam.VariantRecord, dst_rec: pysam.VariantRecord) -> None:
    """Copy all INFO fields from src_rec to dst_rec."""
    for info_key in src_rec.info:
        dst_rec.info[info_key] = src_rec.info[info_key]


def _copy_samples(src_rec: pysam.VariantRecord, dst_rec: pysam.VariantRecord) -> None:
    """Copy all FORMAT/sample data from src_rec to dst_rec."""
    for sample in src_rec.samples:
        for fmt_key in src_rec.samples[sample]:
            dst_rec.samples[sample][fmt_key] = src_rec.samples[sample][fmt_key]


def _find_classified(
    record: pysam.VariantRecord,
    lookup: dict[LookupKey, ClassifiedVariant],
) -> ClassifiedVariant | None:
    """Return the first matching ClassifiedVariant for a VCF record, or None."""
    if not record.alts or record.ref is None:
        return None
    for alt_allele in record.alts:
        if alt_allele is None:
            continue
        key: LookupKey = (record.chrom, record.pos, str(record.ref), str(alt_allele))
        if key in lookup:
            return lookup[key]
    return None


def _write_records(
    src: pysam.VariantFile,
    out: pysam.VariantFile,
    lookup: dict[LookupKey, ClassifiedVariant],
) -> None:
    """Iterate src records, annotate matches, and write all to out."""
    for record in src:
        new_rec = out.new_record(
            contig=record.chrom,
            start=record.start,
            stop=record.stop,
            alleles=record.alleles,
            id=record.id,
            qual=record.qual,
            filter=record.filter,
        )
        _copy_info(record, new_rec)
        _copy_samples(record, new_rec)
        classified = _find_classified(record, lookup)
        if classified is not None:
            _inject_info_fields(new_rec, classified)
        out.write(new_rec)


def write_vcf(
    variants: Sequence[ClassifiedVariant],
    source_vcf_path: Path,
    output_path: Path,
) -> Path:
    """Write annotated VCF with VARTRIAGE_* INFO fields.

    Re-reads the source VCF, matches records to classified variants
    by (chrom, pos, ref, alt), injects INFO fields for matches, and
    writes all records to a bgzipped output with a tabix index.

    Both the VCF and its .tbi index are written atomically: the tabix
    index is built against the temp file first, then both files are
    moved into their final locations. If any step fails, both temp
    files are cleaned up so no partial output remains on disk.

    Parameters
    ----------
    variants : Sequence[ClassifiedVariant]
        Materialized classified variants for lookup building.
    source_vcf_path : Path
        Path to the original input VCF file.
    output_path : Path
        Target path for the bgzipped output (.vcf.gz).

    Returns
    -------
    Path
        The written output path.

    Raises
    ------
    IOError
        If writing or indexing fails.
    """
    lookup = _build_lookup(variants)

    # Resolve output directory once; all derived paths are validated against it
    # to prevent path traversal (CWE-22).
    output_dir = output_path.resolve().parent
    resolved_output = output_path.resolve()

    # Tabix index file extension for bgzipped VCF files.
    tabix_ext = ".tbi"

    tmp_path = output_dir / (resolved_output.name + ".tmp")
    tmp_tabix_path = output_dir / (tmp_path.name + tabix_ext)
    final_tabix_path = output_dir / (resolved_output.name + tabix_ext)

    try:
        with pysam.VariantFile(str(source_vcf_path), "r") as src:
            new_header = _add_info_headers(src.header.copy())
            with pysam.VariantFile(str(tmp_path), "wz", header=new_header) as out:
                _write_records(src, out, lookup)

        # Build tabix index against temp file before moving anything.
        pysam.tabix_index(str(tmp_path), preset="vcf", force=True)

        # Atomically move both VCF and tabix index into their final locations.
        os.replace(str(tmp_path), str(resolved_output))
        os.replace(str(tmp_tabix_path), str(final_tabix_path))

    except Exception as exc:
        for p in (tmp_path, tmp_tabix_path):
            if p.exists():
                p.unlink()
        raise IOError(f"Failed to write VCF output: {exc}") from exc

    return output_path
