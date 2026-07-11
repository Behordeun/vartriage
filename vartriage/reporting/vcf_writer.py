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
        "Description": (
            "Functional consequence assigned by vartriage"
        ),
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
        "Description": (
            "ACMG/AMP classification assigned by vartriage"
        ),
    },
    {
        "ID": "VARTRIAGE_TAGS",
        "Number": "1",
        "Type": "String",
        "Description": (
            "Comma-separated ACMG evidence tags assigned by vartriage"
        ),
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
            '##INFO=<ID={ID},Number={Number},Type={Type},'
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
        tags_str = ",".join(
            sorted(tag.value for tag in classified.evidence_tags)
        )
        record.info["VARTRIAGE_TAGS"] = tags_str


def write_vcf(
    variants: Sequence[ClassifiedVariant],
    source_vcf_path: Path,
    output_path: Path,
) -> Path:
    """Write annotated VCF with VARTRIAGE_* INFO fields.

    Re-reads the source VCF, matches records to classified variants
    by (chrom, pos, ref, alt), injects INFO fields for matches, and
    writes all records to a bgzipped output with a tabix index.

    Uses an atomic write pattern: writes to a temp file first, then
    renames to the target path before generating the tabix index.

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
    tmp_path = output_path.with_suffix(".vcf.gz.tmp")

    try:
        with pysam.VariantFile(str(source_vcf_path), "r") as src:
            new_header = _add_info_headers(src.header.copy())

            with pysam.VariantFile(
                str(tmp_path), "wz", header=new_header
            ) as out:
                for record in src:
                    # pysam stubs type translate() as None incorrectly
                    new_rec: pysam.VariantRecord = record.translate(  # type: ignore[func-returns-value,assignment]  # noqa: E501
                        out.header
                    )
                    alts = record.alts
                    if alts and alts[0] is not None and record.ref is not None:
                        key: LookupKey = (
                            record.chrom,
                            record.pos,
                            str(record.ref),
                            str(alts[0]),
                        )
                        if key in lookup:
                            _inject_info_fields(new_rec, lookup[key])

                    out.write(new_rec)

        os.replace(str(tmp_path), str(output_path))
        pysam.tabix_index(str(output_path), preset="vcf", force=True)

    except Exception as exc:
        if tmp_path.exists():
            tmp_path.unlink()
        raise IOError(
            f"Failed to write VCF output: {exc}"
        ) from exc

    return output_path
