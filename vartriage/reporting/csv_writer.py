"""CSV report writer for classified variant output.

Serializes a sequence of ClassifiedVariant records into an RFC 4180 compliant
CSV file with UTF-8 encoding. Each row represents one variant, with absent
values represented as empty fields.

Supports both Iterator and Sequence inputs — variants are written row by row
as they are consumed from the iterator.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator, Sequence, Union

from vartriage.models.variant import ClassifiedVariant


CSV_FIELDS: list[str] = [
    "chromosome",
    "position",
    "ref_allele",
    "alt_allele",
    "functional_consequence",
    "allele_frequency",
    "composite_rank",
    "clinvar_assertion",
    "acmg_classification",
    "evidence_tags",
]
"""Output field names in the order specified by the report schema."""


def _format_field(value: object) -> str:
    """Convert a field value to its CSV string representation.

    Parameters
    ----------
    value : object
        The value to format. None becomes an empty string; all other values
        are converted via ``str()``.

    Returns
    -------
    str
        The formatted string ready for CSV output.
    """
    if value is None:
        return ""
    return str(value)


def _variant_to_row(variant: ClassifiedVariant) -> list[str]:
    """Extract output fields from a ClassifiedVariant in the canonical order.

    Parameters
    ----------
    variant : ClassifiedVariant
        A fully classified variant record.

    Returns
    -------
    list[str]
        Field values in the order defined by ``CSV_FIELDS``.
    """
    scored = variant.scored
    annotated = scored.annotated
    base = annotated.variant

    consequence_value = (
        annotated.consequence.value
        if annotated.consequence is not None
        else None
    )
    clinvar_value = (
        annotated.clinvar_assertion.value
        if annotated.clinvar_assertion is not None
        else None
    )
    classification_value = (
        variant.classification.value
        if variant.classification is not None
        else None
    )

    evidence_tags_value: str | None
    if variant.evidence_tags:
        evidence_tags_value = ";".join(
            sorted(tag.value for tag in variant.evidence_tags)
        )
    else:
        evidence_tags_value = None

    return [
        _format_field(base.chrom),
        _format_field(base.pos),
        _format_field(base.ref),
        _format_field(base.alt),
        _format_field(consequence_value),
        _format_field(annotated.allele_frequency),
        _format_field(scored.composite_rank),
        _format_field(clinvar_value),
        _format_field(classification_value),
        _format_field(evidence_tags_value),
    ]


def write_csv(
    variants: Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]],
    output_path: Path,
) -> Path:
    """Serialize classified variants to an RFC 4180 compliant CSV file.

    Parameters
    ----------
    variants : Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]]
        The prioritized variant list or iterator to serialize. Variants
        are written row by row as they are consumed.
    output_path : Path
        Destination file path for the CSV output.

    Returns
    -------
    Path
        The path to the written CSV file (same as ``output_path``).

    Raises
    ------
    IOError
        If the file cannot be written due to filesystem or encoding errors.
    """
    with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(CSV_FIELDS)
        for variant in variants:
            writer.writerow(_variant_to_row(variant))

    return output_path
