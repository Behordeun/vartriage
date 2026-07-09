"""JSON report writer for serializing classified variants.

Produces RFC 8259 compliant JSON with UTF-8 encoding. Field ordering is
deterministic and matches the clinical report specification. Round-trip
fidelity is guaranteed: serialize then deserialize yields identical values
and ordering.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

from vartriage.models.variant import ClassifiedVariant


# Output field order as specified in requirements
_OUTPUT_FIELDS: tuple[str, ...] = (
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
)


def _variant_to_dict(variant: ClassifiedVariant) -> dict[str, Any]:
    """Convert a ClassifiedVariant to an ordered dictionary for JSON output.

    Parameters
    ----------
    variant : ClassifiedVariant
        The classified variant to serialize.

    Returns
    -------
    dict[str, Any]
        Dictionary with output fields in specified order, absent values as None.
    """
    scored = variant.scored
    annotated = scored.annotated
    raw = annotated.variant

    consequence: str | None = None
    if annotated.consequence is not None:
        consequence = annotated.consequence.value

    clinvar: str | None = None
    if annotated.clinvar_assertion is not None:
        clinvar = annotated.clinvar_assertion.value

    classification: str | None = None
    if variant.classification is not None:
        classification = variant.classification.value

    evidence: list[str] | None = None
    if variant.evidence_tags:
        evidence = sorted(tag.value for tag in variant.evidence_tags)
    else:
        evidence = []

    record: dict[str, Any] = {}
    record["chromosome"] = raw.chrom
    record["position"] = raw.pos
    record["ref_allele"] = raw.ref
    record["alt_allele"] = raw.alt
    record["functional_consequence"] = consequence
    record["allele_frequency"] = annotated.allele_frequency
    record["composite_rank"] = scored.composite_rank
    record["clinvar_assertion"] = clinvar
    record["acmg_classification"] = classification
    record["evidence_tags"] = evidence

    return record


def write_json(
    variants: Sequence[ClassifiedVariant],
    output_path: Path,
) -> Path:
    """Serialize a list of classified variants to a JSON file.

    Produces RFC 8259 compliant JSON with UTF-8 encoding. All output fields
    are included in the specified order. Absent values are represented as
    JSON ``null``. Round-trip fidelity is guaranteed: deserializing the
    output produces identical values, types, and ordering.

    Parameters
    ----------
    variants : Sequence[ClassifiedVariant]
        The classified variants to serialize, in prioritized rank order.
    output_path : Path
        Destination file path for the JSON output.

    Returns
    -------
    Path
        The path to the written JSON file.

    Raises
    ------
    IOError
        If the file cannot be written due to a filesystem or encoding error.
    """
    records = [_variant_to_dict(v) for v in variants]

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(records, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
    except (OSError, ValueError, TypeError) as exc:
        try:
            if output_path.exists():
                output_path.unlink()
        except OSError:
            pass
        raise IOError(f"Failed to write JSON report: {exc}") from exc

    return output_path
