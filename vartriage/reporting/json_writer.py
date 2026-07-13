"""JSON report writer.

Streams classified variants to RFC 8259 JSON with deterministic field
ordering. Only one variant is in memory at a time (beyond the I/O buffer).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterator, Sequence, Union

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
    """Flatten a ClassifiedVariant into an output-ordered dict."""
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
    variants: Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]],
    output_path: Path,
) -> Path:
    """Write classified variants to a JSON file, streaming one at a time.

    Parameters
    ----------
    variants : Union[Iterator[ClassifiedVariant], Sequence[ClassifiedVariant]]
        Variants in priority order. Iterators are consumed lazily.
    output_path : Path
        Destination file path.

    Returns
    -------
    Path
        The written file path.

    Raises
    ------
    IOError
        If the write fails (filesystem or encoding error).
    """
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("[\n")
            first = True
            for variant in variants:
                if not first:
                    f.write(",\n")
                json.dump(
                    _variant_to_dict(variant),
                    f,
                    ensure_ascii=False,
                    indent=2,
                    allow_nan=False,
                )
                first = False
            f.write("\n]\n")
    except (OSError, ValueError, TypeError) as exc:
        try:
            if output_path.exists():
                output_path.unlink()
        except OSError:
            pass
        raise IOError(f"Failed to write JSON report: {exc}") from exc

    return output_path
