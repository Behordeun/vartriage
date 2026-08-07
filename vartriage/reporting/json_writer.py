"""JSON report writer.

Streams classified variants to RFC 8259 JSON with deterministic field
ordering. Only one variant is in memory at a time (beyond the I/O buffer).
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

from vartriage._internal.path_safety import resolve_path
from vartriage.models.variant import ClassifiedVariant

# Output field order as specified in requirements
_OUTPUT_FIELDS: tuple[str, ...] = (
    "chromosome",
    "position",
    "ref_allele",
    "alt_allele",
    "gene_name",
    "functional_consequence",
    "allele_frequency",
    "revel_score",
    "composite_rank",
    "prioritization_score",
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
    record["gene_name"] = annotated.gene_name
    record["functional_consequence"] = consequence
    record["allele_frequency"] = annotated.allele_frequency
    record["revel_score"] = scored.revel_score
    record["composite_rank"] = scored.composite_rank
    record["prioritization_score"] = scored.prioritization_score
    record["clinvar_assertion"] = clinvar
    record["acmg_classification"] = classification
    record["evidence_tags"] = evidence

    # Gene-disease linkage context (when knowledge base is active)
    if annotated.gene_context is not None:
        ctx = annotated.gene_context
        record["disease_associations"] = [
            {
                "disease_name": a.disease_name,
                "mim_number": a.mim_number,
                "inheritance_mode": a.inheritance_mode,
            }
            for a in ctx.disease_associations
        ]
        record["clingen_validity"] = ctx.clingen_validity
        if ctx.constraint is not None:
            record["gene_constraint"] = {
                "pli": ctx.constraint.pli,
                "loeuf": ctx.constraint.loeuf,
                "mis_z": ctx.constraint.mis_z,
            }
        else:
            record["gene_constraint"] = None
        record["is_actionable"] = ctx.is_actionable
        record["phenotype_match_score"] = ctx.phenotype_match_score

    return record


def write_json(
    variants: Iterator[ClassifiedVariant] | Sequence[ClassifiedVariant],
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
        output_path = resolve_path(output_path)
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
        raise OSError(f"Failed to write JSON report: {exc}") from exc

    return output_path


def write_json_with_mito(
    variants: Iterator[ClassifiedVariant] | Sequence[ClassifiedVariant],
    output_path: Path,
    mito_results: list[Any] | None = None,
) -> Path:
    """Write classified variants to JSON with optional mitochondrial section.

    When mito_results is provided, outputs a structured object:
    {"variants": [...], "mitochondrial_findings": [...]}

    When mito_results is None or empty, delegates to write_json for
    backward-compatible flat array output.

    Parameters
    ----------
    variants
        Nuclear classified variants.
    output_path
        Destination file path.
    mito_results
        Optional list of MitoClassifiedVariant objects.

    Returns
    -------
    Path
        The written file path.
    """
    if not mito_results:
        return write_json(variants, output_path)

    from vartriage.mito.report import build_mito_json_section

    try:
        output_path = resolve_path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        mito_section = build_mito_json_section(mito_results)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write('{\n  "variants": [\n')
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
            f.write("\n  ],\n")
            f.write('  "mitochondrial_findings": ')
            json.dump(mito_section, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write(",\n")
            f.write('  "metadata": {\n')
            f.write(
                '    "mitochondrial_note": '
                '"Mitochondrial variants classified using mtDNA-specific criteria"\n'
            )
            f.write("  }\n}\n")
    except (OSError, ValueError, TypeError) as exc:
        try:
            if output_path.exists():
                output_path.unlink()
        except OSError:
            pass
        raise OSError(f"Failed to write JSON report: {exc}") from exc

    return output_path
