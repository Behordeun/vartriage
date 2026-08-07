"""Mitochondrial variant output serialization.

Converts MitoClassifiedVariant objects into dict/row formats compatible
with the JSON and CSV report writers. Provides a unified interface for
adding mitochondrial findings to both machine-readable and clinical reports.
"""

from __future__ import annotations

from typing import Any

from vartriage.mito.classifier import MitoClassifiedVariant


def mito_variant_to_dict(classified: MitoClassifiedVariant) -> dict[str, Any]:
    """Serialize a classified mitochondrial variant to a JSON-friendly dict.

    Output fields:
    - chromosome, position, ref_allele, alt_allele (standard VCF fields)
    - mt_classification (mtDNA-specific classification)
    - heteroplasmy_level (percentage, 0.0-100.0)
    - heteroplasmy_category (homoplasmic/high/moderate/low/sub_threshold)
    - heteroplasmy_depth (total read depth)
    - mitomap_disease (disease association or null)
    - mitomap_status (Cfrm/Reported/P.M. or null)
    - mt_gene (gene name or null)
    - mt_gene_type (protein_coding/tRNA/rRNA/control_region/intergenic)
    - helix_af (HelixMTdb allele frequency or null)
    - classification_reason (human-readable explanation)
    """
    variant = classified.variant
    record: dict[str, Any] = {
        "chromosome": variant.chrom,
        "position": variant.pos,
        "ref_allele": variant.ref,
        "alt_allele": variant.alt,
        "mt_classification": classified.classification.value,
        "heteroplasmy_level": None,
        "heteroplasmy_category": None,
        "heteroplasmy_depth": None,
        "mitomap_disease": None,
        "mitomap_status": None,
        "mt_gene": classified.gene_context.gene_name,
        "mt_gene_type": classified.gene_context.gene_type,
        "helix_af": classified.helix_af,
        "classification_reason": classified.classification_reason,
    }

    if classified.heteroplasmy is not None:
        record["heteroplasmy_level"] = round(classified.heteroplasmy.percentage, 2)
        record["heteroplasmy_category"] = classified.heteroplasmy.category
        record["heteroplasmy_depth"] = classified.heteroplasmy.depth

    if classified.mitomap_entry is not None:
        record["mitomap_disease"] = classified.mitomap_entry.disease
        record["mitomap_status"] = classified.mitomap_entry.status

    return record


MITO_CSV_FIELDS: list[str] = [
    "chromosome",
    "position",
    "ref_allele",
    "alt_allele",
    "mt_gene",
    "mt_gene_type",
    "mt_classification",
    "heteroplasmy_level",
    "heteroplasmy_category",
    "heteroplasmy_depth",
    "mitomap_disease",
    "mitomap_status",
    "helix_af",
    "classification_reason",
]


def mito_variant_to_row(classified: MitoClassifiedVariant) -> list[str]:
    """Serialize a classified mitochondrial variant to a CSV row.

    Returns field values in the order defined by MITO_CSV_FIELDS.
    """
    d = mito_variant_to_dict(classified)
    row: list[str] = []
    for field in MITO_CSV_FIELDS:
        value = d.get(field)
        row.append("" if value is None else str(value))
    return row


def build_mito_json_section(
    results: list[MitoClassifiedVariant],
) -> list[dict[str, Any]]:
    """Build the mitochondrial findings section for JSON output.

    Parameters
    ----------
    results
        Classified mitochondrial variants from the pipeline.

    Returns
    -------
    list of dict
        Serialized mitochondrial findings ready for JSON output.
    """
    return [mito_variant_to_dict(r) for r in results]
