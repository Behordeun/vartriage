"""Combined SNV + SV findings merger for unified reports.

When both a point-variant VCF and an SV VCF are analyzed in the same
session, this module merges findings into a single ranked output sorted
by clinical significance tier.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

from vartriage.models.variant import ACMGClassification, ClassifiedVariant
from vartriage.structural.models import ClassifiedSV, SVClassification


@dataclass(frozen=True)
class CombinedFinding:
    """A single variant finding (SNV or SV) with unified ranking.

    Parameters
    ----------
    variant_type : str
        "SNV" or "SV" indicating the source.
    tier : int
        Numeric tier for sorting (0=Pathogenic, 4=Benign).
    score : float
        Pathogenicity/composite score for secondary sort.
    data : dict[str, Any]
        Serialized finding data for report output.
    """

    variant_type: str
    tier: int
    score: float
    data: dict[str, Any] = field(default_factory=dict)


# Tier mappings for unified sorting
_SNV_TIER: dict[ACMGClassification, int] = {
    ACMGClassification.PATHOGENIC: 0,
    ACMGClassification.LIKELY_PATHOGENIC: 1,
    ACMGClassification.VUS: 2,
    ACMGClassification.LIKELY_BENIGN: 3,
    ACMGClassification.BENIGN: 4,
}

_SV_TIER: dict[SVClassification, int] = {
    SVClassification.PATHOGENIC: 0,
    SVClassification.LIKELY_PATHOGENIC: 1,
    SVClassification.VUS: 2,
    SVClassification.LIKELY_BENIGN: 3,
    SVClassification.BENIGN: 4,
}


def merge_findings(
    snv_variants: Sequence[ClassifiedVariant],
    sv_variants: Sequence[ClassifiedSV],
) -> list[CombinedFinding]:
    """Merge SNV and SV findings into a single ranked list.

    Findings are sorted by tier (Pathogenic first), then by score
    descending within each tier. This produces a unified view for
    clinical reports where the most actionable findings appear first.

    Parameters
    ----------
    snv_variants : Sequence[ClassifiedVariant]
        Point variant results from the standard pipeline.
    sv_variants : Sequence[ClassifiedSV]
        Structural variant results from the SV pipeline.

    Returns
    -------
    list[CombinedFinding]
        Merged findings sorted by clinical priority.
    """
    findings: list[CombinedFinding] = []

    for snv in snv_variants:
        tier = _SNV_TIER.get(snv.classification, 2)
        score = snv.scored.prioritization_score or 0.0
        data = _serialize_snv(snv)
        findings.append(
            CombinedFinding(variant_type="SNV", tier=tier, score=score, data=data)
        )

    for sv in sv_variants:
        tier = _SV_TIER.get(sv.classification, 2)
        score = sv.scored.pathogenicity_score or 0.0
        data = _serialize_sv(sv)
        findings.append(
            CombinedFinding(variant_type="SV", tier=tier, score=score, data=data)
        )

    # Sort: tier ascending (Pathogenic=0 first), then score descending
    findings.sort(key=lambda f: (f.tier, -f.score))
    return findings


def _serialize_snv(variant: ClassifiedVariant) -> dict[str, Any]:
    """Serialize a classified SNV for unified output."""
    v = variant.scored.annotated.variant
    return {
        "type": "SNV",
        "chrom": v.chrom,
        "pos": v.pos,
        "ref": v.ref,
        "alt": v.alt,
        "gene": variant.scored.annotated.gene_name,
        "consequence": variant.scored.annotated.consequence.value,
        "classification": variant.classification.value,
        "evidence_tags": sorted(t.value for t in variant.evidence_tags),
        "score": variant.scored.prioritization_score,
    }


def _serialize_sv(classified: ClassifiedSV) -> dict[str, Any]:
    """Serialize a classified SV for unified output."""
    sv = classified.scored.annotated.sv
    annotated = classified.scored.annotated
    genes = [o.gene_symbol for o in annotated.gene_overlaps[:10]]
    return {
        "type": "SV",
        "chrom": sv.chrom,
        "start": sv.start,
        "end": sv.end,
        "sv_type": sv.sv_type.value,
        "size": sv.length,
        "consequence": annotated.consequence.value,
        "classification": classified.classification.value,
        "genes": genes,
        "genes_affected": annotated.genes_affected,
        "syndrome_name": classified.syndrome_name,
        "evidence_score": classified.evidence_score,
        "score": classified.scored.pathogenicity_score,
    }
