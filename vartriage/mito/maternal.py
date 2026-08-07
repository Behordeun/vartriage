"""Maternal inheritance verification for mitochondrial variants.

When trio data is available (proband + mother + father), verifies that
mtDNA variants follow maternal inheritance:
- Present in mother = expected (maternal transmission)
- Absent in father = expected (mtDNA is maternally inherited)
- Present in proband but absent in mother = potential de novo mtDNA mutation

De novo mitochondrial mutations are rare but clinically significant.
Their presence strengthens pathogenicity evidence for novel variants.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Literal

from vartriage.mito.classifier import MitoClassifiedVariant
from vartriage.models.variant import Variant

logger = logging.getLogger(__name__)

MaternalStatus = Literal[
    "maternal",
    "de_novo",
    "paternal_unexpected",
    "unknown",
]


@dataclass(frozen=True, slots=True)
class MaternalInheritanceResult:
    """Result of maternal inheritance check for a mtDNA variant.

    Parameters
    ----------
    status
        Classification of the inheritance pattern.
    in_mother
        Whether the variant allele is present in the mother's genotype.
    in_father
        Whether the variant allele is present in the father's genotype.
        Should always be False for true mtDNA variants.
    note
        Human-readable explanation of the finding.
    """

    status: MaternalStatus
    in_mother: bool
    in_father: bool
    note: str


def check_maternal_inheritance(
    variant: Variant,
    proband_name: str,
    mother_name: str,
    father_name: str,
) -> MaternalInheritanceResult | None:
    """Check maternal inheritance pattern for a mitochondrial variant.

    Examines trio genotype data in the variant's info dict to determine
    whether the mtDNA variant follows expected maternal transmission.

    Parameters
    ----------
    variant
        Raw Variant with _pysam_samples in info (from VCFParser with
        extract_samples=True).
    proband_name
        Sample name of the proband/child.
    mother_name
        Sample name of the mother.
    father_name
        Sample name of the father.

    Returns
    -------
    MaternalInheritanceResult or None
        Result of the check, or None if trio genotype data is unavailable.
    """
    sample_data = variant.info.get("_pysam_samples")
    if sample_data is None:
        return None

    mother_gt = _extract_gt_string(sample_data, mother_name)
    father_gt = _extract_gt_string(sample_data, father_name)

    if mother_gt is None and father_gt is None:
        return None

    in_mother = _has_alt_allele(mother_gt) if mother_gt is not None else False
    in_father = _has_alt_allele(father_gt) if father_gt is not None else False

    status = _classify_status(in_mother, in_father)
    note = _build_note(status)

    return MaternalInheritanceResult(
        status=status,
        in_mother=in_mother,
        in_father=in_father,
        note=note,
    )


def annotate_maternal_inheritance(
    results: list[MitoClassifiedVariant],
    proband_name: str,
    mother_name: str,
    father_name: str,
) -> list[tuple[MitoClassifiedVariant, MaternalInheritanceResult | None]]:
    """Annotate a list of classified mito variants with maternal inheritance info.

    Parameters
    ----------
    results
        Classified mitochondrial variants from the MitochondrialPipeline.
    proband_name
        Proband sample name.
    mother_name
        Mother sample name.
    father_name
        Father sample name.

    Returns
    -------
    list of (MitoClassifiedVariant, MaternalInheritanceResult | None)
        Each classified variant paired with its maternal inheritance check,
        or None if trio data was unavailable for that variant.
    """
    annotated: list[tuple[MitoClassifiedVariant, MaternalInheritanceResult | None]] = []
    de_novo_count = 0

    for classified in results:
        inheritance = check_maternal_inheritance(
            classified.variant,
            proband_name,
            mother_name,
            father_name,
        )
        annotated.append((classified, inheritance))

        if inheritance is not None and inheritance.status == "de_novo":
            de_novo_count += 1
            logger.warning(
                "Potential de novo mtDNA mutation at pos %d (%s>%s): absent in mother",
                classified.variant.pos,
                classified.variant.ref,
                classified.variant.alt,
            )

    if de_novo_count > 0:
        logger.info(
            "Maternal inheritance check: %d potential de novo mtDNA mutations",
            de_novo_count,
        )

    return annotated


def _extract_gt_string(sample_data: dict[str, Any], sample_name: str) -> str | None:
    """Extract genotype string from pysam sample data."""
    entry = sample_data.get(sample_name, {})
    gt = entry.get("GT")
    if gt is None:
        return None

    # Format pysam GT tuple into string
    parts = []
    for allele in gt:
        parts.append(str(allele) if allele is not None else ".")
    return "/".join(parts)


def _has_alt_allele(gt_string: str) -> bool:
    """Check if a genotype string contains any alt allele (non-0, non-.)."""
    for allele in gt_string.replace("|", "/").split("/"):
        stripped = allele.strip()
        if stripped not in (".", "0", ""):
            return True
    return False


def _classify_status(in_mother: bool, in_father: bool) -> MaternalStatus:
    """Determine maternal inheritance status from parental genotypes."""
    if in_mother and not in_father:
        return "maternal"
    if not in_mother and not in_father:
        return "de_novo"
    if in_father:
        # mtDNA in father is unexpected (contamination, NUMTs, or error)
        return "paternal_unexpected"
    return "unknown"


def _build_note(status: MaternalStatus) -> str:
    """Generate a human-readable note for the inheritance result."""
    if status == "maternal":
        return "Maternally inherited (present in mother, absent in father)"
    if status == "de_novo":
        return (
            "Potential de novo mtDNA mutation (absent in both parents). "
            "Rare but clinically significant if confirmed."
        )
    if status == "paternal_unexpected":
        return (
            "Unexpected: variant detected in father. "
            "May indicate NUMTs, sample contamination, or data error."
        )
    return "Inheritance pattern could not be determined"
