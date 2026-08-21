"""Validate VarTriage against ClinGen Evidence Repository (Expert Panel) variants.

Extracts expert-panel-curated variants from ClinVar, runs VarTriage's ACMG
classifier on each, and computes concordance metrics (sensitivity, PPV, specificity).

Usage:
    python scripts/validate_clingen_erepo.py [--output-dir validation_results/erepo]

Requires:
    - ClinVar VCF at data/references/clinvar.vcf.gz (with CLNREVSTAT field)
    - gnomAD chr-level VCFs (or remote tabix)
    - REVEL scores (optional, improves PP3)
    - GENCODE GTF

Output:
    - erepo_validation_results.json (per-variant classifications)
    - erepo_validation_metrics.json (summary metrics for paper)
"""

from __future__ import annotations

import json
import logging
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = PROJECT_ROOT / "validation_results" / "erepo"


def _get_vartriage_version() -> str:
    """Get the installed vartriage version."""
    try:
        from vartriage import __version__

        return __version__
    except (ImportError, AttributeError):
        return "unknown"


@dataclass
class ErepoVariant:
    """A variant from the ClinGen Evidence Repository."""

    chrom: str
    pos: int
    ref: str
    alt: str
    expert_assertion: str  # Pathogenic, Likely_pathogenic, Benign, Likely_benign, VUS
    review_status: str


def _is_expert_reviewed(rev_str: str) -> bool:
    return "reviewed_by_expert_panel" in rev_str or "practice_guideline" in rev_str


def _normalize_assertion(sig_str: str) -> str | None:
    """Normalize ClinVar CLNSIG into a small set of buckets for metrics.

    Handles composite assertions (e.g. 'Pathogenic/Likely_pathogenic') by
    normalizing each component and returning the single unique value if they
    all agree, or None if components map to different categories.
    """
    if not sig_str:
        return None

    sig_str = sig_str.strip()

    if "Conflicting_interpretations_of_pathogenicity" in sig_str:
        return None

    term_map = {
        "Pathogenic": "Pathogenic",
        "Likely_pathogenic": "Likely_pathogenic",
        "Benign": "Benign",
        "Likely_benign": "Likely_benign",
        "Uncertain_significance": "VUS",
    }

    parts = [sig_str]
    for sep in ("/", ";"):
        if sep in sig_str:
            parts = sig_str.split(sep)
            break

    normalized_components: set[str] = set()
    for raw in parts:
        term = raw.strip()
        if term in term_map:
            normalized_components.add(term_map[term])
            continue
        for key, normalized in term_map.items():
            if term.startswith(key):
                normalized_components.add(normalized)
                break

    if not normalized_components:
        return None
    if len(normalized_components) == 1:
        return next(iter(normalized_components))
    return None


def _parse_vcf_record(rec) -> ErepoVariant | None:
    rev = rec.info.get("CLNREVSTAT")
    if rev is None:
        return None
    rev_str = ",".join(rev) if isinstance(rev, tuple) else str(rev)
    if not _is_expert_reviewed(rev_str):
        return None

    clnsig = rec.info.get("CLNSIG")
    if clnsig is None:
        return None
    sig_str = ",".join(clnsig) if isinstance(clnsig, tuple) else str(clnsig)
    assertion = _normalize_assertion(sig_str)
    if assertion is None or not rec.alts:
        return None

    return ErepoVariant(
        chrom=rec.chrom,
        pos=rec.pos,
        ref=rec.ref,
        alt=rec.alts[0],
        expert_assertion=assertion,
        review_status=rev_str,
    )


def extract_expert_panel_variants(clinvar_vcf: Path) -> list[ErepoVariant]:
    """Extract variants with 'reviewed_by_expert_panel' or 'practice_guideline' status."""
    logger.info("Extracting expert panel variants from %s...", clinvar_vcf)

    vcf = pysam.VariantFile(str(clinvar_vcf))
    variants = [v for rec in vcf.fetch() if (v := _parse_vcf_record(rec)) is not None]
    vcf.close()

    logger.info("Extracted %d expert panel variants", len(variants))
    dist = Counter(v.expert_assertion for v in variants)
    for k, v in sorted(dist.items()):
        logger.info("  %s: %d", k, v)

    return variants


def classify_variant(variant: ErepoVariant, classifier) -> dict:
    """Classify a single variant using VarTriage's ACMG classifier.

    Returns a dict with classification result and evidence tags.
    """
    from vartriage.models import Variant

    v = Variant(
        chrom=variant.chrom,
        pos=variant.pos,
        id=None,
        ref=variant.ref,
        alt=variant.alt,
        qual=None,
        filter_status=".",
    )

    try:
        result = classifier.classify(v)
        return {
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": variant.ref,
            "alt": variant.alt,
            "expert_assertion": variant.expert_assertion,
            "vartriage_classification": result.classification.value,
            "evidence_tags": [t.name for t in result.evidence_tags],
            "success": True,
        }
    except Exception as e:
        return {
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": variant.ref,
            "alt": variant.alt,
            "expert_assertion": variant.expert_assertion,
            "vartriage_classification": "ERROR",
            "evidence_tags": [],
            "error": str(e),
            "success": False,
        }


def compute_metrics(results: list[dict]) -> dict:
    """Compute concordance metrics from classification results."""
    successful = [r for r in results if r["success"]]

    # Expert P/LP
    expert_path = [
        r
        for r in successful
        if r["expert_assertion"] in ("Pathogenic", "Likely_pathogenic")
    ]
    expert_benign = [
        r for r in successful if r["expert_assertion"] in ("Benign", "Likely_benign")
    ]
    expert_vus = [r for r in successful if r["expert_assertion"] == "VUS"]

    # VarTriage P/LP
    vt_path = [
        r
        for r in expert_path
        if r["vartriage_classification"] in ("Pathogenic", "Likely_Pathogenic")
    ]
    vt_benign_correct = [
        r
        for r in expert_benign
        if r["vartriage_classification"] in ("Benign", "Likely_Benign")
    ]

    # False positives: expert B/LB but VarTriage says P/LP
    false_pos_path = [
        r
        for r in expert_benign
        if r["vartriage_classification"] in ("Pathogenic", "Likely_Pathogenic")
    ]
    # False positives benign: expert P/LP but VarTriage says B/LB
    false_pos_benign = [
        r
        for r in expert_path
        if r["vartriage_classification"] in ("Benign", "Likely_Benign")
    ]

    path_sensitivity = len(vt_path) / len(expert_path) if expert_path else 0
    path_ppv = (
        len(vt_path) / (len(vt_path) + len(false_pos_path))
        if (len(vt_path) + len(false_pos_path)) > 0
        else 0
    )
    benign_sensitivity = (
        len(vt_benign_correct) / len(expert_benign) if expert_benign else 0
    )
    benign_ppv = (
        len(vt_benign_correct) / (len(vt_benign_correct) + len(false_pos_benign))
        if (len(vt_benign_correct) + len(false_pos_benign)) > 0
        else 0
    )

    # Evidence tag distribution
    tag_dist = Counter()
    for r in successful:
        for tag in r["evidence_tags"]:
            tag_dist[tag] += 1

    return {
        "total_variants": len(results),
        "successful_classifications": len(successful),
        "failed_classifications": len(results) - len(successful),
        "expert_pathogenic_count": len(expert_path),
        "expert_benign_count": len(expert_benign),
        "expert_vus_count": len(expert_vus),
        "pathogenic_sensitivity": round(path_sensitivity, 4),
        "pathogenic_ppv": round(path_ppv, 4),
        "pathogenic_true_positives": len(vt_path),
        "pathogenic_false_negatives": len(expert_path) - len(vt_path),
        "pathogenic_false_positives": len(false_pos_path),
        "benign_sensitivity": round(benign_sensitivity, 4),
        "benign_ppv": round(benign_ppv, 4),
        "benign_true_positives": len(vt_benign_correct),
        "benign_false_positives_from_path": len(false_pos_benign),
        "evidence_tag_distribution": dict(tag_dist.most_common()),
        "classification_distribution": dict(
            Counter(r["vartriage_classification"] for r in successful)
        ),
    }


def main():
    output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    clinvar_vcf = PROJECT_ROOT / "data" / "references" / "clinvar.vcf.gz"
    if not clinvar_vcf.exists():
        logger.error("ClinVar VCF not found at %s", clinvar_vcf)
        sys.exit(1)

    # Step 1: Extract expert panel variants
    variants = extract_expert_panel_variants(clinvar_vcf)

    # Step 2: Initialize VarTriage classifier
    logger.info("Initializing VarTriage ACMG classifier...")
    from vartriage.classification.acmg import ACMGClassifier

    classifier = ACMGClassifier()
    logger.info("Classifier initialized")

    # Step 3: Classify each variant
    logger.info("Classifying %d variants...", len(variants))
    results = []
    start = time.perf_counter()

    for i, variant in enumerate(variants):
        result = classify_variant(variant, classifier)
        results.append(result)

        if (i + 1) % 1000 == 0:
            elapsed = time.perf_counter() - start
            rate = (i + 1) / elapsed
            logger.info(
                "  Progress: %d/%d (%.1f variants/sec)", i + 1, len(variants), rate
            )

    elapsed = time.perf_counter() - start
    logger.info(
        "Classification complete: %d variants in %.1f seconds (%.1f/sec)",
        len(results),
        elapsed,
        len(results) / elapsed,
    )

    # Step 4: Compute metrics
    metrics = compute_metrics(results)
    metrics["runtime_seconds"] = round(elapsed, 1)
    metrics["vartriage_version"] = _get_vartriage_version()

    # Step 5: Save outputs
    results_path = output_dir / "erepo_validation_results.json"
    metrics_path = output_dir / "erepo_validation_metrics.json"

    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    # Print summary
    logger.info("=" * 60)
    logger.info("eRepo VALIDATION RESULTS")
    logger.info("=" * 60)
    logger.info("Total variants: %d", metrics["total_variants"])
    logger.info("Successful: %d", metrics["successful_classifications"])
    logger.info("Expert P/LP: %d", metrics["expert_pathogenic_count"])
    logger.info("Expert B/LB: %d", metrics["expert_benign_count"])
    logger.info(
        "Pathogenic sensitivity: %.1f%%", metrics["pathogenic_sensitivity"] * 100
    )
    logger.info("Pathogenic PPV: %.1f%%", metrics["pathogenic_ppv"] * 100)
    logger.info("Benign sensitivity: %.1f%%", metrics["benign_sensitivity"] * 100)
    logger.info("Benign PPV: %.1f%%", metrics["benign_ppv"] * 100)
    logger.info("Results: %s", results_path)
    logger.info("Metrics: %s", metrics_path)


if __name__ == "__main__":
    main()
