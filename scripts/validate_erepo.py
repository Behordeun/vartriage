#!/usr/bin/env python3
"""End-to-end ClinGen Expert Panel (eRepo) validation for VarTriage.

Extracts expert-panel variants from ClinVar, filters REVEL scores using
an efficient set-based lookup (fixing the coverage gap from the original
awk approach), creates a synthetic VCF, runs the VarTriage pipeline, and
computes concordance metrics.

Requirements:
    - vartriage >= 0.17.2 (pip install vartriage)
    - pysam
    - ClinVar VCF: data/references/clinvar.vcf.gz (+ .tbi)
    - REVEL genome-wide: data/references/revel_genome_wide.tsv
    - GENCODE: data/references/gencode.v46.annotation.gtf
    - Internet access for remote gnomAD queries

Usage:
    cd /path/to/vartriage
    .venv/bin/python3 scripts/validate_erepo.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import pysam

# --- Configuration -----------------------------------------------------------

VARTRIAGE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = VARTRIAGE_DIR / "data" / "references"
OUTPUT_DIR = VARTRIAGE_DIR / "validation_results" / "erepo"

CLINVAR_VCF = DATA_DIR / "clinvar.vcf.gz"
REVEL_TSV = DATA_DIR / "revel_genome_wide.tsv"
GENCODE_GTF = DATA_DIR / "gencode.v46.annotation.gtf"
CLINVAR_SIMPLE = DATA_DIR / "clinvar.tsv"

EREPO_POSITIONS = OUTPUT_DIR / "erepo_positions.txt"
REVEL_FILTERED = OUTPUT_DIR / "revel_erepo.tsv"
EREPO_VCF = OUTPUT_DIR / "erepo_variants.vcf.gz"
CLASSIFICATIONS_JSON = OUTPUT_DIR / "erepo_classifications.json"
METRICS_JSON = OUTPUT_DIR / "erepo_validation_metrics.json"


def step(msg: str) -> None:
    print(f"\n{'=' * 60}\n{msg}\n{'=' * 60}")


# --- Step 1: Extract expert-panel variants from ClinVar ---------------------


def extract_erepo_variants() -> list[tuple[str, int, str, str, str]]:
    """Return (chrom, pos, ref, alt, clnsig) for expert-panel variants."""
    step("Step 1: Extracting expert-panel variants from ClinVar VCF")

    vcf = pysam.VariantFile(str(CLINVAR_VCF))
    variants: list[tuple[str, int, str, str, str]] = []

    for rec in vcf.fetch():
        rev = rec.info.get("CLNREVSTAT")
        if rev is None:
            continue
        rev_str = ",".join(rev) if isinstance(rev, tuple) else str(rev)
        if (
            "reviewed_by_expert_panel" not in rev_str
            and "practice_guideline" not in rev_str
        ):
            continue
        if rec.alts is None or len(rec.alts) == 0:
            continue
        clnsig = rec.info.get("CLNSIG")
        if clnsig is None:
            continue
        sig_str = ",".join(clnsig) if isinstance(clnsig, tuple) else str(clnsig)
        chrom = rec.chrom if rec.chrom.startswith("chr") else f"chr{rec.chrom}"
        variants.append((chrom, rec.pos, rec.ref, rec.alts[0], sig_str))

    vcf.close()

    # Write positions file
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(EREPO_POSITIONS, "w") as f:
        for chrom, pos, ref, alt, _ in variants:
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

    print(f"  Extracted {len(variants)} expert-panel variants")
    snvs = sum(1 for _, _, r, a, _ in variants if len(r) == 1 and len(a) == 1)
    indels = len(variants) - snvs
    print(f"  SNVs: {snvs}, Indels: {indels}")
    return variants


# --- Step 2: Filter REVEL scores (THE FIX) ----------------------------------


def _compute_revel_fingerprint(
    variants: list[tuple[str, int, str, str, str]],
) -> str:
    """Hash variant keys + REVEL source mtime for cache invalidation."""
    import hashlib

    h = hashlib.sha256()
    # Hash the SNV key set (sorted for determinism)
    snv_keys = sorted(
        f"{c}:{p}:{r}:{a}" for c, p, r, a, _ in variants if len(r) == 1 and len(a) == 1
    )
    h.update(str(len(snv_keys)).encode())
    for k in snv_keys[:100]:  # sample first 100 for speed
        h.update(k.encode())
    h.update(snv_keys[-1].encode() if snv_keys else b"")
    # Include REVEL source file mtime
    if REVEL_TSV.exists():
        h.update(str(int(REVEL_TSV.stat().st_mtime)).encode())
    return h.hexdigest()[:16]


def filter_revel(
    variants: list[tuple[str, int, str, str, str]], force: bool = False
) -> int:
    """Filter genome-wide REVEL to eRepo positions using set-based lookup.

    This fixes the coverage gap from the original awk-based approach.
    The awk command scanned 80M rows but was correct -- the real issue is
    that REVEL only scores missense variants (~53% of eRepo SNVs).

    This Python implementation uses a hash set for O(1) lookups and streams
    the 80M-row file once. Takes ~10-15 minutes on the 1.9GB file.

    If the filtered file already exists and is non-empty, it is reused
    (the REVEL reference data doesn't change). Pass force=True to re-filter.
    """
    step("Step 2: Filtering REVEL scores for eRepo positions")

    # Compute fingerprint early so it's available for both reuse check and write
    fingerprint_file = REVEL_FILTERED.with_suffix(".fingerprint")
    current_fingerprint = _compute_revel_fingerprint(variants)

    # Reuse existing filtered file only if inputs haven't changed
    if (
        not force
        and REVEL_FILTERED.exists()
        and REVEL_FILTERED.stat().st_size > 100
        and fingerprint_file.exists()
        and fingerprint_file.read_text().strip() == current_fingerprint
    ):
        with open(REVEL_FILTERED) as f:
            line_count = sum(1 for _ in f) - 1  # minus header
        print(f"  Reusing existing {REVEL_FILTERED.name} ({line_count} scores)")
        print("  (pass --force to re-filter from genome-wide REVEL)")

        snv_count = sum(1 for _, _, r, a, _ in variants if len(r) == 1 and len(a) == 1)
        coverage_pct = line_count / snv_count * 100 if snv_count else 0
        print(f"  Coverage: {line_count} / {snv_count} SNVs ({coverage_pct:.1f}%)")
        print("  (REVEL only scores missense variants)")
        return line_count

    # Build lookup set of eRepo SNV keys: (chrom, pos, ref, alt)
    snv_keys: set[tuple[str, str, str, str]] = set()
    for chrom, pos, ref, alt, _ in variants:
        if len(ref) == 1 and len(alt) == 1:
            snv_keys.add((chrom, str(pos), ref, alt))

    print(f"  eRepo SNVs to look up: {len(snv_keys)}")
    print(f"  Streaming {REVEL_TSV.name} (1.9 GB, ~80M rows)...")
    print("  This takes 10-15 minutes. Use existing file next time (--no-force).")

    # Stream REVEL TSV, match against set
    matched = 0
    with open(REVEL_TSV) as fin, open(REVEL_FILTERED, "w") as fout:
        fout.write("chrom\tpos\tref\talt\tscore\n")
        fin.readline()  # skip header
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            key = (parts[0], parts[1], parts[2], parts[3])
            if key in snv_keys:
                fout.write(line)
                matched += 1

    coverage_pct = matched / len(snv_keys) * 100 if snv_keys else 0
    print(f"  REVEL matches: {matched} / {len(snv_keys)} SNVs ({coverage_pct:.1f}%)")
    print("  (REVEL only scores missense variants -- remaining are")
    print("   synonymous, nonsense, splice-site, or non-coding)")
    fingerprint_file.write_text(current_fingerprint + "\n")
    return matched


# --- Step 3: Create synthetic VCF -------------------------------------------


def create_erepo_vcf(variants: list[tuple[str, int, str, str, str]]) -> None:
    """Create synthetic VCF with QUAL=99 to pass the quality filter.

    VarTriage's QualityFilter drops variants with missing QUAL (None).
    Since expert-panel variants are curated (not machine-called), we assign
    QUAL=99 as a pass-through value -- this is standard practice for
    synthetic validation VCFs.
    """
    step("Step 3: Creating eRepo synthetic VCF")

    # Exclude chrM variants -- they use a separate mito pipeline and
    # would cause the JSON writer to use structured output format
    nuclear_variants = [v for v in variants if not v[0].startswith("chrM")]
    mito_skipped = len(variants) - len(nuclear_variants)

    header = pysam.VariantHeader()
    header.add_sample("EREPO")
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    contigs: set[str] = set()
    for chrom, _, _, _, _ in nuclear_variants:
        contigs.add(chrom)

    for c in sorted(contigs):
        header.add_line(f"##contig=<ID={c},length=300000000>")

    output_path = str(EREPO_VCF)
    with pysam.VariantFile(output_path, "wz", header=header) as out:
        for chrom, pos, ref, alt, _ in sorted(
            nuclear_variants, key=lambda x: (x[0], x[1])
        ):
            rec = out.new_record()
            rec.contig = chrom
            rec.pos = pos
            rec.alleles = (ref, alt)
            rec.qual = 99
            rec.samples["EREPO"]["GT"] = (0, 1)
            out.write(rec)

    pysam.tabix_index(output_path, preset="vcf", force=True)
    print(
        f"  Created VCF with {len(nuclear_variants)} nuclear variants (skipped {mito_skipped} chrM)"
    )


# --- Step 4: Run VarTriage pipeline -----------------------------------------


def run_pipeline() -> None:
    step("Step 4: Running VarTriage pipeline")

    vartriage_bin = VARTRIAGE_DIR / ".venv" / "bin" / "vartriage"
    cmd = [
        str(vartriage_bin),
        "--vcf",
        str(EREPO_VCF),
        "--gene-annotation",
        str(GENCODE_GTF),
        "--clinvar",
        str(CLINVAR_SIMPLE),
        "--revel-scores",
        str(REVEL_FILTERED),
        "--gnomad-remote",
        "gnomad-exomes-v4-grch38",
        "--output",
        str(CLASSIFICATIONS_JSON),
        "--output-format",
        "json",
        "--skip-mito",
        "--no-confirm",
    ]

    print("  Command: vartriage --vcf ... --gnomad-remote ... --skip-mito --no-confirm")
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        print(f"  FAILED (exit {result.returncode})")
        print(f"  stderr: {result.stderr[:1000]}")
        sys.exit(1)

    print(f"  Completed in {elapsed:.1f}s")


# --- Step 5: Compute concordance metrics ------------------------------------


def normalize_assertion(sig_str: str) -> str | None:
    """Normalize ClinVar CLNSIG to a 5-tier classification."""
    sig = sig_str.lower().replace("_", " ").replace("/", ",")

    # Handle composite assertions
    if "pathogenic" in sig and "benign" in sig:
        return None  # conflicting
    if "pathogenic" in sig and "likely" not in sig:
        return "Pathogenic"
    if "likely pathogenic" in sig or "likely_pathogenic" in sig:
        return "Likely_pathogenic"
    if "benign" in sig and "likely" not in sig:
        return "Benign"
    if "likely benign" in sig or "likely_benign" in sig:
        return "Likely_benign"
    if "uncertain" in sig:
        return "VUS"
    return None


def compute_metrics(
    variants: list[tuple[str, int, str, str, str]],
) -> dict:
    step("Step 5: Computing concordance metrics")

    # Load pipeline output
    raw = json.loads(CLASSIFICATIONS_JSON.read_text())
    # Handle both formats: flat list (--skip-mito) or structured dict
    if isinstance(raw, list):
        data = raw
    elif isinstance(raw, dict):
        data = raw.get("variants", [])
    else:
        print(f"  ERROR: unexpected output format: {type(raw)}")
        sys.exit(1)
    print(f"  Pipeline classified: {len(data)} variants")

    # Build expert assertion lookup
    expert: dict[tuple[str, int, str, str], str] = {}
    for chrom, pos, ref, alt, sig_str in variants:
        assertion = normalize_assertion(sig_str)
        if assertion:
            expert[(chrom, pos, ref, alt)] = assertion

    # Match pipeline output to expert assertions
    matched: list[dict] = []
    for v in data:
        key = (
            v.get("chromosome", ""),
            v.get("position", 0),
            v.get("ref_allele", ""),
            v.get("alt_allele", ""),
        )
        if key in expert:
            matched.append(
                {
                    "vt": v.get("acmg_classification", "Unknown"),
                    "expert": expert[key],
                    "tags": v.get("evidence_tags", []),
                    "cons": v.get("functional_consequence", ""),
                }
            )

    print(f"  Matched to expert assertions: {len(matched)}")

    # Core metrics
    exp_path = [
        m for m in matched if m["expert"] in ("Pathogenic", "Likely_pathogenic")
    ]
    exp_benign = [m for m in matched if m["expert"] in ("Benign", "Likely_benign")]
    exp_vus = [m for m in matched if m["expert"] == "VUS"]

    tp = [m for m in exp_path if m["vt"] in ("Pathogenic", "Likely_Pathogenic")]
    fp = [m for m in exp_benign if m["vt"] in ("Pathogenic", "Likely_Pathogenic")]
    bn_correct = [m for m in exp_benign if m["vt"] in ("Benign", "Likely_Benign")]
    fn = [m for m in exp_path if m["vt"] not in ("Pathogenic", "Likely_Pathogenic")]

    sens = len(tp) / len(exp_path) if exp_path else 0
    ppv = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) else 0
    bn_sens = len(bn_correct) / len(exp_benign) if exp_benign else 0

    # Per-consequence breakdown
    cons_stats: dict[str, dict] = {}
    for cons in sorted(set(m["cons"] for m in exp_path)):
        cv = [m for m in exp_path if m["cons"] == cons]
        ct = [m for m in cv if m["vt"] in ("Pathogenic", "Likely_Pathogenic")]
        if cv:
            cons_stats[cons] = {
                "total": len(cv),
                "tp": len(ct),
                "sensitivity": round(len(ct) / len(cv), 4),
            }

    # VarTriage classification distribution for P/LP expert variants
    vt_dist = Counter(m["vt"] for m in exp_path)

    # Evidence tag distribution
    tag_dist = Counter(t for m in matched for t in m["tags"])

    # REVEL coverage among missense P/LP
    missense_path = [m for m in exp_path if m["cons"] == "Missense"]
    missense_with_pp3 = [
        m
        for m in missense_path
        if "PP3" in m["tags"]
        or "PP3_MODERATE" in m["tags"]
        or "PP3_Moderate" in m["tags"]
    ]

    metrics = {
        "vartriage_version": "0.17.2",
        "dataset": "ClinGen Expert Panel (ClinVar CLNREVSTAT)",
        "total_erepo_variants": len(variants),
        "total_classified": len(data),
        "total_matched": len(matched),
        "expert_pathogenic_lp": len(exp_path),
        "expert_benign_lb": len(exp_benign),
        "expert_vus": len(exp_vus),
        "pathogenic_sensitivity": round(sens, 4),
        "pathogenic_ppv": round(ppv, 4),
        "pathogenic_tp": len(tp),
        "pathogenic_fp": len(fp),
        "pathogenic_fn": len(fn),
        "benign_sensitivity": round(bn_sens, 4),
        "benign_tp": len(bn_correct),
        "per_consequence": cons_stats,
        "vt_classification_of_expert_plp": dict(vt_dist.most_common()),
        "top_evidence_tags": dict(tag_dist.most_common(15)),
        "missense_plp_total": len(missense_path),
        "missense_plp_with_pp3": len(missense_with_pp3),
        "revel_coverage_note": (
            "REVEL scores ~53% of eRepo SNVs. Coverage gap is inherent: "
            "REVEL only predicts pathogenicity for missense variants. "
            "Nonsense, frameshift, splice-site, and synonymous variants "
            "are not scored by REVEL."
        ),
    }

    METRICS_JSON.write_text(json.dumps(metrics, indent=2))

    # Print summary
    print(f"\n{'=' * 60}")
    print("  eRepo VALIDATION RESULTS (VarTriage v0.17.2)")
    print(f"{'=' * 60}")
    print(f"  Dataset: {len(variants)} expert-panel variants from ClinVar")
    print(f"  Classified: {len(data)} | Matched: {len(matched)}")
    print("")
    print(f"  Expert P/LP: {len(exp_path)}")
    print(f"  Expert B/LB: {len(exp_benign)}")
    print(f"  Expert VUS:  {len(exp_vus)}")
    print("")
    print(f"  PATHOGENIC SENSITIVITY: {sens * 100:.1f}% ({len(tp)}/{len(exp_path)})")
    print(f"  PATHOGENIC PPV:         {ppv * 100:.1f}% ({len(tp)}/{len(tp) + len(fp)})")
    print(
        f"  BENIGN SENSITIVITY:     {bn_sens * 100:.1f}% ({len(bn_correct)}/{len(exp_benign)})"
    )
    print("")
    print("  Per-consequence sensitivity (P/LP):")
    for cons, s in sorted(cons_stats.items(), key=lambda x: -x[1]["total"]):
        print(f"    {cons:30s} {s['total']:5d} variants  {s['sensitivity'] * 100:.1f}%")
    print("")
    print("  VarTriage calls for expert P/LP variants:")
    for cls, count in vt_dist.most_common():
        print(f"    {cls:20s} {count:5d}  ({count / len(exp_path) * 100:.1f}%)")
    print("")
    print("  REVEL coverage:")
    print(
        f"    Missense P/LP with PP3/PP3_Moderate: {len(missense_with_pp3)}/{len(missense_path)}"
    )
    if missense_path:
        miss_tp = [
            m for m in missense_path if m["vt"] in ("Pathogenic", "Likely_Pathogenic")
        ]
        print(
            f"    Missense P/LP sensitivity: {len(miss_tp) / len(missense_path) * 100:.1f}%"
        )
    print("")
    print(f"  Metrics saved: {METRICS_JSON}")
    print(f"{'=' * 60}")

    return metrics


# --- Main --------------------------------------------------------------------


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="VarTriage eRepo validation")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-filter REVEL from genome-wide file (slow, ~15 min)",
    )
    cli_args = parser.parse_args()

    t0 = time.time()
    print("VarTriage eRepo Validation")
    print(f"Working directory: {VARTRIAGE_DIR}")
    print(f"Output directory:  {OUTPUT_DIR}")

    # Verify prerequisites
    for path, label in [
        (CLINVAR_VCF, "ClinVar VCF"),
        (REVEL_TSV, "REVEL genome-wide TSV"),
        (GENCODE_GTF, "GENCODE GTF"),
        (CLINVAR_SIMPLE, "ClinVar simple TSV"),
    ]:
        if not path.exists():
            print(f"ERROR: {label} not found at {path}")
            sys.exit(1)

    variants = extract_erepo_variants()
    filter_revel(variants, force=cli_args.force)
    create_erepo_vcf(variants)
    run_pipeline()
    compute_metrics(variants)

    elapsed = time.time() - t0
    print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
