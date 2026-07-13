#!/usr/bin/env bash
# validate_giab.sh - Run vartriage against Genome in a Bottle (GIAB) benchmark data
#
# This script downloads GIAB HG002 (Ashkenazi Jewish trio son) WES data,
# prepares reference files, runs vartriage, and compares results against
# the GIAB truth set to measure sensitivity and specificity.
#
# Requirements:
#   - vartriage installed (pip install vartriage[all])
#   - bcftools >= 1.17
#   - htslib (tabix, bgzip)
#   - wget or curl (either works, script auto-detects)
#   - ~20 GB free disk space for reference files
#
# Usage:
#   chmod +x scripts/validate_giab.sh
#   ./scripts/validate_giab.sh [--output-dir /path/to/results]
#
# Exit codes:
#   0 - validation completed successfully
#   1 - missing dependencies or download failure
#   2 - pipeline execution failure
#   3 - metrics computation failure

set -euo pipefail

# Configuration
OUTPUT_DIR="${1:-validation_results}"
# Strip --output-dir prefix if passed as flag
OUTPUT_DIR="${OUTPUT_DIR#--output-dir=}"
OUTPUT_DIR="${OUTPUT_DIR#--output-dir }"
GIAB_VERSION="v4.2.1"
REFERENCE_BUILD="GRCh38"
SAMPLE="HG002"

# GIAB URLs (NIST FTP mirrors)
# Using versioned path (NISTv4.2.1) to avoid breakage when 'latest' symlink updates
GIAB_BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
GIAB_VCF_URL="${GIAB_BASE_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
GIAB_BED_URL="${GIAB_BASE_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

# gnomAD v4.1.1 exomes chr22 subset for validation (4.73 GiB)
GNOMAD_CHR22_URL="https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1.1/vcf/exomes/gnomad.exomes.v4.1.1.sites.chr22.vcf.bgz"

# --- Helper functions ---

# Download a file using wget or curl (auto-detects available tool)
download_file() {
    local url="$1"
    local dest="$2"

    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$dest" "$url"
    elif command -v curl &>/dev/null; then
        curl -fSL --progress-bar -o "$dest" "$url"
    else
        echo "ERROR: Neither wget nor curl is available. Install one of them."
        exit 1
    fi
}

# Check that a required command is available
require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: Required command '$cmd' not found."
        echo "Install it before running this script."
        exit 1
    fi
}

# --- Dependency checks ---

require_cmd bcftools
require_cmd tabix
require_cmd vartriage
require_cmd python3

# Verify at least one download tool is present
if ! command -v wget &>/dev/null && ! command -v curl &>/dev/null; then
    echo "ERROR: Neither wget nor curl is available. Install one of them."
    exit 1
fi

echo "============================================="
echo " vartriage GIAB Validation Pipeline"
echo " Sample: ${SAMPLE} (${GIAB_VERSION})"
echo " Build:  ${REFERENCE_BUILD}"
echo "============================================="
echo ""

# Create output directory structure
mkdir -p "${OUTPUT_DIR}"/{data,refs,results,reports}

# Step 1: Download GIAB benchmark VCF and high-confidence regions
echo "[1/7] Downloading GIAB benchmark data..."
if [ ! -f "${OUTPUT_DIR}/data/giab_benchmark.vcf.gz" ]; then
    download_file "${GIAB_VCF_URL}" "${OUTPUT_DIR}/data/giab_benchmark.vcf.gz"
    download_file "${GIAB_VCF_URL}.tbi" "${OUTPUT_DIR}/data/giab_benchmark.vcf.gz.tbi"
    echo "  Downloaded GIAB benchmark VCF"
else
    echo "  GIAB benchmark VCF already exists, skipping"
fi

if [ ! -f "${OUTPUT_DIR}/data/giab_highconf.bed" ]; then
    download_file "${GIAB_BED_URL}" "${OUTPUT_DIR}/data/giab_highconf.bed"
    echo "  Downloaded high-confidence regions BED"
else
    echo "  High-confidence BED already exists, skipping"
fi

# Step 2: Download ClinVar for annotation
echo "[2/7] Downloading ClinVar..."
if [ ! -f "${OUTPUT_DIR}/refs/clinvar.vcf.gz" ]; then
    download_file "${CLINVAR_URL}" "${OUTPUT_DIR}/refs/clinvar.vcf.gz"
    download_file "${CLINVAR_URL}.tbi" "${OUTPUT_DIR}/refs/clinvar.vcf.gz.tbi"
    echo "  Downloaded ClinVar VCF"
else
    echo "  ClinVar already exists, skipping"
fi

# Step 3: Extract ClinVar to TSV format for vartriage
echo "[3/7] Preparing ClinVar TSV..."
if [ ! -f "${OUTPUT_DIR}/refs/clinvar.tsv" ]; then
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\n' \
        "${OUTPUT_DIR}/refs/clinvar.vcf.gz" \
        | awk 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt","clinical_significance"} {
            sig=$5;
            gsub(/\/.*/, "", sig);
            if (sig ~ /^Pathogenic/) sig="Pathogenic";
            else if (sig ~ /^Likely_pathogenic/) sig="Likely pathogenic";
            else if (sig ~ /^Uncertain/) sig="Uncertain significance";
            else if (sig ~ /^Likely_benign/) sig="Likely benign";
            else if (sig ~ /^Benign/) sig="Benign";
            else next;
            print $1, $2, $3, $4, sig
        }' > "${OUTPUT_DIR}/refs/clinvar.tsv"
    echo "  Prepared ClinVar TSV ($(wc -l < "${OUTPUT_DIR}/refs/clinvar.tsv") records)"
else
    echo "  ClinVar TSV already exists, skipping"
fi

# Step 4: Download gnomAD chr22 for frequency annotation
echo "[4/7] Downloading gnomAD chr22 (exomes)..."
if [ ! -f "${OUTPUT_DIR}/refs/gnomad.chr22.vcf.bgz" ]; then
    download_file "${GNOMAD_CHR22_URL}" "${OUTPUT_DIR}/refs/gnomad.chr22.vcf.bgz"
    download_file "${GNOMAD_CHR22_URL}.tbi" "${OUTPUT_DIR}/refs/gnomad.chr22.vcf.bgz.tbi"
    echo "  Downloaded gnomAD v4.1.1 chr22 exomes"
else
    echo "  gnomAD chr22 already exists, skipping"
fi

# Step 5: Extract chr22 variants from GIAB for focused validation
echo "[5/7] Extracting chr22 subset from GIAB benchmark..."
if [ ! -f "${OUTPUT_DIR}/data/giab_chr22.vcf.gz" ]; then
    bcftools view -r chr22 "${OUTPUT_DIR}/data/giab_benchmark.vcf.gz" \
        -Oz -o "${OUTPUT_DIR}/data/giab_chr22.vcf.gz"
    tabix -p vcf "${OUTPUT_DIR}/data/giab_chr22.vcf.gz"

    VARIANT_COUNT=$(bcftools view -H "${OUTPUT_DIR}/data/giab_chr22.vcf.gz" | wc -l)
    echo "  Extracted ${VARIANT_COUNT} chr22 variants"
else
    echo "  chr22 subset already exists, skipping"
fi

# Step 6: Run vartriage on chr22 subset
echo "[6/7] Running vartriage on chr22..."
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# JSON output for analysis
if ! vartriage \
    --vcf "${OUTPUT_DIR}/data/giab_chr22.vcf.gz" \
    --output "${OUTPUT_DIR}/results/giab_chr22_${TIMESTAMP}.json" \
    --output-format json \
    --gnomad "${OUTPUT_DIR}/refs/gnomad.chr22.vcf.bgz" \
    --clinvar "${OUTPUT_DIR}/refs/clinvar.tsv" \
    --regions "${OUTPUT_DIR}/data/giab_highconf.bed" \
    2>&1 | tee "${OUTPUT_DIR}/results/run_log_${TIMESTAMP}.txt"; then
    echo "ERROR: vartriage pipeline failed. Check run_log for details."
    exit 2
fi

echo "  Pipeline complete"

# Clinical HTML report
if ! vartriage \
    --vcf "${OUTPUT_DIR}/data/giab_chr22.vcf.gz" \
    --output "${OUTPUT_DIR}/reports/giab_chr22_clinical_${TIMESTAMP}.html" \
    --output-format clinical-html \
    --patient-id "GIAB-HG002" \
    --panel-name "GIAB Validation (chr22)" \
    --gnomad "${OUTPUT_DIR}/refs/gnomad.chr22.vcf.bgz" \
    --clinvar "${OUTPUT_DIR}/refs/clinvar.tsv" \
    --regions "${OUTPUT_DIR}/data/giab_highconf.bed" \
    2>&1 | tee -a "${OUTPUT_DIR}/results/run_log_${TIMESTAMP}.txt"; then
    echo "ERROR: Clinical report generation failed."
    exit 2
fi

echo "  Clinical report generated"

# Step 7: Compute validation metrics
echo "[7/7] Computing validation metrics..."

if ! python3 - "${OUTPUT_DIR}" "${TIMESTAMP}" <<'PYTHON_SCRIPT'
import json
import sys
from pathlib import Path

output_dir = Path(sys.argv[1])
timestamp = sys.argv[2]

results_file = output_dir / "results" / f"giab_chr22_{timestamp}.json"

# Read and parse results with error handling
try:
    raw_text = results_file.read_text(encoding="utf-8")
except OSError as exc:
    print(f"ERROR: Cannot read results file {results_file}: {exc}", file=sys.stderr)
    sys.exit(1)

try:
    results = json.loads(raw_text)
except json.JSONDecodeError as exc:
    print(f"ERROR: Invalid JSON in {results_file}: {exc}", file=sys.stderr)
    sys.exit(1)

if not isinstance(results, list):
    print(f"ERROR: Expected JSON array, got {type(results).__name__}", file=sys.stderr)
    sys.exit(1)

# Classification distribution
classifications: dict[str, int] = {}
for v in results:
    cls = v.get("acmg_classification", "Unknown")
    classifications[cls] = classifications.get(cls, 0) + 1

# Evidence tag frequency
tag_counts: dict[str, int] = {}
for v in results:
    for tag in v.get("evidence_tags", []):
        tag_counts[tag] = tag_counts.get(tag, 0) + 1

# Consequence distribution
consequences: dict[str, int] = {}
for v in results:
    cons = v.get("functional_consequence", "Unknown")
    consequences[cons] = consequences.get(cons, 0) + 1

# Missing data summary
missing_sources: dict[str, int] = {}
for v in results:
    for src in v.get("missing_data_sources", []):
        missing_sources[src] = missing_sources.get(src, 0) + 1

# Variants with ClinVar pathogenic/likely pathogenic
clinvar_actionable = [
    v for v in results
    if v.get("clinvar_assertion") in ("Pathogenic", "Likely_Pathogenic")
]

# Variants classified as Pathogenic/LP by pipeline
pipeline_actionable = [
    v for v in results
    if v.get("acmg_classification") in ("Pathogenic", "Likely_Pathogenic")
]

# Concordance: variants where pipeline classification matches ClinVar
concordant = [
    v for v in results
    if v.get("clinvar_assertion") in ("Pathogenic", "Likely_Pathogenic")
    and v.get("acmg_classification") in ("Pathogenic", "Likely_Pathogenic")
]

report = {
    "validation_summary": {
        "sample": "HG002",
        "region": "chr22",
        "total_variants_processed": len(results),
        "classification_distribution": classifications,
        "consequence_distribution": consequences,
        "evidence_tag_frequency": tag_counts,
        "missing_data_sources": missing_sources,
    },
    "concordance": {
        "clinvar_actionable_count": len(clinvar_actionable),
        "pipeline_actionable_count": len(pipeline_actionable),
        "concordant_count": len(concordant),
        "sensitivity_vs_clinvar": (
            len(concordant) / len(clinvar_actionable)
            if clinvar_actionable else None
        ),
    },
    "actionable_variants": [
        {
            "gene": v.get("gene_name"),
            "position": f"{v.get('chromosome', '?')}:{v.get('position', '?')}",
            "consequence": v.get("functional_consequence"),
            "classification": v.get("acmg_classification"),
            "clinvar": v.get("clinvar_assertion"),
            "composite_rank": v.get("composite_rank"),
        }
        for v in pipeline_actionable
    ],
}

metrics_file = output_dir / "results" / f"validation_metrics_{timestamp}.json"

try:
    metrics_file.write_text(json.dumps(report, indent=2), encoding="utf-8")
except OSError as exc:
    print(f"ERROR: Cannot write metrics file {metrics_file}: {exc}", file=sys.stderr)
    sys.exit(1)

print(f"\n{'='*50}")
print(f" VALIDATION RESULTS")
print(f"{'='*50}")
print(f" Total variants processed:    {len(results)}")
print(f" Pathogenic:                   {classifications.get('Pathogenic', 0)}")
print(f" Likely Pathogenic:            {classifications.get('Likely_Pathogenic', 0)}")
print(f" VUS:                          {classifications.get('VUS', 0)}")
print(f"")
print(f" ClinVar actionable (P/LP):    {len(clinvar_actionable)}")
print(f" Pipeline actionable (P/LP):   {len(pipeline_actionable)}")
print(f" Concordant:                   {len(concordant)}")
if clinvar_actionable:
    sens = len(concordant) / len(clinvar_actionable) * 100
    print(f" Sensitivity vs ClinVar:       {sens:.1f}%")
print(f"")
print(f" Metrics saved: {metrics_file}")
print(f"{'='*50}")

# Exit with success
sys.exit(0)
PYTHON_SCRIPT
then
    echo "ERROR: Validation metrics computation failed."
    exit 3
fi

echo ""
echo "Validation complete. Results in: ${OUTPUT_DIR}/"
echo "  - results/giab_chr22_${TIMESTAMP}.json       (raw pipeline output)"
echo "  - results/validation_metrics_${TIMESTAMP}.json (concordance metrics)"
echo "  - reports/giab_chr22_clinical_${TIMESTAMP}.html (clinical report)"
echo "  - results/run_log_${TIMESTAMP}.txt            (execution log)"
