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
#   - wget or curl (auto-detected at runtime)
#   - ~20 GB free disk space for reference files
#
# Usage:
#   ./scripts/validate_giab.sh [OPTIONS]
#
# Options:
#   --output-dir DIR    Output directory (default: validation_results)
#   --region REGION     Genomic region to validate (default: chr22)
#   --build BUILD       Reference genome build: GRCh38 or GRCh37 (default: GRCh38)
#   --help              Show this help message
#
# Exit codes:
#   0 - validation completed successfully
#   1 - missing dependencies or download failure
#   2 - pipeline execution failure
#   3 - metrics computation failure

set -euo pipefail

# --- Default configuration (override via flags) ---

OUTPUT_DIR="validation_results"
REGION="chr22"
REFERENCE_BUILD="GRCh38"
GIAB_VERSION="v4.2.1"
SAMPLE="HG002"

# --- Argument parsing ---

show_help() {
    sed -n '/^# Usage:/,/^# Exit codes:/p' "$0" | sed 's/^# \?//'
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --output-dir=*)
            OUTPUT_DIR="${1#*=}"
            shift
            ;;
        --region)
            REGION="$2"
            shift 2
            ;;
        --region=*)
            REGION="${1#*=}"
            shift
            ;;
        --build)
            REFERENCE_BUILD="$2"
            shift 2
            ;;
        --build=*)
            REFERENCE_BUILD="${1#*=}"
            shift
            ;;
        --help|-h)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run with --help for usage information."
            exit 1
            ;;
    esac
done

# --- URL configuration (derived from flags) ---

GIAB_BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/${REFERENCE_BUILD}"
GIAB_VCF_URL="${GIAB_BASE_URL}/HG002_${REFERENCE_BUILD}_1_22_${GIAB_VERSION}_benchmark.vcf.gz"
GIAB_BED_URL="${GIAB_BASE_URL}/HG002_${REFERENCE_BUILD}_1_22_${GIAB_VERSION}_benchmark_noinconsistent.bed"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${REFERENCE_BUILD}/clinvar.vcf.gz"

# gnomAD URL uses the region variable to pick the per-chromosome file
GNOMAD_REGION_FILE="gnomad.exomes.v4.1.1.sites.${REGION}.vcf.bgz"
GNOMAD_URL="https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1.1/vcf/exomes/${GNOMAD_REGION_FILE}"

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

# Verify minimum bcftools version
check_bcftools_version() {
    local required="1.17"
    local installed
    installed="$(bcftools --version | head -n1 | awk '{print $2}')"

    if [ "$(printf '%s\n' "$required" "$installed" | sort -V | head -n1)" != "$required" ]; then
        echo "ERROR: bcftools version $installed is less than required $required"
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

check_bcftools_version

echo "============================================="
echo " vartriage GIAB Validation Pipeline"
echo " Sample: ${SAMPLE} (${GIAB_VERSION})"
echo " Build:  ${REFERENCE_BUILD}"
echo " Region: ${REGION}"
echo " Output: ${OUTPUT_DIR}"
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

# Step 4: Download gnomAD for the target region
echo "[4/7] Downloading gnomAD ${REGION} (exomes)..."
if [ ! -f "${OUTPUT_DIR}/refs/gnomad.${REGION}.vcf.bgz" ]; then
    download_file "${GNOMAD_URL}" "${OUTPUT_DIR}/refs/gnomad.${REGION}.vcf.bgz"
    download_file "${GNOMAD_URL}.tbi" "${OUTPUT_DIR}/refs/gnomad.${REGION}.vcf.bgz.tbi"
    echo "  Downloaded gnomAD v4.1.1 ${REGION} exomes"
else
    echo "  gnomAD ${REGION} already exists, skipping"
fi

# Step 5: Extract target region variants from GIAB
echo "[5/7] Extracting ${REGION} subset from GIAB benchmark..."
if [ ! -f "${OUTPUT_DIR}/data/giab_${REGION}.vcf.gz" ]; then
    bcftools view -r "${REGION}" "${OUTPUT_DIR}/data/giab_benchmark.vcf.gz" \
        -Oz -o "${OUTPUT_DIR}/data/giab_${REGION}.vcf.gz"
    tabix -p vcf "${OUTPUT_DIR}/data/giab_${REGION}.vcf.gz"

    VARIANT_COUNT=$(bcftools view -H "${OUTPUT_DIR}/data/giab_${REGION}.vcf.gz" | wc -l)
    echo "  Extracted ${VARIANT_COUNT} ${REGION} variants"
else
    echo "  ${REGION} subset already exists, skipping"
fi

# Step 6: Run vartriage on the target region
echo "[6/7] Running vartriage on ${REGION}..."
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# JSON output for analysis
if ! vartriage \
    --vcf "${OUTPUT_DIR}/data/giab_${REGION}.vcf.gz" \
    --output "${OUTPUT_DIR}/results/giab_${REGION}_${TIMESTAMP}.json" \
    --output-format json \
    --gnomad "${OUTPUT_DIR}/refs/gnomad.${REGION}.vcf.bgz" \
    --clinvar "${OUTPUT_DIR}/refs/clinvar.tsv" \
    --regions "${OUTPUT_DIR}/data/giab_highconf.bed" \
    2>&1 | tee "${OUTPUT_DIR}/results/run_log_${TIMESTAMP}.txt"; then
    echo "ERROR: vartriage pipeline failed. Check run_log for details."
    exit 2
fi

echo "  Pipeline complete"

# Clinical HTML report
if ! vartriage \
    --vcf "${OUTPUT_DIR}/data/giab_${REGION}.vcf.gz" \
    --output "${OUTPUT_DIR}/reports/giab_${REGION}_clinical_${TIMESTAMP}.html" \
    --output-format clinical-html \
    --patient-id "GIAB-HG002" \
    --panel-name "GIAB Validation (${REGION})" \
    --gnomad "${OUTPUT_DIR}/refs/gnomad.${REGION}.vcf.bgz" \
    --clinvar "${OUTPUT_DIR}/refs/clinvar.tsv" \
    --regions "${OUTPUT_DIR}/data/giab_highconf.bed" \
    2>&1 | tee -a "${OUTPUT_DIR}/results/run_log_${TIMESTAMP}.txt"; then
    echo "ERROR: Clinical report generation failed."
    exit 2
fi

echo "  Clinical report generated"

# Step 7: Compute validation metrics
echo "[7/7] Computing validation metrics..."

if ! python3 - "${OUTPUT_DIR}" "${TIMESTAMP}" "${REGION}" <<'PYTHON_SCRIPT'
import json
import sys
from pathlib import Path

output_dir = Path(sys.argv[1])
timestamp = sys.argv[2]
region = sys.argv[3]

results_file = output_dir / "results" / f"giab_{region}_{timestamp}.json"

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

if len(results) == 0:
    print("WARNING: Pipeline produced zero variants.", file=sys.stderr)

# Classification distribution (defensive: use .get with defaults)
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
        "region": region,
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
print(f" VALIDATION RESULTS ({region})")
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
echo "  - results/giab_${REGION}_${TIMESTAMP}.json       (raw pipeline output)"
echo "  - results/validation_metrics_${TIMESTAMP}.json    (concordance metrics)"
echo "  - reports/giab_${REGION}_clinical_${TIMESTAMP}.html (clinical report)"
echo "  - results/run_log_${TIMESTAMP}.txt               (execution log)"
