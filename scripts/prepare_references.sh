#!/usr/bin/env bash
# prepare_references.sh - Normalize reference files for vartriage
#
# Ensures all reference files use consistent chromosome naming (chr-prefixed)
# and validates file formats before pipeline execution.
#
# Run from the project root after downloading reference files:
#   ./scripts/prepare_references.sh [--refs-dir data/references]
#
# What this script does:
#   1. Detects chromosome naming convention in each TSV file
#   2. Adds 'chr' prefix to any file using bare chromosome names (1, 2, ..., X, Y)
#   3. Validates required column headers exist
#   4. Reports any format issues that would cause pipeline failures

set -euo pipefail

REFS_DIR="${1:-data/references}"

echo "============================================="
echo " vartriage Reference File Preparation"
echo " Directory: ${REFS_DIR}"
echo "============================================="
echo ""

if [ ! -d "${REFS_DIR}" ]; then
    echo "ERROR: Reference directory not found: ${REFS_DIR}"
    echo "Run scripts/download_test_data.sh first."
    exit 1
fi

# Track issues
ISSUES=0

# --- Helper functions ---

detect_chr_style() {
    # Returns "chr" if file uses chr-prefixed names, "bare" if not, "empty" if no data
    local file="$1"
    local first_chrom
    first_chrom=$(awk -F'\t' 'NR>1 && $0 !~ /^#/ {print $1; exit}' "$file" 2>/dev/null)
    
    if [ -z "$first_chrom" ]; then
        echo "empty"
    elif [[ "$first_chrom" == chr* ]]; then
        echo "chr"
    else
        echo "bare"
    fi
}

add_chr_prefix() {
    local file="$1"
    local tmp="${file}.tmp"
    
    echo "    Fixing: adding 'chr' prefix to chromosome names..."
    awk -F'\t' 'BEGIN{OFS="\t"} 
        NR==1 || /^#/ {print; next}
        {$1="chr"$1; print}
    ' "$file" > "$tmp" && mv "$tmp" "$file"
}

validate_header() {
    local file="$1"
    shift
    local expected_cols=("$@")
    
    # Get first non-comment line as header
    local header
    header=$(awk '/^[^#]/{print; exit}' "$file")
    
    for col in "${expected_cols[@]}"; do
        if ! echo "$header" | grep -qi "$col"; then
            echo "    WARNING: Missing expected column '$col' in header"
            echo "    Header found: $header"
            ((ISSUES++))
            return 1
        fi
    done
    return 0
}

# --- Process ClinVar ---

echo "[1/4] ClinVar TSV"
CLINVAR_FILE="${REFS_DIR}/clinvar.tsv"
if [ -f "$CLINVAR_FILE" ]; then
    style=$(detect_chr_style "$CLINVAR_FILE")
    if [ "$style" = "bare" ]; then
        add_chr_prefix "$CLINVAR_FILE"
        echo "    Fixed: chromosome names now chr-prefixed"
    elif [ "$style" = "chr" ]; then
        echo "    OK: already chr-prefixed"
    else
        echo "    WARNING: file appears empty or unreadable"
        ((ISSUES++))
    fi
    validate_header "$CLINVAR_FILE" "chrom" "pos" "ref" "alt" "clinical_significance" || true
else
    echo "    SKIP: file not found"
fi

echo ""

# --- Process gnomAD ---

echo "[2/4] gnomAD TSV"
GNOMAD_FILE="${REFS_DIR}/gnomad_chr22.tsv"
if [ -f "$GNOMAD_FILE" ]; then
    style=$(detect_chr_style "$GNOMAD_FILE")
    if [ "$style" = "bare" ]; then
        add_chr_prefix "$GNOMAD_FILE"
        echo "    Fixed: chromosome names now chr-prefixed"
    elif [ "$style" = "chr" ]; then
        echo "    OK: already chr-prefixed"
    else
        echo "    WARNING: file appears empty or unreadable"
        ((ISSUES++))
    fi
    validate_header "$GNOMAD_FILE" "chrom" "pos" "ref" "alt" "af" || true
else
    echo "    SKIP: file not found"
fi

echo ""

# --- Process CADD ---

echo "[3/4] CADD scores TSV"
# Check both possible CADD file names
CADD_FILE=""
for candidate in "${REFS_DIR}/cadd_chr22_giab.tsv" "${REFS_DIR}/cadd_chr22_sample.tsv"; do
    if [ -f "$candidate" ]; then
        CADD_FILE="$candidate"
        break
    fi
done

if [ -n "$CADD_FILE" ]; then
    # CADD files use # comment lines for header
    style=$(detect_chr_style "$CADD_FILE")
    if [ "$style" = "bare" ]; then
        add_chr_prefix "$CADD_FILE"
        echo "    Fixed: chromosome names now chr-prefixed"
    elif [ "$style" = "chr" ]; then
        echo "    OK: already chr-prefixed (${CADD_FILE})"
    else
        echo "    WARNING: file appears empty or unreadable"
        ((ISSUES++))
    fi
    
    # CADD uses #chrom header style
    line_count=$(grep -cv "^#" "$CADD_FILE" 2>/dev/null || echo "0")
    echo "    Entries: ${line_count}"
    if [ "$line_count" -lt 100 ]; then
        echo "    NOTE: CADD file has few entries. For meaningful results,"
        echo "          download the full chr22 CADD scores from:"
        echo "          https://cadd.gs.washington.edu/download"
    fi
else
    echo "    SKIP: no CADD file found"
fi

echo ""

# --- Process GENCODE GTF ---

echo "[4/4] GENCODE GTF"
GTF_FILE=""
for candidate in "${REFS_DIR}/gencode_chr22.gtf" "${REFS_DIR}/gencode.v46.annotation.gtf"; do
    if [ -f "$candidate" ]; then
        GTF_FILE="$candidate"
        break
    fi
done

if [ -n "$GTF_FILE" ]; then
    # GTF uses chr prefix by default from GENCODE
    first_chr=$(awk '$0 !~ /^#/ {print $1; exit}' "$GTF_FILE")
    if [[ "$first_chr" == chr* ]]; then
        echo "    OK: chr-prefixed (${GTF_FILE})"
    else
        echo "    WARNING: GTF does not use chr prefix. GENCODE GRCh38 should use chr-prefixed names."
        echo "    If using Ensembl GTF (bare names), run: sed -i '' 's/^\\([0-9XY]\\)/chr\\1/' ${GTF_FILE}"
        ((ISSUES++))
    fi
else
    echo "    SKIP: no GTF file found"
fi

echo ""

# --- Summary ---

echo "============================================="
if [ "$ISSUES" -gt 0 ]; then
    echo " DONE with ${ISSUES} warning(s). Review above."
else
    echo " DONE. All reference files validated."
fi
echo "============================================="
echo ""
echo "You can now run the pipeline:"
echo "  vartriage --vcf your_file.vcf.gz \\"
echo "    --output results.json \\"
echo "    --gene-annotation ${REFS_DIR}/gencode_chr22.gtf \\"
echo "    --gnomad ${REFS_DIR}/gnomad_chr22.tsv \\"
echo "    --clinvar ${REFS_DIR}/clinvar.tsv \\"
echo "    --cadd-scores ${REFS_DIR}/cadd_chr22_giab.tsv"

exit $ISSUES
