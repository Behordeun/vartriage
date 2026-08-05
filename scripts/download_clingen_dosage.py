"""Download and process ClinGen dosage sensitivity data.

Fetches the latest ClinGen gene-level dosage sensitivity curation
from the ClinGen FTP site, extracts haploinsufficiency (HI) and
triplosensitivity (TS) scores, and writes a clean TSV suitable
for use with the SV triage pipeline.

Usage:
    python scripts/download_clingen_dosage.py [--output path/to/output.tsv]

Output TSV columns:
    gene_symbol    HI score (0-3)    TS score (0-3)

Score scale:
    0 = No evidence
    1 = Little evidence
    2 = Emerging evidence
    3 = Sufficient evidence
    40 = Dosage sensitivity unlikely
"""

from __future__ import annotations

import argparse
import csv
import sys
import urllib.request
from pathlib import Path

CLINGEN_DOSAGE_URL = (
    "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
)

DEFAULT_OUTPUT = Path("vartriage/data/clingen_dosage.tsv")


def _find_header_idx(lines: list[str]) -> int:
    """Return the index of the header line in the raw TSV lines."""
    for i, line in enumerate(lines):
        if "Gene Symbol" in line or "gene_symbol" in line.lower():
            return i
    return 0


def _build_col_map(header: list[str]) -> dict[str, int]:
    """Map logical column names to indices from the header row."""
    col_map: dict[str, int] = {}
    for idx, col_name in enumerate(header):
        col_lower = col_name.strip().lower()
        if "gene symbol" in col_lower or col_lower == "gene_symbol":
            col_map["gene_symbol"] = idx
        elif "haploinsufficiency" in col_lower and "score" in col_lower:
            col_map["hi_score"] = idx
        elif "triplosensitivity" in col_lower and "score" in col_lower:
            col_map["ts_score"] = idx
    return col_map


def _write_dosage_tsv(
    output_path: Path,
    lines: list[str],
    header_idx: int,
    col_map: dict[str, int],
) -> int:
    """Write processed dosage rows to output_path, return gene count."""
    gene_idx = col_map["gene_symbol"]
    gene_count = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_symbol", "hi_score", "ts_score"])
        for line in lines[header_idx + 1 :]:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if gene_idx >= len(fields):
                continue
            symbol = fields[gene_idx].strip()
            if not symbol:
                continue
            hi_score = _extract_score(fields, col_map.get("hi_score"))
            ts_score = _extract_score(fields, col_map.get("ts_score"))
            if hi_score == "" and ts_score == "":
                continue
            writer.writerow([symbol, hi_score, ts_score])
            gene_count += 1
    return gene_count


def fetch_dosage_data(output_path: Path) -> None:
    """Download and process ClinGen dosage sensitivity data."""
    # Validate output path does not escape expected directory
    output_path.resolve()
    if ".." in str(output_path):
        print("Error: output path must not contain '..'", file=sys.stderr)
        sys.exit(2)

    print(f"Downloading ClinGen dosage data from:\n  {CLINGEN_DOSAGE_URL}")

    try:
        with urllib.request.urlopen(CLINGEN_DOSAGE_URL, timeout=30) as response:
            raw_data = response.read().decode("utf-8")
    except Exception as exc:
        print(f"Error downloading: {exc}", file=sys.stderr)
        sys.exit(1)

    lines = raw_data.strip().split("\n")
    header_idx = _find_header_idx(lines)
    header = lines[header_idx].lstrip("#").strip().split("\t")
    col_map = _build_col_map(header)

    if "gene_symbol" not in col_map:
        print(
            "Error: could not find gene symbol column in ClinGen data", file=sys.stderr
        )
        sys.exit(1)

    gene_count = _write_dosage_tsv(output_path, lines, header_idx, col_map)
    print(f"Wrote {gene_count} genes to {output_path}")


def _extract_score(fields: list[str], col_idx: int | None) -> str:
    """Extract a numeric score from a field, handling edge cases."""
    if col_idx is None or col_idx >= len(fields):
        return ""
    raw = fields[col_idx].strip()
    if not raw or raw in ("N/A", "Not yet evaluated", "-"):
        return ""
    # ClinGen uses 0, 1, 2, 3, or 40 (dosage sensitivity unlikely)
    try:
        val = int(raw)
        # Treat 40 as 0 (unlikely = no sensitivity)
        if val == 40:
            return "0"
        return str(val)
    except ValueError:
        return ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download ClinGen dosage sensitivity data"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()
    fetch_dosage_data(args.output)


if __name__ == "__main__":
    main()
