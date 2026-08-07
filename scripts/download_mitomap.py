"""Download and process MITOMAP pathogenic mutations table.

Fetches the MITOMAP disease-associated mtDNA mutations from the
publicly available resources, parses confirmed and reported pathogenic
variants, and writes a standardized TSV for use by the mitochondrial
pipeline.

Usage:
    python scripts/download_mitomap.py [--output path/to/output.tsv]

Output TSV columns:
    position    ref    alt    disease    status    locus

Status values:
    Cfrm = Confirmed pathogenic
    Reported = Reported association
    P.M. = Point mutation (functional evidence)
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
from pathlib import Path

MITOMAP_URL = "https://www.mitomap.org/foswiki/bin/view/MITOMAP/ConfirmedMutations"

DEFAULT_OUTPUT = Path("vartriage/data/mito/mitomap_pathogenic.tsv")


def _parse_position(pos_str: str) -> int | None:
    """Extract numeric position from MITOMAP position strings."""
    match = re.search(r"(\d+)", pos_str.strip())
    if match:
        return int(match.group(1))
    return None


def _parse_nucleotide_change(change_str: str) -> tuple[str, str] | None:
    """Parse a nucleotide change like 'A-G' or 'T>C' into (ref, alt).

    MITOMAP uses various separators: '-', '>', '/', 'to'.
    """
    cleaned = change_str.strip().upper()
    for sep in ["-", ">", "/", " TO "]:
        if sep in cleaned:
            parts = cleaned.split(sep, 1)
            ref = parts[0].strip()
            alt = parts[1].strip()
            if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                return ref, alt
    return None


def _download_mitomap_html(url: str) -> str:
    """Fetch the MITOMAP mutations page as raw HTML."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "vartriage-data-updater/0.15.0"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def _parse_mitomap_table(html: str) -> list[dict[str, str]]:
    """Parse MITOMAP HTML table rows into structured records.

    Extracts position, nucleotide change, disease, status, and locus
    from the confirmed mutations table. Falls back to regex-based
    extraction since MITOMAP doesn't provide a stable machine-readable
    download format.
    """
    entries: list[dict[str, str]] = []

    # Match table rows with mtDNA mutation data
    # Pattern targets the typical MITOMAP table structure
    row_pattern = re.compile(r"<tr[^>]*>.*?</tr>", re.DOTALL | re.IGNORECASE)
    cell_pattern = re.compile(r"<td[^>]*>(.*?)</td>", re.DOTALL | re.IGNORECASE)

    for row_match in row_pattern.finditer(html):
        row_html = row_match.group(0)
        cells = cell_pattern.findall(row_html)

        if len(cells) < 5:
            continue

        # Clean HTML tags from cell content
        clean_cells = [re.sub(r"<[^>]+>", "", c).strip() for c in cells]

        # Try to parse as a mutation row
        pos = _parse_position(clean_cells[0])
        if pos is None or pos < 1 or pos > 16569:
            continue

        nucleotide_change = _parse_nucleotide_change(clean_cells[1])
        if nucleotide_change is None:
            continue

        ref, alt = nucleotide_change
        disease = clean_cells[2] if len(clean_cells) > 2 else ""
        status = clean_cells[3] if len(clean_cells) > 3 else "Reported"
        locus = clean_cells[4] if len(clean_cells) > 4 else ""

        # Normalize status
        if "cfrm" in status.lower() or "confirmed" in status.lower():
            status = "Cfrm"
        elif "reported" in status.lower():
            status = "Reported"
        elif "p.m." in status.lower() or "point" in status.lower():
            status = "P.M."

        entries.append(
            {
                "position": str(pos),
                "ref": ref,
                "alt": alt,
                "disease": disease,
                "status": status,
                "locus": locus,
            }
        )

    return entries


def _write_tsv(entries: list[dict[str, str]], output_path: Path) -> None:
    """Write parsed entries to the output TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["position", "ref", "alt", "disease", "status", "locus"],
            delimiter="\t",
        )
        writer.writeheader()
        # Sort by position for consistent output
        sorted_entries = sorted(entries, key=lambda e: int(e["position"]))
        writer.writerows(sorted_entries)


def main(argv: list[str] | None = None) -> None:
    """Download MITOMAP data and write processed TSV."""
    parser = argparse.ArgumentParser(
        description="Download MITOMAP pathogenic mutations and produce a TSV"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--url",
        type=str,
        default=MITOMAP_URL,
        help="MITOMAP source URL (override for testing)",
    )
    args = parser.parse_args(argv)

    print(f"Downloading MITOMAP confirmed mutations from: {args.url}")
    try:
        html = _download_mitomap_html(args.url)
    except Exception as exc:
        print(f"Error: failed to download MITOMAP data: {exc}", file=sys.stderr)
        sys.exit(1)

    entries = _parse_mitomap_table(html)
    if not entries:
        print(
            "Warning: no entries parsed from MITOMAP HTML. "
            "The page format may have changed. No output written.",
            file=sys.stderr,
        )
        sys.exit(1)

    _write_tsv(entries, args.output)
    print(f"Wrote {len(entries)} MITOMAP entries to: {args.output}")


if __name__ == "__main__":
    main()
