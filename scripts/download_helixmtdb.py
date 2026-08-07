"""Download and process HelixMTdb mitochondrial frequency data.

Fetches mtDNA variant frequency data from the HelixMTdb public
dataset, extracts position/ref/alt/AF/allele_count, and writes a
standardized TSV for use by the mitochondrial pipeline.

Usage:
    python scripts/download_helixmtdb.py [--output path/to/output.tsv]

Output TSV columns:
    position    ref    alt    af    allele_count

Data source:
    HelixMTdb is maintained by Helix and provides allele frequencies
    for mtDNA variants observed in their WGS cohort (~200k samples).
    The VCF is available from:
    https://www.helix.com/pages/mitochondrial-variant-database
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import sys
import tempfile
import urllib.request
from pathlib import Path

# HelixMTdb VCF download URL (sites-only VCF)
HELIXMTDB_URL = (
    "https://helix-research-public.s3.us-east-1.amazonaws.com/"
    "mito/HelixMTdb_20200327.vcf.gz"
)

DEFAULT_OUTPUT = Path("vartriage/data/mito/helixmtdb_frequency.tsv")


def _download_file(url: str, dest: Path) -> None:
    """Download a file from a URL to a local path."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "vartriage-data-updater/0.15.0"},
    )
    with (
        urllib.request.urlopen(request, timeout=120) as response,
        open(dest, "wb") as f,
    ):
        while True:
            chunk = response.read(65536)
            if not chunk:
                break
            f.write(chunk)


def _parse_vcf_line(line: str) -> dict[str, str] | None:
    """Parse a single VCF data line into a frequency record.

    Extracts AF and AC from the INFO field. Returns None for
    header lines or unparseable records.
    """
    if line.startswith("#"):
        return None

    fields = line.strip().split("\t")
    if len(fields) < 8:
        return None

    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    info = fields[7]

    # Skip multi-allelic (should be rare for mtDNA)
    if "," in alt:
        return None

    # Parse AF from INFO field
    af = None
    ac = None
    for kv in info.split(";"):
        if kv.startswith("AF="):
            af = kv.split("=", 1)[1]
        elif kv.startswith("AC="):
            ac = kv.split("=", 1)[1]

    if af is None:
        return None

    try:
        af_float = float(af)
    except ValueError:
        return None

    allele_count = "0"
    if ac is not None:
        with contextlib.suppress(ValueError):
            allele_count = str(int(ac))

    return {
        "position": pos,
        "ref": ref.upper(),
        "alt": alt.upper(),
        "af": f"{af_float:.6g}",
        "allele_count": allele_count,
    }


def _parse_vcf(vcf_path: Path) -> list[dict[str, str]]:
    """Parse a gzipped VCF file into frequency records."""
    entries: list[dict[str, str]] = []

    open_func = gzip.open if str(vcf_path).endswith(".gz") else open

    with open_func(vcf_path, "rt", encoding="utf-8") as f:  # type: ignore[call-overload]
        for line in f:
            record = _parse_vcf_line(line)
            if record is not None:
                entries.append(record)

    return entries


def _write_tsv(entries: list[dict[str, str]], output_path: Path) -> None:
    """Write parsed entries to the output TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["position", "ref", "alt", "af", "allele_count"],
            delimiter="\t",
        )
        writer.writeheader()
        # Sort by position
        sorted_entries = sorted(entries, key=lambda e: int(e["position"]))
        writer.writerows(sorted_entries)


def main(argv: list[str] | None = None) -> None:
    """Download HelixMTdb VCF and produce a frequency TSV."""
    parser = argparse.ArgumentParser(
        description="Download HelixMTdb and produce a frequency TSV"
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
        default=HELIXMTDB_URL,
        help="HelixMTdb VCF URL (override for testing)",
    )
    parser.add_argument(
        "--input-vcf",
        type=Path,
        default=None,
        help="Use a local VCF file instead of downloading",
    )
    args = parser.parse_args(argv)

    if args.input_vcf is not None:
        vcf_path = args.input_vcf
        print(f"Using local VCF: {vcf_path}")
    else:
        print(f"Downloading HelixMTdb VCF from: {args.url}")
        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tmp:
            vcf_path = Path(tmp.name)
        try:
            _download_file(args.url, vcf_path)
        except Exception as exc:
            print(
                f"Error: failed to download HelixMTdb: {exc}",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Downloaded to: {vcf_path}")

    entries = _parse_vcf(vcf_path)
    if not entries:
        print(
            "Warning: no entries parsed from VCF. File may be empty or format changed.",
            file=sys.stderr,
        )
        sys.exit(1)

    _write_tsv(entries, args.output)
    print(f"Wrote {len(entries)} HelixMTdb entries to: {args.output}")


if __name__ == "__main__":
    main()
