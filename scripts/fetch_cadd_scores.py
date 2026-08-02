"""Fetch CADD scores for GIAB chr22 SNVs via the CADD REST API.

Queries the CADD v1.7 GRCh38 API one position at a time for each SNV
in the GIAB VCF. Writes a vartriage-compatible TSV.

For large variant sets (>5000), this takes a while due to rate limits.
The script saves progress incrementally so you can resume if interrupted.

Usage:
    python scripts/fetch_cadd_scores.py \
        --vcf validation_results/data/giab_chr22.vcf.gz \
        --output data/references/cadd_chr22_full.tsv

Requirements:
    - pysam (already installed with vartriage)
    - requests: pip install requests
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    print("ERROR: requests required. Install with: pip install requests")
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("ERROR: pysam required. Install with: pip install pysam")
    sys.exit(1)


CADD_API_BASE = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7"
REQUESTS_PER_SECOND = 2
BATCH_SAVE_INTERVAL = 500


def extract_snvs(vcf_path: Path) -> list[tuple[str, int, str, str]]:
    """Extract single-nucleotide variants from a VCF file."""
    snvs = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            ref = record.ref
            for alt in record.alts or []:
                if len(ref) == 1 and len(alt) == 1 and ref != alt:
                    snvs.append((record.chrom, record.pos, ref, alt))
    print(f"Extracted {len(snvs)} SNVs from {vcf_path}")
    return snvs


def load_progress(output_path: Path) -> dict[str, float]:
    """Load previously fetched scores to enable resume."""
    scores = {}
    if output_path.exists():
        with open(output_path) as f:
            f.readline()  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 5:
                    key = f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"
                    scores[key] = float(parts[4])
    if scores:
        print(f"Resuming: {len(scores)} scores already fetched")
    return scores


def fetch_single_score(chrom: str, pos: int, ref: str, alt: str) -> float | None:
    """Query CADD API for a single variant. Returns Phred score or None."""
    chrom_clean = chrom.replace("chr", "")
    url = f"{CADD_API_BASE}/{chrom_clean}:{pos}_{ref}_{alt}"

    max_retries = 5
    backoff = 5.0

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0:
                    return float(data[0].get("PHRED", 0))
                return None
            elif response.status_code == 429:
                time.sleep(backoff * (2 ** attempt))
                continue
            else:
                return None
        except (requests.RequestException, ValueError, KeyError):
            return None

    return None


def write_scores(
    scores: dict[str, float], output_path: Path
) -> None:
    """Write all scores to the output TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("chrom\tpos\tref\talt\tscore\n")
        for key, score in sorted(scores.items()):
            parts = key.split(":")
            f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{score:.4f}\n")


def _print_progress(
    done: int, total: int, fetched: int, failed: int,
    start_time: float,
) -> None:
    """Print a progress line every 100 variants."""
    if done % 100 != 0:
        return
    elapsed = time.time() - start_time
    rate = fetched / elapsed if elapsed > 0 else 0
    eta = (total - done) / (REQUESTS_PER_SECOND * 60) if REQUESTS_PER_SECOND > 0 else 0
    print(
        f"  [{done}/{total}] fetched={fetched} failed={failed} "
        f"rate={rate:.1f}/s ETA={eta:.0f}min"
    )


def _fetch_loop(
    snvs: list,
    scores: dict[str, float],
    output_path: Path,
    start_time: float,
) -> tuple[int, int, int]:
    """Fetch CADD scores for all SNVs; return (fetched, skipped, failed)."""
    total = len(snvs)
    fetched = skipped = failed = 0
    for chrom, pos, ref, alt in snvs:
        key = f"{chrom}:{pos}:{ref}:{alt}"
        if key in scores:
            skipped += 1
        else:
            score = fetch_single_score(chrom, pos, ref, alt)
            if score is not None:
                scores[key] = score
                fetched += 1
            else:
                failed += 1
        time.sleep(1.0 / REQUESTS_PER_SECOND)
        done = skipped + fetched + failed
        _print_progress(done, total, fetched, failed, start_time)
        if fetched > 0 and fetched % BATCH_SAVE_INTERVAL == 0:
            write_scores(scores, output_path)
    return fetched, skipped, failed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch CADD scores for GIAB SNVs via the CADD REST API"
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        default=Path("validation_results/data/giab_chr22.vcf.gz"),
        help="Input VCF file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/references/cadd_chr22_full.tsv"),
        help="Output TSV path",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Max variants to fetch (0 = all). Use for testing.",
    )
    args = parser.parse_args()

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}")
        sys.exit(1)

    snvs = extract_snvs(args.vcf)
    if args.limit > 0:
        snvs = snvs[:args.limit]

    # Load existing progress
    scores = load_progress(args.output)

    total = len(snvs)
    start_time = time.time()

    print(f"Fetching CADD scores for {total} SNVs (rate: {REQUESTS_PER_SECOND}/sec)")
    print(f"Estimated time: {total / REQUESTS_PER_SECOND / 60:.0f} minutes")
    print()

    _, _, failed = _fetch_loop(snvs, scores, args.output, start_time)

    # Final save
    write_scores(scores, args.output)

    elapsed = time.time() - start_time
    print(f"\nDone in {elapsed:.0f}s. Scored: {len(scores)}/{total} "
          f"({100 * len(scores) / total:.1f}%)")
    print(f"Failed: {failed}. Output: {args.output}")


if __name__ == "__main__":
    main()
