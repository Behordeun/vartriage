"""Build a ClinVar expert-curated truth set for validation benchmarking.

Implements the truth-set construction stage of the validation-benchmarks
spec: filter a ClinVar GRCh38 VCF to the two highest-confidence review
tiers, exclude conflicting interpretations, split into a pathogenic
(P/LP) positive set and a benign (B/LB) negative set, normalize each
variant through the pipeline's own normalizer, and emit a truth VCF plus
a labels TSV.

The ClinVar release is pinned by filename (clinvar_YYYYMMDD.vcf.gz) and
recorded in the provenance sidecar for reproducibility.

Usage:
    python scripts/build_truth_set.py \
        --clinvar /path/to/clinvar_20260928.vcf.gz \
        --out-dir validation_results/truth_set
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from vartriage._internal.normalizer import VariantNormalizer

# The two review tiers that constitute an expert-curated truth set.
_TRUTH_REVIEW_STATUSES: frozenset[str] = frozenset(
    {
        "reviewed_by_expert_panel",
        "practice_guideline",
    }
)

# CLNSIG tokens that place a variant in the positive (pathogenic) set.
_POSITIVE_SIG: frozenset[str] = frozenset(
    {
        "Pathogenic",
        "Likely_pathogenic",
        "Pathogenic/Likely_pathogenic",
    }
)

# CLNSIG tokens that place a variant in the negative (benign) set.
_NEGATIVE_SIG: frozenset[str] = frozenset(
    {
        "Benign",
        "Likely_benign",
        "Benign/Likely_benign",
    }
)


@dataclass
class _Counters:
    """Running tally of filtering decisions for the provenance record."""

    total: int = 0
    wrong_review_status: int = 0
    conflicting: int = 0
    unmapped_significance: int = 0
    positive: int = 0
    negative: int = 0
    normalization_failed: int = 0
    dropped_reasons: dict[str, int] = field(default_factory=dict)


def _parse_info(info: str) -> dict[str, str]:
    """Parse a VCF INFO column into a flat key->value dict (flags map to '')."""
    fields: dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            key, _, value = item.partition("=")
            fields[key] = value
        else:
            fields[item] = ""
    return fields


def _review_status_ok(info: dict[str, str]) -> bool:
    """True when CLNREVSTAT is one of the truth-set review tiers."""
    return info.get("CLNREVSTAT", "") in _TRUTH_REVIEW_STATUSES


def _classify_significance(info: dict[str, str]) -> str | None:
    """Map CLNSIG to 'positive', 'negative', or None (not a truth label).

    A record carrying a conflicting classification is rejected here, since
    the review-status filter already excludes the conflicting tier but a
    CLNSIGCONF field is a second guard.
    """
    if info.get("CLNSIGCONF"):
        return None
    sig = info.get("CLNSIG", "")
    if sig in _POSITIVE_SIG:
        return "positive"
    if sig in _NEGATIVE_SIG:
        return "negative"
    return None


def _gene_symbol(info: dict[str, str]) -> str:
    """Extract the first gene symbol from GENEINFO (symbol:id|symbol:id)."""
    geneinfo = info.get("GENEINFO", "")
    if not geneinfo:
        return ""
    return geneinfo.split("|")[0].split(":")[0]


def _open_maybe_gzip(path: Path):
    """Open a plain or gzipped text file transparently."""
    if path.suffix == ".gz" or path.suffix == ".bgz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def build_truth_set(clinvar_path: Path, out_dir: Path, fasta_path: Path) -> _Counters:
    """Construct the truth set from a ClinVar VCF; write VCF + labels + provenance."""
    out_dir.mkdir(parents=True, exist_ok=True)
    normalizer = VariantNormalizer(fasta_path)
    counters = _Counters()

    vcf_out = out_dir / "truth_set.vcf"
    labels_out = out_dir / "truth_labels.tsv"

    with (
        _open_maybe_gzip(clinvar_path) as source,
        open(vcf_out, "w", encoding="utf-8") as vcf_fh,
        open(labels_out, "w", encoding="utf-8") as labels_fh,
    ):
        vcf_fh.write("##fileformat=VCFv4.2\n")
        vcf_fh.write(f"##source=build_truth_set.py from {clinvar_path.name}\n")
        vcf_fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        labels_fh.write("chrom\tpos\tref\talt\tgene\tclinvar_class\tlabel\n")

        for line in source:
            if line.startswith("#"):
                continue
            counters.total += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]

            # ClinVar uses bare chromosome names (1, 2, ..., X, Y, MT); the
            # pipeline and the reference FASTA use chr-prefixed names. Prefix
            # here so normalization fetches the right contig and the truth-set
            # keys match pipeline output.
            if not chrom.startswith("chr"):
                chrom = "chrM" if chrom == "MT" else f"chr{chrom}"

            # ClinVar uses '.' for ALT on some records (e.g. CNVs); skip them.
            if alt in (".", "") or "," in alt:
                counters.dropped_reasons["multiallelic_or_no_alt"] = (
                    counters.dropped_reasons.get("multiallelic_or_no_alt", 0) + 1
                )
                continue

            info = _parse_info(cols[7])

            if not _review_status_ok(info):
                counters.wrong_review_status += 1
                continue

            if info.get("CLNSIGCONF"):
                counters.conflicting += 1
                continue

            label = _classify_significance(info)
            if label is None:
                counters.unmapped_significance += 1
                continue

            try:
                norm_chrom, norm_pos, norm_ref, norm_alt = normalizer.normalize(
                    chrom, int(pos), ref, alt
                )
            except (ValueError, TypeError, KeyError):
                counters.normalization_failed += 1
                continue

            gene = _gene_symbol(info)
            clnsig = info.get("CLNSIG", "")

            vcf_fh.write(
                f"{norm_chrom}\t{norm_pos}\t{vid}\t{norm_ref}\t{norm_alt}\t"
                f".\tPASS\tCLNSIG={clnsig};CLNREVSTAT={info.get('CLNREVSTAT', '')}\n"
            )
            labels_fh.write(
                f"{norm_chrom}\t{norm_pos}\t{norm_ref}\t{norm_alt}\t"
                f"{gene}\t{clnsig}\t{label}\n"
            )

            if label == "positive":
                counters.positive += 1
            else:
                counters.negative += 1

    _write_provenance(clinvar_path, out_dir, counters)
    return counters


def _write_provenance(clinvar_path: Path, out_dir: Path, counters: _Counters) -> None:
    """Record the pinned release and exact filters for reproducibility."""
    provenance = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "clinvar_source": clinvar_path.name,
        "reference": "GRCh38",
        "filters": {
            "review_status": sorted(_TRUTH_REVIEW_STATUSES),
            "positive_significance": sorted(_POSITIVE_SIG),
            "negative_significance": sorted(_NEGATIVE_SIG),
            "conflicting_excluded": True,
            "multiallelic_excluded": True,
        },
        "counts": {
            "total_records": counters.total,
            "positive": counters.positive,
            "negative": counters.negative,
            "truth_set_size": counters.positive + counters.negative,
            "rejected_wrong_review_status": counters.wrong_review_status,
            "rejected_conflicting": counters.conflicting,
            "rejected_unmapped_significance": counters.unmapped_significance,
            "rejected_normalization_failed": counters.normalization_failed,
            "dropped_reasons": counters.dropped_reasons,
        },
    }
    (out_dir / "truth_set_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a ClinVar truth set.")
    parser.add_argument(
        "--clinvar", required=True, type=Path, help="Pinned ClinVar GRCh38 VCF."
    )
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path("data/references/GRCh38.primary_assembly.genome.fa"),
        help="Indexed GRCh38 reference FASTA used for left-alignment.",
    )
    args = parser.parse_args(argv)

    if not args.clinvar.exists():
        print(f"error: ClinVar file not found: {args.clinvar}", file=sys.stderr)
        return 1
    if not args.fasta.exists():
        print(f"error: reference FASTA not found: {args.fasta}", file=sys.stderr)
        return 1

    counters = build_truth_set(args.clinvar, args.out_dir, args.fasta)
    size = counters.positive + counters.negative
    print(f"Truth set written to {args.out_dir}")
    print(f"  positive (P/LP): {counters.positive}")
    print(f"  negative (B/LB): {counters.negative}")
    print(f"  total:           {size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
