"""Quick QC analysis of vartriage output."""

import json
from collections import Counter
from pathlib import Path

output_path = Path("data/test_output.json")

print(f"Loading {output_path} ...")
with open(output_path) as f:
    data = json.load(f)

print(f"Total variants: {len(data):,}\n")

# Variant type distribution
variant_types = Counter()
for v in data:
    ref, alt = v["ref_allele"], v["alt_allele"]
    if len(ref) == 1 and len(alt) == 1:
        variant_types["SNV"] += 1
    elif len(ref) > len(alt):
        variant_types["Deletion"] += 1
    elif len(alt) > len(ref):
        variant_types["Insertion"] += 1
    else:
        variant_types["Complex"] += 1

print("Variant type distribution:")
for vtype, count in variant_types.most_common():
    print(f"  {vtype}: {count:,} ({100 * count / len(data):.1f}%)")

# Chromosome distribution
chroms = Counter(v["chromosome"] for v in data)
print("Top 5 chromosomes by variant count:")
for chrom, count in chroms.most_common(5):
    print(f"  {chrom}: {count:,}")

# Ti/Tv ratio (transitions vs transversions for SNVs)
transitions = {"A>G", "G>A", "C>T", "T>C"}
ti, tv = 0, 0
for v in data:
    if len(v["ref_allele"]) == 1 and len(v["alt_allele"]) == 1:
        change = f"{v['ref_allele']}>{v['alt_allele']}"
        if change in transitions:
            ti += 1
        else:
            tv += 1

if tv > 0:
    print(f"\nTi/Tv ratio: {ti / tv:.2f} (expected ~2.0-2.1 for WGS)")
    print(f"  Transitions: {ti:,}")
    print(f"  Transversions: {tv:,}")
else:
    print("\nNo SNVs found for Ti/Tv calculation")
