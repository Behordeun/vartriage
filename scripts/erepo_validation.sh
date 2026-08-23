#!/bin/bash
# eRepo Validation Pipeline for VarTriage Paper 1
# Runs in ~/Documents/DevProjects/personal_projects/Bioinformatics_Libraries/vartriage
set -euo pipefail

VARTRIAGE_DIR="$HOME/Documents/DevProjects/personal_projects/Bioinformatics_Libraries/vartriage"
cd "$VARTRIAGE_DIR"

REFS="data/references"
OUTPUT="validation_results/erepo"
mkdir -p "$OUTPUT"

echo "=== Step 1: Extract eRepo variant positions from ClinVar ==="
.venv/bin/python3 << 'PYEOF'
import pysam
from pathlib import Path

vcf = pysam.VariantFile("data/references/clinvar.vcf.gz")
positions = []
for rec in vcf.fetch():
    rev = rec.info.get("CLNREVSTAT")
    if rev is None:
        continue
    rev_str = ",".join(rev) if isinstance(rev, tuple) else str(rev)
    if "reviewed_by_expert_panel" not in rev_str and "practice_guideline" not in rev_str:
        continue
    if rec.alts is None or len(rec.alts) == 0:
        continue
    chrom = rec.chrom if rec.chrom.startswith("chr") else f"chr{rec.chrom}"
    positions.append(f"{chrom}\t{rec.pos}\t{rec.ref}\t{rec.alts[0]}")
vcf.close()

Path("validation_results/erepo/erepo_positions.txt").write_text("\n".join(positions) + "\n")
print(f"Extracted {len(positions)} eRepo variant positions")
PYEOF

echo "=== Step 2: Filter REVEL to eRepo positions ==="
echo -e "chrom\tpos\tref\talt\tscore" > "$OUTPUT/revel_erepo.tsv"
awk 'NR==FNR {keys[$1"\t"$2"\t"$3"\t"$4]=1; next} FNR==1{next} ($1"\t"$2"\t"$3"\t"$4) in keys' \
    "$OUTPUT/erepo_positions.txt" \
    "$REFS/revel_genome_wide.tsv" >> "$OUTPUT/revel_erepo.tsv"
echo "Filtered REVEL: $(wc -l < "$OUTPUT/revel_erepo.tsv") rows"

echo "=== Step 3: Create eRepo VCF ==="
.venv/bin/python3 << 'PYEOF'
import pysam
from pathlib import Path

vcf_in = pysam.VariantFile("data/references/clinvar.vcf.gz")
output_path = Path("validation_results/erepo/erepo_variants.vcf.gz")

header = pysam.VariantHeader()
header.add_sample("EREPO")
header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

# Collect contigs
contigs = set()
variants = []
for rec in vcf_in.fetch():
    rev = rec.info.get("CLNREVSTAT")
    if rev is None:
        continue
    rev_str = ",".join(rev) if isinstance(rev, tuple) else str(rev)
    if "reviewed_by_expert_panel" not in rev_str and "practice_guideline" not in rev_str:
        continue
    if rec.alts is None or len(rec.alts) == 0:
        continue
    clnsig = rec.info.get("CLNSIG")
    if clnsig is None:
        continue
    chrom = rec.chrom if rec.chrom.startswith("chr") else f"chr{rec.chrom}"
    contigs.add(chrom)
    variants.append((chrom, rec.pos, rec.ref, rec.alts[0]))
vcf_in.close()

for c in sorted(contigs):
    header.add_line(f'##contig=<ID={c},length=300000000>')

with pysam.VariantFile(str(output_path), "wz", header=header) as out:
    for chrom, pos, ref, alt in sorted(variants):
        rec = out.new_record()
        rec.contig = chrom
        rec.pos = pos
        rec.alleles = (ref, alt)
        rec.samples["EREPO"]["GT"] = (0, 1)
        out.write(rec)

pysam.tabix_index(str(output_path), preset="vcf", force=True)
print(f"Created eRepo VCF: {len(variants)} variants")
PYEOF

echo "=== Step 4: Run VarTriage on eRepo VCF ==="
time .venv/bin/vartriage \
    --vcf "$OUTPUT/erepo_variants.vcf.gz" \
    --gene-annotation "$REFS/gencode.v46.annotation.gtf" \
    --clinvar "$REFS/clinvar.tsv" \
    --revel-scores "$OUTPUT/revel_erepo.tsv" \
    --gnomad-remote gnomad-exomes-v4-grch38 \
    --output "$OUTPUT/erepo_classifications.json" \
    --output-format json \
    --no-confirm

echo "=== Step 5: Compute metrics ==="
.venv/bin/python3 << 'PYEOF'
import json, pysam
from collections import Counter
from pathlib import Path

# Load classifications
data = json.loads(Path("validation_results/erepo/erepo_classifications.json").read_text())
print(f"Classified {len(data)} variants")

# Build expert assertion lookup from ClinVar
vcf = pysam.VariantFile("data/references/clinvar.vcf.gz")
expert = {}
for rec in vcf.fetch():
    rev = rec.info.get("CLNREVSTAT")
    if rev is None:
        continue
    rev_str = ",".join(rev) if isinstance(rev, tuple) else str(rev)
    if "reviewed_by_expert_panel" not in rev_str and "practice_guideline" not in rev_str:
        continue
    if rec.alts is None or len(rec.alts) == 0:
        continue
    clnsig = rec.info.get("CLNSIG")
    if clnsig is None:
        continue
    sig_str = ",".join(clnsig) if isinstance(clnsig, tuple) else str(clnsig)
    if "Pathogenic" in sig_str and "Likely" not in sig_str:
        assertion = "Pathogenic"
    elif "Likely_pathogenic" in sig_str:
        assertion = "Likely_pathogenic"
    elif "Benign" in sig_str and "Likely" not in sig_str:
        assertion = "Benign"
    elif "Likely_benign" in sig_str:
        assertion = "Likely_benign"
    elif "Uncertain" in sig_str:
        assertion = "VUS"
    else:
        continue
    chrom = rec.chrom if rec.chrom.startswith("chr") else f"chr{rec.chrom}"
    expert[(chrom, rec.pos, rec.ref, rec.alts[0])] = assertion
vcf.close()

# Match
matched = []
for v in data:
    key = (v.get("chromosome",""), v.get("position",0), v.get("ref_allele",""), v.get("alt_allele",""))
    if key in expert:
        matched.append({"vt": v.get("acmg_classification","Unknown"), "expert": expert[key],
                        "tags": v.get("evidence_tags",[]), "cons": v.get("functional_consequence","")})

print(f"Matched: {len(matched)}")
exp_path = [m for m in matched if m["expert"] in ("Pathogenic","Likely_pathogenic")]
exp_benign = [m for m in matched if m["expert"] in ("Benign","Likely_benign")]

tp = [m for m in exp_path if m["vt"] in ("Pathogenic","Likely_Pathogenic")]
fp_path = [m for m in exp_benign if m["vt"] in ("Pathogenic","Likely_Pathogenic")]
bn_correct = [m for m in exp_benign if m["vt"] in ("Benign","Likely_Benign")]

sens = len(tp)/len(exp_path) if exp_path else 0
ppv = len(tp)/(len(tp)+len(fp_path)) if (len(tp)+len(fp_path)) else 0
bn_sens = len(bn_correct)/len(exp_benign) if exp_benign else 0

# Per-consequence
cons_stats = {}
for cons in set(m["cons"] for m in exp_path):
    cv = [m for m in exp_path if m["cons"]==cons]
    ct = [m for m in cv if m["vt"] in ("Pathogenic","Likely_Pathogenic")]
    cons_stats[cons] = {"total":len(cv), "tp":len(ct), "sens":len(ct)/len(cv) if cv else 0}

metrics = {
    "vartriage_version":"0.17.2", "total_classified":len(data), "matched":len(matched),
    "expert_plp":len(exp_path), "expert_blb":len(exp_benign),
    "pathogenic_sensitivity":round(sens,4), "pathogenic_ppv":round(ppv,4),
    "pathogenic_tp":len(tp), "pathogenic_fp":len(fp_path),
    "benign_sensitivity":round(bn_sens,4), "benign_tp":len(bn_correct),
    "per_consequence":cons_stats,
    "tag_distribution":dict(Counter(t for m in matched for t in m["tags"]).most_common()),
}
Path("validation_results/erepo/erepo_validation_metrics.json").write_text(json.dumps(metrics, indent=2))

print(f"\n{'='*50}")
print(f"eRepo VALIDATION (v0.17.2)")
print(f"{'='*50}")
print(f"Expert P/LP: {len(exp_path)}")
print(f"Expert B/LB: {len(exp_benign)}")
print(f"Pathogenic sensitivity: {sens*100:.1f}%")
print(f"Pathogenic PPV: {ppv*100:.1f}%")
print(f"Benign sensitivity: {bn_sens*100:.1f}%")
for cons, s in sorted(cons_stats.items(), key=lambda x:-x[1]["total"]):
    print(f"  {cons}: {s['total']} variants, {s['sens']*100:.1f}% sensitivity")
PYEOF

echo "=== DONE ==="
