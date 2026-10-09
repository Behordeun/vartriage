"""Derive the PVS1 LoF-mechanism gene set from the ClinGen gene curation list.

Loss of function reaches PVS1 Very Strong only when LoF is an established disease
mechanism for the gene (ClinGen PVS1, Abou Tayoun et al. 2018). ClinGen's dosage
curation encodes that mechanism in the Haploinsufficiency Score column:

  3   Sufficient evidence for dosage pathogenicity  -> established dominant haploinsufficiency
  30  Gene associated with autosomal recessive phenotype -> established recessive LoF mechanism

Both classes are LoF-disease-mechanism genes. Scores 0, 1, 2, and 40 are not
sufficient mechanism evidence and are excluded, so a gene off this list fails
closed to PVS1 Strong rather than Very Strong.
"""

from __future__ import annotations

import sys
from pathlib import Path

LOF_HI_SCORES = {"3", "30"}


def derive(src: Path, dest: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    with src.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            gene, hi_score = cols[0].strip(), cols[4].strip()
            if gene and hi_score in LOF_HI_SCORES and gene not in seen:
                seen.add(gene)
                rows.append((gene, hi_score))
    rows.sort()
    header = (
        "# vartriage PVS1 LoF-mechanism gene set\n"
        "# Derived from the ClinGen gene curation list (ftp.clinicalgenome.org,\n"
        "# ClinGen_gene_curation_list_GRCh38.tsv). A gene is included when its\n"
        "# Haploinsufficiency Score is 3 (sufficient evidence for dominant\n"
        "# haploinsufficiency) or 30 (gene associated with an autosomal recessive\n"
        "# phenotype). Both encode an established loss-of-function disease\n"
        "# mechanism per the ClinGen PVS1 decision tree (Abou Tayoun et al. 2018).\n"
        "# A gene NOT on this list fails closed to PVS1 Strong, never Very Strong.\n"
        "# ClinGen dosage curation is not exhaustive: some genuine LoF-disease\n"
        "# genes (e.g. GAA, ITGA2B, RPGR) are not yet dosage-curated and are\n"
        "# therefore absent here by design, pending a supplementary mechanism source.\n"
        "gene_symbol\thi_score\n"
    )
    dest.write_text(
        header + "\n".join(f"{g}\t{s}" for g, s in rows) + "\n", encoding="utf-8"
    )
    return rows


if __name__ == "__main__":
    src = Path(sys.argv[1])
    dest = Path(sys.argv[2])
    rows = derive(src, dest)
    genes = {g for g, _ in rows}
    print(f"wrote {len(rows)} LoF-mechanism genes to {dest}")
    for check in [
        "BRCA1",
        "BRCA2",
        "CFTR",
        "MSH2",
        "MLH1",
        "MSH6",
        "PAH",
        "GAA",
        "CDH1",
        "RUNX1",
    ]:
        print(f"  {check}: {'present' if check in genes else 'ABSENT'}")
