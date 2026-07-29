#!/usr/bin/env python3
"""Validate vartriage gene-disease linkage with a simulated affected patient.

Takes the healthy GIAB HG002 chr22 VCF and spikes in known pathogenic
NF2 variants from ClinVar, creating a realistic "needle in haystack"
scenario. Then runs the full pipeline with phenotype terms matching
neurofibromatosis type 2.

Clinical scenario:
    Patient presents with bilateral vestibular schwannomas, hearing loss,
    and multiple meningiomas — classic NF2 presentation. Whole exome
    sequencing reveals ~50k variants on chr22. The pipeline should surface
    the causative NF2 nonsense/frameshift variants at the top of the
    prioritized list.

Usage:
    python scripts/validate_affected_patient.py
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

# Pathogenic NF2 variants from ClinVar to spike into the healthy VCF.
# These are real variants reported in neurofibromatosis type 2 patients.
# Format: (pos, ref, alt, description)
NF2_SPIKE_VARIANTS = [
    # Frameshift deletions — classic loss-of-function
    (29604004, "CG", "C", "NF2 exon1 frameshift del"),
    (29604027, "GC", "G", "NF2 exon1 frameshift del c.76delC"),
    (29604040, "CA", "C", "NF2 exon1 frameshift del c.89delA"),
    # Nonsense (stop gain)
    (29604041, "A", "T", "NF2 exon1 nonsense c.90A>T p.Arg30*"),
    (29604056, "A", "T", "NF2 exon1 nonsense"),
]

# HPO terms for neurofibromatosis type 2
NF2_HPO_TERMS = [
    "HP:0009588",  # Vestibular schwannoma (bilateral)
    "HP:0002858",  # Meningioma
    "HP:0000365",  # Hearing impairment
    "HP:0009592",  # Astrocytoma (spinal)
]


def build_spike_vcf(source_vcf: Path, output_vcf: Path) -> int:
    """Inject pathogenic NF2 variants into the GIAB VCF.

    Returns the number of spiked variants.
    """
    import pysam

    vcf_in = pysam.VariantFile(str(source_vcf))
    try:
        header = vcf_in.header.copy()

        # Add INFO field for spike-in tracking
        header.add_line(
            '##INFO=<ID=SPIKE,Number=0,Type=Flag,Description="Spiked-in pathogenic variant">'
        )

        vcf_out = pysam.VariantFile(str(output_vcf), "wz", header=header)
        try:
            spike_positions = {pos for pos, _, _, _ in NF2_SPIKE_VARIANTS}
            spiked_count = 0
            spike_inserted = False

            for rec in vcf_in:
                # Insert spike variants just before we pass their position
                if not spike_inserted and rec.pos >= min(spike_positions):
                    for pos, ref, alt, desc in sorted(NF2_SPIKE_VARIANTS):
                        new_rec = vcf_out.new_record()
                        new_rec.contig = "chr22"
                        new_rec.pos = pos
                        new_rec.alleles = (ref, alt)
                        new_rec.qual = 99
                        new_rec.filter.add("PASS")
                        new_rec.info["SPIKE"] = True
                        # Set genotype: heterozygous (typical for AD condition)
                        new_rec.samples["HG002"]["GT"] = (0, 1)
                        vcf_out.write(new_rec)
                        spiked_count += 1
                    spike_inserted = True

                vcf_out.write(rec)
        finally:
            vcf_out.close()
    finally:
        vcf_in.close()

    # Index the output
    pysam.tabix_index(str(output_vcf), preset="vcf", force=True)

    return spiked_count


def update_knowledge_base(knowledge_dir: Path) -> None:
    """Add NF2 to the knowledge base TSV files for this validation run."""
    omim_path = knowledge_dir / "omim_gene_disease.tsv"
    content = omim_path.read_text()
    if "NF2" not in content:
        with open(omim_path, "a") as f:
            f.write("NF2\tNeurofibromatosis type 2\t101000\tAD\n")
            f.write("NF2\tMeningioma, familial\t607174\tAD\n")
            f.write("NF2\tSchwannomatosis\t162091\tAD\n")

    hpo_path = knowledge_dir / "hpo_gene_annotations.tsv"
    content = hpo_path.read_text()
    if "NF2" not in content:
        with open(hpo_path, "a") as f:
            f.write("NF2\tHP:0009588;HP:0002858;HP:0000365;HP:0009592\n")

    validity_path = knowledge_dir / "clingen_validity.tsv"
    content = validity_path.read_text()
    if "NF2" not in content:
        with open(validity_path, "a") as f:
            f.write("NF2\tDefinitive\n")

    constraint_path = knowledge_dir / "gnomad_constraint.tsv"
    content = constraint_path.read_text()
    if "NF2" not in content:
        with open(constraint_path, "a") as f:
            f.write("NF2\t1.00\t0.12\t3.21\n")

    action_path = knowledge_dir / "clingen_actionability.tsv"
    content = action_path.read_text()
    if "NF2" not in content:
        with open(action_path, "a") as f:
            f.write("NF2\tsurveillance\n")


def run_pipeline(
    vcf_path: Path,
    output_path: Path,
    references_dir: Path,
    knowledge_dir: Path,
) -> Path:
    """Run vartriage with gene-disease linkage on the spiked VCF."""
    from vartriage.knowledge.config import KnowledgeBaseConfig
    from vartriage.models.config import (
        AnnotationConfig,
        PipelineConfig,
        PrioritizationConfig,
        ReportConfig,
    )
    from vartriage.pipeline import Pipeline

    hpo_terms = frozenset(NF2_HPO_TERMS)

    print(f"\n{'='*70}")
    print("RUNNING PIPELINE")
    print(f"{'='*70}")
    print(f"  VCF: {vcf_path.name}")
    print(f"  HPO terms: {', '.join(NF2_HPO_TERMS)}")
    print(f"  Knowledge dir: {knowledge_dir}")

    annotation_config = AnnotationConfig(
        gene_annotation_path=references_dir / "gencode_chr22.gtf",
        gnomad_path=references_dir / "gnomad_chr22.tsv",
        clinvar_path=references_dir / "clinvar.tsv",
    )

    knowledge_config = KnowledgeBaseConfig(
        data_dir=knowledge_dir,
        hpo_terms=hpo_terms,
    )

    pipeline_config = PipelineConfig(
        vcf_path=vcf_path,
        output_path=output_path,
        annotation=annotation_config,
        prioritization=PrioritizationConfig(),
        report=ReportConfig(output_format="json"),
        knowledge=knowledge_config,
    )

    pipeline = Pipeline(pipeline_config)
    result_path = pipeline.run()
    print(f"  Output: {result_path}")
    return result_path


def _print_check(label_pass: str, label_fail: str, passed: bool) -> int:
    """Print a PASS/FAIL check line and return 1 if passed, 0 otherwise."""
    print(f"  [{'PASS' if passed else 'FAIL'}] {label_pass if passed else label_fail}")
    return 1 if passed else 0


def _print_variant_detail(v: dict) -> None:  # type: ignore[type-arg]
    """Print a single pathogenic variant's details."""
    gene = v.get("gene_name", "Intergenic")
    constraint = v.get("gene_constraint")
    diseases = v.get("disease_associations", [])
    print(f"\n  {gene} {v.get('chromosome')}:{v.get('position')}")
    print(f"    Consequence: {v.get('functional_consequence')}")
    print(f"    ClinVar: {v.get('clinvar_assertion')}")
    print(f"    Phenotype match: {v.get('phenotype_match_score')}")
    print(f"    Actionable: {v.get('is_actionable')}")
    if constraint:
        print(f"    Constraint: pLI={constraint['pli']}, LOEUF={constraint['loeuf']}")
    for d in diseases[:3]:
        print(f"    Disease: {d['disease_name']} ({d['inheritance_mode']})")


def _format_nf2_summary(v: dict) -> str:  # type: ignore[type-arg]
    """Format a single NF2 variant as a compact one-line summary."""
    return (
        f"chr22:{v.get('position')} "
        f"{v.get('ref_allele')}>{v.get('alt_allele')} | "
        f"{v.get('functional_consequence')} | "
        f"{v.get('acmg_classification')} | "
        f"ClinVar:{v.get('clinvar_assertion')} | "
        f"HPO:{v.get('phenotype_match_score')}"
    )


def _print_classification_distribution(results: list) -> None:  # type: ignore[type-arg]
    classifications: dict[str, int] = {}
    for r in results:
        cls = r.get("acmg_classification", "Unknown")
        classifications[cls] = classifications.get(cls, 0) + 1
    print("Classification distribution:")
    for cls, count in sorted(classifications.items()):
        print(f"  {cls}: {count}")


def _print_pathogenic_findings(results: list) -> None:  # type: ignore[type-arg]
    pathogenic = [
        r for r in results
        if r.get("acmg_classification") in ("Pathogenic", "Likely_Pathogenic")
    ]
    print(f"\n{'='*70}")
    print(f"PATHOGENIC / LIKELY PATHOGENIC FINDINGS: {len(pathogenic)}")
    print(f"{'='*70}")
    for v in pathogenic[:10]:
        _print_variant_detail(v)


def _print_nf2_detail(nf2_variants: list) -> None:  # type: ignore[type-arg]
    if not nf2_variants:
        return
    print(f"\n{'='*70}")
    print("NF2 VARIANTS DETAIL")
    print(f"{'='*70}")
    for v in nf2_variants:
        print(f"  {_format_nf2_summary(v)}")


def _build_nf2_checks(nf2_variants: list) -> list[tuple[str, str, bool]]:  # type: ignore[type-arg]
    nf2_with_disease = [v for v in nf2_variants if v.get("disease_associations")]
    nf2_with_pheno = [
        v for v in nf2_variants
        if v.get("phenotype_match_score") and v["phenotype_match_score"] > 0
    ]
    nf2_actionable = [v for v in nf2_variants if v.get("is_actionable")]
    nf2_constrained = [v for v in nf2_variants if v.get("gene_constraint")]

    score_detail = f" ({nf2_with_pheno[0]['phenotype_match_score']:.2f})" if nf2_with_pheno else ""
    pli_detail = f" (pLI={nf2_constrained[0]['gene_constraint']['pli']})" if nf2_constrained else ""

    return [
        ("NF2 variants detected in output",
         "No NF2 variants found",
         len(nf2_variants) > 0),
        ("NF2 variants have disease associations attached",
         "NF2 variants missing disease associations",
         bool(nf2_with_disease)),
        (f"NF2 variants have phenotype match score > 0{score_detail}",
         "NF2 variants have zero phenotype match",
         bool(nf2_with_pheno)),
        ("NF2 variants flagged as actionable",
         "NF2 variants not flagged actionable",
         bool(nf2_actionable)),
        (f"NF2 constraint metrics present{pli_detail}",
         "NF2 missing constraint metrics",
         bool(nf2_constrained)),
    ]


def analyze_results(output_path: Path) -> None:
    """Analyze pipeline output and validate NF2 variants surfaced correctly."""
    with open(output_path) as f:
        results = json.load(f)

    print(f"\n{'='*70}")
    print("RESULTS ANALYSIS")
    print(f"{'='*70}")
    print(f"Total classified variants: {len(results)}")

    with_context = [r for r in results if "disease_associations" in r]
    print(f"Variants with gene-disease context: {len(with_context)}")

    nf2_variants = [r for r in results if r.get("gene_name") == "NF2"]
    print(f"NF2 variants found: {len(nf2_variants)}")

    _print_classification_distribution(results)
    _print_pathogenic_findings(results)
    _print_nf2_detail(nf2_variants)

    checks = _build_nf2_checks(nf2_variants)

    print(f"\n{'='*70}")
    print("VALIDATION")
    print(f"{'='*70}")

    checks_passed = sum(
        _print_check(pass_msg, fail_msg, condition)
        for pass_msg, fail_msg, condition in checks
    )

    print(f"\n  Result: {checks_passed}/{len(checks)} checks passed")

    if checks_passed == len(checks):
        print("\n  Gene-disease linkage validation PASSED")
    else:
        print("\n  Gene-disease linkage validation FAILED")
        sys.exit(1)


def main() -> None:
    """Run the full spike-in validation pipeline."""
    # Suppress noisy missing-data warnings during batch processing
    warnings.filterwarnings("ignore", category=UserWarning)

    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data"
    references_dir = data_dir / "references"
    source_vcf = data_dir / "giab_chr22.vcf.gz"

    if not source_vcf.exists():
        print(f"Error: source VCF not found: {source_vcf}")
        sys.exit(1)

    with tempfile.TemporaryDirectory(prefix="vartriage_validate_") as tmp:
        tmp_path = Path(tmp)

        # Copy knowledge base to temp dir so we can add NF2
        knowledge_dir = tmp_path / "knowledge"
        shutil.copytree(
            project_root / "vartriage" / "data" / "knowledge",
            knowledge_dir,
        )
        update_knowledge_base(knowledge_dir)

        # Build spike-in VCF
        spiked_vcf = tmp_path / "patient_nf2.vcf.gz"
        print(f"{'='*70}")
        print("SPIKE-IN VCF CONSTRUCTION")
        print(f"{'='*70}")
        print(f"Source: {source_vcf} (HG002, healthy)")
        print(f"Target: {spiked_vcf}")
        print(f"Injecting {len(NF2_SPIKE_VARIANTS)} pathogenic NF2 variants...")

        spiked_count = build_spike_vcf(source_vcf, spiked_vcf)
        print(f"Spiked {spiked_count} variants into {spiked_vcf.name}")

        # Run pipeline
        output_path = tmp_path / "results.json"
        actual_output = run_pipeline(
            vcf_path=spiked_vcf,
            output_path=output_path,
            references_dir=references_dir,
            knowledge_dir=knowledge_dir,
        )

        # Analyze results
        analyze_results(actual_output)


if __name__ == "__main__":
    main()
