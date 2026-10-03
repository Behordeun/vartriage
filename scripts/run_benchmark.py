"""Run the classifier over a truth-set VCF and report accuracy metrics.

Driver for the validation-benchmarks runner. Builds a pipeline config from
the given reference files, runs the truth VCF through classification, and
feeds the result to BenchmarkRunner. Writes a JSON metrics summary.

Scope note: with chr22-only gnomAD/CADD this is a HARNESS CHECK, not a
release figure. The release figures come from the full-genome run with
genome-wide references.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from vartriage.models.config import (
    AnnotationConfig,
    PipelineConfig,
    PrioritizationConfig,
)
from vartriage.pipeline import Pipeline
from vartriage.validation.benchmark import BenchmarkReport, BenchmarkRunner


def _report_to_dict(report: BenchmarkReport, scope: str) -> dict:
    def metric(m):
        return {
            "value": round(m.value, 4),
            "ci_low": round(m.low, 4),
            "ci_high": round(m.high, 4),
            "n": m.denominator,
            "underpowered": m.underpowered,
        }

    return {
        "scope": scope,
        "truth_total": report.truth_total,
        "matched": report.matched,
        "unmatched_truth": report.unmatched_truth,
        "metrics": {
            "ppv": metric(report.ppv),
            "sensitivity": metric(report.sensitivity),
            "specificity": metric(report.specificity),
            "npv": metric(report.npv),
        },
        "cohen_kappa": round(report.cohen_kappa, 4),
        "counts": report.counts,
        "confusion": report.confusion,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Benchmark against a truth set.")
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--gnomad", type=Path, default=None)
    parser.add_argument("--revel", type=Path, default=None)
    parser.add_argument("--cadd", type=Path, default=None)
    parser.add_argument("--fasta", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--scope", default="harness")
    args = parser.parse_args(argv)

    for label, path in [
        ("truth VCF", args.truth_vcf),
        ("labels", args.labels),
        ("GTF", args.gtf),
    ]:
        if not path.exists():
            print(f"error: {label} not found: {path}", file=sys.stderr)
            return 1

    config = PipelineConfig(
        vcf_path=args.truth_vcf,
        output_path=args.out.with_suffix(".classified.json"),
        annotation=AnnotationConfig(
            gene_annotation_path=args.gtf,
            gnomad_path=args.gnomad,
            reference_fasta_path=args.fasta,
        ),
        prioritization=PrioritizationConfig(
            revel_scores_path=args.revel,
            cadd_scores_path=args.cadd,
        ),
    )

    pipeline = Pipeline(config)
    runner = BenchmarkRunner(args.labels)
    report = runner.evaluate(pipeline.run_to_classification())

    result = _report_to_dict(report, args.scope)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    ppv = report.ppv
    sens = report.sensitivity
    print(f"[{args.scope}] matched {report.matched}/{report.truth_total}")
    print(
        f"  PPV         {ppv.value:.3f}  95% CI [{ppv.low:.3f}, {ppv.high:.3f}]"
        f"  n={ppv.denominator}"
    )
    print(
        f"  sensitivity {sens.value:.3f}  95% CI [{sens.low:.3f}, {sens.high:.3f}]"
        f"  n={sens.denominator}"
    )
    print(f"  Cohen kappa {report.cohen_kappa:.3f}")
    print(f"  counts      {report.counts}")
    print(f"  written to  {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
