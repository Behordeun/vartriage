"""Run full annotated pipeline on chr22 with timing and memory tracking."""

import time
import resource
import warnings
from pathlib import Path

from vartriage import (
    Pipeline, PipelineConfig, AnnotationConfig,
    PrioritizationConfig, ReportConfig,
)

warnings.filterwarnings("ignore")

config = PipelineConfig(
    vcf_path=Path("data/giab_chr22.vcf.gz"),
    output_path=Path("data/chr22_annotated.json"),
    annotation=AnnotationConfig(
        gene_annotation_path=Path("data/references/gencode.v46.annotation.gtf"),
        gnomad_path=Path("data/references/gnomad_chr22.tsv"),
        clinvar_path=Path("data/references/clinvar.tsv"),
    ),
    prioritization=PrioritizationConfig(
        max_allele_frequency=0.01,
        cadd_scores_path=Path("data/references/cadd_chr22_sample.tsv"),
    ),
    report=ReportConfig(output_format="json"),
)

pipeline = Pipeline(config)

start = time.time()
result = pipeline.run()
elapsed = time.time() - start

peak_rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
peak_rss_mb = peak_rss_bytes / (1024 * 1024)

print("--- Chr22 Annotated Benchmark ---")
print(f"Wall time:    {elapsed:.1f}s")
print(f"Peak RSS:     {peak_rss_mb:.0f} MB")
print(f"Output:       {result}")
print(f"Output size:  {result.stat().st_size / (1024*1024):.1f} MB")
print(f"Warnings:     {pipeline.warning_accumulator.total_count}")
print(f"Memory OK:    {'YES' if peak_rss_mb < 2048 else 'NO'}")
