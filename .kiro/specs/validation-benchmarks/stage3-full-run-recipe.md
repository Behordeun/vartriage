# Stage 3 — full-genome release run: staging requirement and recipe

The benchmark code (truth set, `BenchmarkRunner`, `scripts/run_benchmark.py`)
is complete and proven on the chr22 harness. The full-genome release run is
**data-gated**, not code-gated: it needs genome-wide reference data that is
not yet materialized locally. This file records exactly what to stage and
how to run, so the run is resumable without re-deriving anything.

## Reference data status (checked 2026-10-03)

| Source | Needed for | Status |
| --- | --- | --- |
| ClinVar `clinvar_20260928` | truth set | present, pinned |
| GRCh38 FASTA + GENCODE v46 | consequence, normalization | present (`data/references/`) |
| REVEL genome-wide | PP3 / BP4 (primary computational) | OneDrive placeholder, **hydrating** (`Bioref/revel/revel_genome_wide.tsv`, blocks=0 -> dehydrated) |
| gnomAD constraint + knowledge | PVS1 strength, PM1, gene-disease | present (`data/references/gnomad_constraint_full.tsv`, `knowledge_full/`) |
| gnomAD allele frequencies, genome-wide | PM2 / BA1 / BS1 / BS2 | **chr22 only** — needs genome-wide staging |
| CADD, genome-wide | BP4 / PP3 fallback when REVEL absent | **chr22 only** — needs genome-wide staging |

## Staging steps (hours-to-overnight; runway is a week)

1. Hydrate REVEL: `cat Bioref/revel/revel_genome_wide.tsv > /dev/null` (2 GB
   OneDrive pull) — triggered in the background on 2026-10-03.
2. gnomAD v4.1 genome-wide allele frequencies: either download the per-chrom
   gnomAD VCFs (large, ~500 GB exomes) OR use the remote-tabix backend
   (`--remote`) to range-query gnomAD over HTTP without a full download — the
   remote backend already exists and is the pragmatic path for a 18 K-variant
   truth set, since only the truth coordinates are queried.
3. CADD genome-wide SNV scores (optional, REVEL fallback only): download the
   CADD v1.7 whole-genome SNV TSV + tabix index, or leave absent and let the
   run record CADD as a missing source (BP4 falls back cleanly).

## Run recipe (genome-wide, once staged)

```
uv run --extra all python scripts/run_benchmark.py \
  --truth-vcf validation_results/truth_set/truth_set.vcf \
  --labels   validation_results/truth_set/truth_labels.tsv \
  --gtf      data/references/gencode.v46.annotation.gtf \
  --fasta    data/references/GRCh38.primary_assembly.genome.fa \
  --revel    <genome-wide REVEL tsv> \
  --gnomad   <genome-wide gnomAD tsv, or use --remote in a config> \
  --out      validation_results/benchmark_full/metrics_full.json \
  --scope    full-genome-release
```

The full run is 18,166 variants and is externally-latency-bound if gnomAD is
queried remotely; run it detached and poll the output JSON, do not block a
single turn on it.

## Interpreting the result

Apply `accuracy-decision-rule.md` exactly: judge PPV on its Wilson CI lower
bound, report every metric with its CI, record any absent source (gnomAD-freq,
CADD) as measured-missing per FR-7, and position the release by the committed
ladder. Do not retune thresholds against the truth set.
