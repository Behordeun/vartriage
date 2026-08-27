# VCF Quality Control

Compute sample-level quality metrics from a VCF and validate them against expected ranges before annotation runs. QC catches contamination, sample swaps, and variant-caller artifacts early, so a bad sample is flagged before it reaches manual review.

A syntactically valid VCF can still come from a compromised sample. Excess heterozygosity signals contamination, a low Ti/Tv ratio points to caller problems, and a variant count far outside the expected range for the assay hints at a sample swap. QC computes these summary statistics in a single streaming pass and reports a PASS/WARN/FAIL verdict per metric.

## Quick start

Standalone QC check, no annotation:

```bash
vartriage qc --vcf sample.vcf.gz --sample SAMPLE1 --assay-type wes
```

QC as a gate on the full pipeline:

```bash
# Halt before annotation if any metric fails
vartriage --vcf sample.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv \
  --assay-type wgs --strict-qc

# Skip QC entirely for a pre-validated file
vartriage --vcf sample.vcf.gz --output results.json \
  --gene-annotation gencode.gtf --gnomad gnomad.tsv --skip-qc
```

## Metrics

| Metric | What it measures | Signal when abnormal |
| ------ | ---------------- | -------------------- |
| Ti/Tv ratio | Transitions (A-G, C-T and their reverse) over transversions across biallelic SNVs | A low ratio indicates false-positive SNVs from a miscalibrated caller |
| Het/Hom ratio | Heterozygous (0/1) over homozygous-alt (1/1) calls for the target sample | High ratio signals contamination; low ratio signals consanguinity or a sample issue |
| Total variants | Count of all records | Far outside the assay range points to a sample swap or wrong reference |
| Ins/Del ratio | Insertions over deletions | Extreme deviation flags systematic caller bias |
| Per-chromosome counts | Variant count per contig | Reported for downstream density inspection |

Het/Hom is only computed when a sample is available. Pass `--sample` for multi-sample VCFs; a single-sample VCF is detected automatically. When a VCF has no indels or no genotype calls, the affected ratio is reported but not counted against the verdict.

## Thresholds

Each assay type carries its own expected ranges. A metric is WARN when it falls outside the warn range and FAIL when it falls outside the wider fail range.

| Assay | Ti/Tv warn | Ti/Tv fail | Het/Hom warn | Het/Hom fail | Variant count (warn) |
| ----- | ---------- | ---------- | ------------ | ------------ | -------------------- |
| `wgs` | 1.8 to 2.5 | 1.5 to 3.0 | 1.0 to 3.0 | 0.5 to 5.0 | 3.5M to 6.0M |
| `wes` | 1.8 to 2.5 | 1.5 to 3.0 | 1.0 to 3.0 | 0.5 to 5.0 | 20K to 150K |
| `panel` | 1.5 to 3.0 | 1.2 to 3.5 | 0.8 to 4.0 | 0.3 to 6.0 | no constraint |

Ins/Del ratio warn range is 0.5 to 1.5 for wgs and wes (fail 0.3 to 2.0), and 0.3 to 2.0 for panel (fail 0.2 to 3.0).

Variant count FAIL triggers below 0.5x the warn minimum or above 2x the warn maximum. Panel assays apply no count constraint because panel sizes vary widely.

Override the warn ranges from the CLI:

```bash
vartriage qc --vcf sample.vcf.gz --sample S1 --assay-type wes \
  --expected-titv 1.9,2.3 --expected-het-hom 1.2,2.5
```

Or from a TOML config under `~/.vartriage/config.toml`:

```toml
[qc]
expected_titv = [1.9, 2.3]
expected_het_hom = [1.2, 2.5]
```

CLI values take precedence over TOML values.

## Verdicts

The overall verdict is the worst status across all metrics:

- **PASS**: every metric is within its warn range
- **WARN**: at least one metric is borderline, none critical
- **FAIL**: at least one metric is critically out of range

Without `--strict-qc`, a WARN or FAIL is printed to stderr and the pipeline proceeds. With `--strict-qc`, a FAIL halts the pipeline before annotation and exits with code 3, naming the failing metrics.

## Output

QC prints a table to stderr before annotation begins:

```text
======================================================================
  VCF Quality Control  |  Assay: WES
======================================================================
  Metric                      Value         Expected     Status
  --------------------------------------------------------------
  Ti/Tv Ratio                  2.03          1.8-2.5     ✓ PASS
  Total Variants             44,747   20,000-150,000     ✓ PASS
  Ins/Del Ratio                1.01          0.5-1.5     ✓ PASS
  Het/Hom Ratio                1.62          1.0-3.0     ✓ PASS
  --------------------------------------------------------------
  Overall QC Verdict: ✓ PASS
======================================================================
```

The `qc` subcommand can also write a machine-readable report with `--output-json`:

```bash
vartriage qc --vcf sample.vcf.gz --sample S1 --assay-type wes \
  --output-json qc_report.json
```

The JSON contains the raw metrics, per-metric check results with expected ranges, and the overall verdict.

## CLI reference

### Standalone QC

```bash
vartriage qc [OPTIONS]
```

| Flag | Default | Description |
| ---- | ------- | ----------- |
| `--vcf` | required | Input VCF (`.vcf` or `.vcf.gz`) |
| `--sample` | None | Sample name for het/hom extraction |
| `--assay-type` | wes | `wgs`, `wes`, or `panel` |
| `--output-json` | None | Write the QC report as JSON to this path |
| `--strict` | false | Exit with code 3 if any metric reaches FAIL |
| `--expected-titv` | None | Override Ti/Tv warn range as `MIN,MAX` |
| `--expected-het-hom` | None | Override het/hom warn range as `MIN,MAX` |

### QC on the full pipeline

| Flag | Default | Description |
| ---- | ------- | ----------- |
| `--assay-type` | wes | Assay type for threshold selection |
| `--strict-qc` | false | Halt with exit code 3 on any FAIL before annotation |
| `--skip-qc` | false | Skip QC entirely |
| `--expected-titv` | None | Override Ti/Tv warn range as `MIN,MAX` |
| `--expected-het-hom` | None | Override het/hom warn range as `MIN,MAX` |

## Python API

```python
from pathlib import Path

from vartriage.qc.config import QCConfig
from vartriage.qc.metrics import compute_qc_metrics
from vartriage.qc.validator import QCStatus, QCValidator

metrics = compute_qc_metrics(Path("sample.vcf.gz"), sample_id="SAMPLE1")

config = QCConfig(assay_type="wes", strict=True)
report = QCValidator(config).validate(metrics)

print(report.overall_status)  # QCStatus.PASS / WARN / FAIL
for check in report.checks:
    print(check.metric_name, check.status.value, check.message)

if report.overall_status is QCStatus.FAIL:
    ...  # handle the failing sample
```

Run QC as part of a pipeline by setting the `qc` field on `PipelineConfig`:

```python
from vartriage.models.config import PipelineConfig
from vartriage.pipeline import Pipeline
from vartriage.qc.config import QCConfig

config = PipelineConfig(
    vcf_path=Path("sample.vcf.gz"),
    output_path=Path("results.json"),
    qc=QCConfig(assay_type="wgs", strict=True),
)
pipeline = Pipeline(config)
pipeline.run()

report = pipeline.qc_report  # QCReport from the pre-flight pass
```

## Performance

QC runs as a single streaming pass over the VCF with no annotation lookups and no random access. Memory use is bounded by the number of contigs. The pass is independent of the annotation pass, so a QC-gated run reads the file twice.

## Configuration reference

### QCConfig

| Field | Type | Default | Notes |
| ----- | ---- | ------- | ----- |
| `assay_type` | str | "wes" | `wgs`, `wes`, or `panel` |
| `strict` | bool | False | FAIL halts the pipeline with exit code 3 |
| `skip` | bool | False | Bypass QC entirely |
| `expected_ti_tv` | tuple | None | Override Ti/Tv warn range `(min, max)` |
| `expected_het_hom` | tuple | None | Override het/hom warn range `(min, max)` |
| `sample_id` | str | None | Sample for het/hom; single-sample VCFs auto-detect |

Raises `ValueError` when `assay_type` is unrecognized or an override range has `min >= max`.
