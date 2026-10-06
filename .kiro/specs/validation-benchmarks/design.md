# Validation and Benchmarks (v1.0.0) - Design

## Validation Package

```
vartriage/validation/              # or scripts/validation/
├── __init__.py
├── truth_set.py                   # Download, filter, normalize ClinVar truth set
├── benchmark.py                   # BenchmarkRunner: run pipeline, compute metrics
├── metrics.py                     # Sensitivity, specificity, kappa, Wilson CI
├── consequence_comparison.py      # Local vs API mode consequence concordance
└── report_generator.py            # Generate validation_report.md from results
```

## BenchmarkRunner

```python
class BenchmarkRunner:
    def __init__(self, truth_vcf: Path, truth_labels: Path,
                 pipeline_config: PipelineConfig) -> None:
        self._truth_vcf = truth_vcf
        self._labels = self._load_labels(truth_labels)
        self._config = pipeline_config

    def run(self) -> BenchmarkReport:
        # 1. Run pipeline on truth VCF
        pipeline = Pipeline(self._config)
        pipeline.run()

        # 2. Load pipeline output
        predictions = self._load_predictions(self._config.output_path)

        # 3. Match predictions to truth labels by coordinate
        matched = self._match(predictions, self._labels)

        # 4. Compute metrics
        return self._compute_metrics(matched)


@dataclass
class BenchmarkReport:
    total_variants: int
    matched_variants: int
    sensitivity: MetricWithCI      # TP / (TP + FN)
    specificity: MetricWithCI      # TN / (TN + FP)
    ppv: MetricWithCI              # TP / (TP + FP)
    npv: MetricWithCI              # TN / (TN + FN)
    cohens_kappa: MetricWithCI
    confusion_matrix: dict[str, dict[str, int]]  # predicted -> truth -> count
    per_criterion: dict[str, CriterionMetrics]
    consequence_concordance: float | None  # if comparison run


@dataclass
class MetricWithCI:
    value: float
    ci_lower: float
    ci_upper: float
    n: int
```

## Wilson Confidence Interval

```python
def wilson_ci(successes: int, total: int, confidence: float = 0.95
              ) -> tuple[float, float]:
    """Wilson score interval for binomial proportion.

    More accurate than normal approximation for small samples
    or proportions near 0 or 1.
    """
    from math import sqrt
    z = 1.96 if confidence == 0.95 else ...
    p_hat = successes / total
    denom = 1 + z**2 / total
    center = (p_hat + z**2 / (2 * total)) / denom
    margin = z * sqrt((p_hat * (1 - p_hat) + z**2 / (4 * total)) / total) / denom
    return (center - margin, center + margin)
```

## Truth Set Construction

```python
class TruthSetBuilder:
    def __init__(self, clinvar_vcf_url: str, reference_fasta: Path,
                 output_dir: Path) -> None: ...

    def build(self) -> tuple[Path, Path]:
        """Download, filter, normalize. Returns (truth.vcf, labels.tsv)."""
        # 1. Download ClinVar VCF
        # 2. Filter: review_status >= "reviewed by expert panel"
        # 3. Filter: no conflicting interpretations
        # 4. Split: CLNSIG in (Pathogenic, Likely_pathogenic) -> positive
        #           CLNSIG in (Benign, Likely_benign) -> negative
        # 5. Normalize all variants
        # 6. Write filtered VCF + labels TSV
        ...
```

## Consequence Comparison

```python
class ConsequenceComparison:
    """Compare local vs API consequence calls for the same variant set."""

    def run(self, vcf_path: Path, fasta_path: Path) -> ComparisonReport:
        # 1. Run local mode with codon resolver
        # 2. Run API mode with VEP
        # 3. Pair results by coordinate
        # 4. Tabulate: concordant, local-only-correct, api-only-correct, both-wrong
        ...

@dataclass
class ComparisonReport:
    total: int
    concordant: int
    concordance_rate: float
    discordant_by_type: dict[str, int]  # e.g., "local=Synonymous,api=Missense": 15
```

## Validation Report Template

```markdown
# VarTriage Validation Report

## Truth Set
- Source: ClinVar {release_date}
- Filter: review_status >= "reviewed by expert panel"
- Positive (P/LP): {n_positive} variants
- Negative (B/LB): {n_negative} variants

## Classification Accuracy
| Metric | Value | 95% CI |
|--------|-------|--------|
| Sensitivity | {sens} | [{ci_l}, {ci_u}] |
| Specificity | {spec} | [{ci_l}, {ci_u}] |
| PPV | {ppv} | [{ci_l}, {ci_u}] |
| NPV | {npv} | [{ci_l}, {ci_u}] |
| Cohen's Kappa | {kappa} | [{ci_l}, {ci_u}] |

## Confusion Matrix
...

## Per-Criterion Analysis
...
```
