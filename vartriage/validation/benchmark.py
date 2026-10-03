"""Benchmark the classifier against an expert-curated truth set.

Implements the benchmark-runner stage of the validation-benchmarks spec:
match classifications to truth labels by normalized coordinate and compute
sensitivity, specificity, PPV, NPV (each with a Wilson 95% confidence
interval), Cohen's kappa, and a 5x5 confusion matrix.

Per the pre-registered decision rule, the primary figure is PPV on the
P/LP call, judged on its CI lower bound. A VUS prediction is a non-call:
on a positive truth variant it is a miss (false negative for sensitivity),
on a negative truth variant it is neither a true negative nor a false
positive, and it is always counted and reported separately.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from vartriage.models.variant import ACMGClassification, ClassifiedVariant

_POSITIVE = "positive"
_NEGATIVE = "negative"

_PATHOGENIC_CALL: frozenset[ACMGClassification] = frozenset(
    {ACMGClassification.PATHOGENIC, ACMGClassification.LIKELY_PATHOGENIC}
)
_BENIGN_CALL: frozenset[ACMGClassification] = frozenset(
    {ACMGClassification.BENIGN, ACMGClassification.LIKELY_BENIGN}
)

_TIERS: tuple[ACMGClassification, ...] = (
    ACMGClassification.PATHOGENIC,
    ACMGClassification.LIKELY_PATHOGENIC,
    ACMGClassification.VUS,
    ACMGClassification.LIKELY_BENIGN,
    ACMGClassification.BENIGN,
)

_WILSON_Z = 1.959963984540054  # z for a two-sided 95% interval


@dataclass(frozen=True)
class MetricWithCI:
    """A proportion with its Wilson 95% confidence interval and sample size."""

    name: str
    value: float
    low: float
    high: float
    numerator: int
    denominator: int

    @property
    def underpowered(self) -> bool:
        """True when the CI half-width exceeds 5 points (per the decision rule)."""
        return (self.high - self.low) / 2.0 > 0.05


def wilson_ci(successes: int, total: int) -> tuple[float, float, float]:
    """Point estimate and Wilson 95% CI for a binomial proportion.

    Returns (point, low, high). For total == 0 returns (0, 0, 0); the caller
    reports such a stratum as having no denominator, not a real proportion.
    """
    if total == 0:
        return (0.0, 0.0, 0.0)
    p = successes / total
    z = _WILSON_Z
    z2 = z * z
    denom = 1.0 + z2 / total
    center = (p + z2 / (2 * total)) / denom
    margin = z * math.sqrt(p * (1 - p) / total + z2 / (4 * total * total)) / denom
    return (p, max(0.0, center - margin), min(1.0, center + margin))


def _metric(name: str, successes: int, total: int) -> MetricWithCI:
    point, low, high = wilson_ci(successes, total)
    return MetricWithCI(name, point, low, high, successes, total)


@dataclass
class BenchmarkReport:
    """Collected benchmark metrics for one run against a truth set."""

    truth_total: int
    matched: int
    unmatched_truth: int
    sensitivity: MetricWithCI
    specificity: MetricWithCI
    ppv: MetricWithCI
    npv: MetricWithCI
    cohen_kappa: float
    confusion: dict[str, dict[str, int]]
    counts: dict[str, int] = field(default_factory=dict)


def _load_labels(labels_path: Path) -> dict[tuple[str, int, str, str], str]:
    """Load truth_labels.tsv into a coordinate-keyed label map."""
    labels: dict[tuple[str, int, str, str], str] = {}
    with open(labels_path, encoding="utf-8") as handle:
        handle.readline()  # header
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            labels[(cols[0], int(cols[1]), cols[2], cols[3])] = cols[6]
    return labels


def _coordinate_key(classified: ClassifiedVariant) -> tuple[str, int, str, str]:
    v = classified.scored.annotated.variant
    return (v.chrom, v.pos, v.ref, v.alt)


def _truth_tier(label: str) -> str:
    """Represent a binary truth label at its canonical five-tier endpoint."""
    return (
        ACMGClassification.PATHOGENIC.value
        if label == _POSITIVE
        else ACMGClassification.BENIGN.value
    )


class BenchmarkRunner:
    """Run the classifier over a truth set and compute accuracy metrics."""

    def __init__(self, labels_path: Path) -> None:
        self._labels = _load_labels(labels_path)

    def evaluate(
        self, classified_stream: Iterator[ClassifiedVariant]
    ) -> BenchmarkReport:
        """Match a classified-variant stream to truth labels and score it.

        Coordinate matching uses the normalized chrom:pos:ref:alt key, the
        same normalization the truth set was built with, so a key that
        misses is a genuine absence rather than a representation mismatch.
        """
        tp = fp = tn = fn = 0
        vus_on_positive = vus_on_negative = 0
        matched = 0
        seen: set[tuple[str, int, str, str]] = set()
        confusion: dict[str, dict[str, int]] = {
            t.value: {u.value: 0 for u in _TIERS} for t in _TIERS
        }
        kappa_pairs: list[tuple[str, str]] = []

        for classified in classified_stream:
            key = _coordinate_key(classified)
            label = self._labels.get(key)
            if label is None:
                continue
            if key in seen:
                continue  # one classification per truth coordinate
            matched += 1
            seen.add(key)

            predicted = classified.classification
            truth_tier = _truth_tier(label)
            confusion[predicted.value][truth_tier] += 1
            kappa_pairs.append((predicted.value, truth_tier))

            tp, fp, tn, fn, vus_on_positive, vus_on_negative = _update_counts(
                label, predicted, tp, fp, tn, fn, vus_on_positive, vus_on_negative
            )

        return BenchmarkReport(
            truth_total=len(self._labels),
            matched=matched,
            unmatched_truth=len(self._labels) - len(seen),
            sensitivity=_metric("sensitivity", tp, tp + fn),
            specificity=_metric("specificity", tn, tn + fp),
            ppv=_metric("ppv", tp, tp + fp),
            npv=_metric("npv", tn, tn + fn),
            cohen_kappa=_cohen_kappa(kappa_pairs),
            confusion=confusion,
            counts={
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "vus_on_positive": vus_on_positive,
                "vus_on_negative": vus_on_negative,
            },
        )


def _update_counts(
    label: str,
    predicted: ACMGClassification,
    tp: int,
    fp: int,
    tn: int,
    fn: int,
    vus_on_positive: int,
    vus_on_negative: int,
) -> tuple[int, int, int, int, int, int]:
    """Update TP/FP/TN/FN and VUS counts for one matched variant."""
    is_path = predicted in _PATHOGENIC_CALL
    is_benign = predicted in _BENIGN_CALL
    if label == _POSITIVE:
        if is_path:
            tp += 1
        else:
            fn += 1
            if not is_benign:
                vus_on_positive += 1
    else:
        if is_path:
            fp += 1
        elif is_benign:
            tn += 1
        else:
            vus_on_negative += 1
    return tp, fp, tn, fn, vus_on_positive, vus_on_negative


def _cohen_kappa(pairs: list[tuple[str, str]]) -> float:
    """Cohen's kappa for agreement between predicted and truth tier labels."""
    if not pairs:
        return 0.0
    n = len(pairs)
    categories = {c for pair in pairs for c in pair}
    observed = sum(1 for a, b in pairs if a == b) / n
    pred_counts = Counter(a for a, _ in pairs)
    truth_counts = Counter(b for _, b in pairs)
    expected = sum(
        (pred_counts.get(c, 0) / n) * (truth_counts.get(c, 0) / n) for c in categories
    )
    if expected >= 1.0:
        return 1.0
    return (observed - expected) / (1.0 - expected)
