# Pre-registered accuracy decision rule (v1.0.0)

Written **before** the v1.0.0 benchmark is run, and committed, so the
interpretation of the result is fixed while the result is still unknown.
This removes the pull to move the goalposts once numbers are in hand. The
figures produced by the benchmark are reported exactly as measured; this
document decides only how the release is *positioned* given those figures.

Status: pre-registered. Do not edit the thresholds below after the
benchmark has been run. If a threshold turns out to be wrong, that is
itself a finding to report, not a number to quietly change.

## Scope of the claim

- The release figures are the **full-genome** run against the pinned
  ClinVar `clinvar_20260928` expert-curated truth set (18,166 variants:
  12,545 P/LP positive, 5,621 B/LB negative), review tiers
  reviewed-by-expert-panel and practice-guideline, conflicting excluded.
- A chromosome-scoped (chr22) run is a **harness check only**. Its numbers
  validate that the measurement works; they are never reported as the
  library's accuracy and never trigger the decision below.
- Metrics are computed with Wilson confidence intervals. A metric whose
  95% CI is wider than ±5 points is reported as **underpowered** for its
  stratum and does not, on its own, move the decision.

## Primary metric and the clinical bar

The primary metric is **PPV on the pathogenic call** (of the variants the
tool calls Pathogenic or Likely Pathogenic, the fraction that are truly
P/LP in the expert truth set). US-1 names the clinical adoption bar as
**95% PPV**. Pathogenic sensitivity (recall of true P/LP) is the primary
secondary metric. Benign sensitivity and specificity are reported but are
known-weak (BP1/BP3/BP6 unimplemented) and do not gate the clinical-vs-
research decision.

A variant classified VUS counts as a non-call: it is neither a true nor a
false pathogenic call. PPV and pathogenic sensitivity are computed on the
P/LP-call cell of the confusion matrix, with VUS handled explicitly and
its rate reported.

## The decision ladder (committed in advance)

The release is positioned by the **lower bound of the 95% Wilson CI** on
PPV, not the point estimate, so a small or lucky sample cannot buy a
stronger claim than the data supports.

- **Clinical-grade** — PPV lower-CI-bound ≥ 0.95 AND pathogenic
  sensitivity point estimate ≥ 0.70. The release may describe itself as
  meeting the stated clinical PPV bar, with the measured sensitivity
  stated plainly beside it. Benign-side weakness is still disclosed.

- **Research-grade (expected)** — PPV lower-CI-bound ≥ 0.80 but below the
  clinical bar, OR sensitivity below 0.70. This is the honest prior given
  the last v0.17.x run (86% PPV, ~70% pathogenic sensitivity, 6% benign
  sensitivity). The release is positioned as a **research / computational-
  triage tool, not a diagnostic one**: the README and validation report
  lead with "computational triage, not a diagnostic classifier," state
  every metric with its CI, and keep the Development Status at Beta rather
  than claiming Production/Stable on the strength of the numbers.

- **Below research-grade** — PPV lower-CI-bound < 0.80. The numbers are
  reported in full and the release is **held**: v1.0.0 does not ship on
  these figures. We investigate the mechanism (a systematic criterion
  firing falsely on benign variants is the first suspect, per the
  per-criterion concordance in FR-3) before any release decision. We do
  not reach the bar by retuning thresholds against the truth set — that
  would be fitting to the test and is forbidden.

## Hard constraints (non-negotiable)

- The combining-path thresholds (PP3/BP4 REVEL cut-offs, the SVI point
  values, the PM2/BA1/BS1 frequency gates) are **not** tuned to improve
  the benchmark result. They are fixed by their published calibrations
  (Pejaver et al. 2022; Tavtigian et al. 2018/2020). Fitting them to this
  truth set would make the validation measure nothing.
- No validation number is fabricated, rounded favourably, or presented
  without its confidence interval. An unavailable comparison (e.g. an
  InterVar run that cannot be completed) is recorded as "not performed,"
  never estimated.
- The 1.0.0 Development-Status bump to "Production/Stable" is gated on the
  clinical-grade branch above. The research-grade branch ships as 1.0.0
  but keeps Beta status, because a 1.0 version number is a stability
  statement about the API, not a clinical claim about the classifications.

## What the report states regardless of branch

- Truth-set source, pinned release date, size, and composition.
- Every metric with its 95% Wilson CI and the N behind it.
- The per-criterion false-fire table (FR-3): which criteria fire on truly
  benign variants, since that is the most actionable weakness signal.
- The REVEL-coverage impact (FR-7): how much of the missense sensitivity
  depends on REVEL being present.
- A plain-language limitations section naming the benign-side gap.
