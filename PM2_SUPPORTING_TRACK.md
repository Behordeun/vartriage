# PM2 Supporting default: engine decision track

Status: resolved (ship Supporting, correct release note)
Owner: Muhammad
Opened: 2026-10-08
Scope: engine and release only.

## The question

v1.1.0 changes the shipped default strength of PM2 (absence from population databases) from Moderate (2 points) to Supporting (1 point), aligning with ClinGen SVI 2020. The merge is already on main as bfbb11c (PR #80). The release is not yet on PyPI: the latest published tag is v1.0.2, and pyproject.toml is bumped to 1.1.0 on the merged line but not tagged or uploaded.

So this is a catch-before-users decision, not a post-ship cleanup. The one open question: is Supporting the right shipped default to publish, and does the release note describe its effect correctly?

Status update 2026-10-08: the reconciliation below is DONE and it reversed an earlier wrong conclusion. The changelog's roughly-8-percent downgrade figure is correct. An earlier measurement that reported no detection loss was run on a narrower cohort with a different tag join and did not describe the shipped engine on the full eRepo set.

## What is settled

- Supporting is the guideline-concordant weight. ClinGen SVI 2020 recommends PM2 at Supporting because absence from population databases is weak standalone evidence: most rare benign and uncertain variants are also absent. This direction is not in dispute.
- The change ships with an auditable escape hatch. score_points, classify_by_points, and combine_evidence take an optional strength_overrides mapping, so a lab with a validated gene-specific VCEP specification can apply a different PM2 weight per run with the applied weight traceable to either the documented default map or the caller's specification. This is the right design: the default is calibrated, the deviation is explicit.
- The default lives in one place: EvidenceTag.PM2 -> EvidenceStrength.SUPPORTING in vartriage/models/variant.py, consumed through EVIDENCE_STRENGTH_MAP.

## Reconciliation (resolved 2026-10-08)

Two measurements of the same change disagreed, and the disagreement was a cohort artifact. Settled by classifying every variant's committed evidence tags through the real engine's classify_by_points at both PM2 strengths, on the full point-system classification set (data/erepo/erepo_classifications_pointsystem.json), and joining to the ClinGen expert truth table (data/erepo/clingen_erepo.tsv). Result written to vartriage-acmg-concordance/results/pm2_reconciliation.json.

On 20,160 evaluable variants (18,151 carry PM2), PM2 Moderate -> Supporting moves:

- 1,384 Likely Pathogenic -> VUS
- 430 VUS -> Likely Benign
- 1,814 total, 9.00 percent of evaluable variants (9.99 percent of PM2 carriers), every move one tier toward uncertainty. Pathogenic and Benign tiers unchanged.

Joined to expert truth (12,203 expert Pathogenic/Likely Pathogenic variants matched):

- pathogenic sensitivity 0.5785 at Moderate, 0.5019 at Supporting
- delta -0.0766 (935 fewer expert-pathogenic variants detected at Supporting)

So the merged commit message and the 1.1.0 changelog are correct: the change is a real one-tier-conservative shift of roughly 8 to 9 percent of variants, and it does lose pathogenic detection. The earlier pm2_clean_engine_measure.json result (0 detection flips, only 390 Pathogenic -> Likely Pathogenic) was measured on a narrower 21,506-record superset with a different tag join; it does not describe the shipped engine on the full eRepo set and should not be cited.

Changelog wording note: the merged commit body says 8.43 percent and the published changelog entry says roughly 8 percent. The clean full-set figure is 9.00 percent of evaluable variants. The changelog's "roughly 8%" is defensible as a round-down but is on the low side; "roughly 9 percent (1,814 of 20,160 variants: 1,384 LP to VUS, 430 VUS to LB)" is the precise statement and is what the changelog is corrected to.

## Decision

Ship Supporting as the default (option A), with the release note corrected to the exact transition counts and the sensitivity delta stated plainly. Rationale:

- Supporting is the ClinGen SVI 2020 calibrated weight. The sensitivity drop is the intended effect of removing over-weighted rarity evidence, not a regression. A variant whose Likely Pathogenic call rested partly on PM2 at Moderate was over-confident under the guideline, and VUS is the honest tier for it.
- Every move is one step toward uncertainty. No variant becomes a more confident pathogenic call, and the confident Pathogenic and Benign tiers are untouched. The clinically dangerous direction (a true pathogenic call flipping to benign) does not occur.
- The strength_overrides primitive lets a lab with a validated gene-specific VCEP specification restore a higher PM2 weight per run, traceably, where that is justified.

Release-note honesty requirement before tagging: state the sensitivity cost (0.58 -> 0.50 on this cohort, -7.7 pp) in the user-facing notes, not only the reclassification percentage. Users running the tool for triage must know the default got more conservative in a way that lowers recall.

## Next action

Correct the 1.1.0 changelog entry to the precise reconciled numbers plus the sensitivity delta (done alongside this update), then tag and publish v1.1.0 when ready.
