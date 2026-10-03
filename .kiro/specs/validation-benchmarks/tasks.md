# Validation and Benchmarks (v1.0.0) - Implementation Tasks

## Truth set construction

- [x] 1. Script: download ClinVar VCF from NCBI FTP (GRCh38, latest release) — `scripts/build_truth_set.py`; pinned release `clinvar_20260928` (GRCh38, 4,555,206 records)
- [x] 2. Filter to review_status in ("reviewed by expert panel", "practice guideline")
- [x] 3. Exclude records with conflicting interpretations
- [x] 4. Split into positive (P/LP) and negative (B/LB) sets — 12,545 positive / 5,621 negative
- [x] 5. Normalize all variants using VariantNormalizer (left-align against GRCh38 FASTA)
- [x] 6. Output: truth_set.vcf + truth_labels.tsv with (chrom, pos, ref, alt, gene, clinvar_class)
- [x] 7. Document: pin ClinVar release date, record exact filters applied — `truth_set_provenance.json`

## Benchmark runner

- [x] 8. Implement BenchmarkRunner: run pipeline on truth VCF, load output, match to labels — `vartriage/validation/benchmark.py`
- [x] 9. Implement coordinate matching (normalized chrom:pos:ref:alt as key)
- [x] 10. Handle unmatched variants (truth variants not in pipeline output and vice versa)
- [x] 11. Implement MetricWithCI dataclass for all reported metrics
- [x] 12. Implement wilson_ci() for binomial confidence intervals
- [x] 13. Compute sensitivity: TP / (TP + FN) where positive = P/LP in truth
- [x] 14. Compute specificity: TN / (TN + FP) where negative = B/LB in truth
- [x] 15. Compute PPV, NPV
- [x] 16. Compute Cohen's kappa for 5-class agreement
- [x] 17. Generate 5x5 confusion matrix (P, LP, VUS, LB, B predicted vs truth)

> Benchmark runner proven on the chr22 harness (43/43 matched). The
> **full-genome release figures are not yet produced** — the run is blocked on
> staging genome-wide reference data (REVEL is a dehydrated OneDrive
> placeholder; genome-wide gnomAD frequency and CADD are not local). See
> `stage3-full-run-recipe.md`. No release-grade number exists until that run
> completes, and the v1.0.0 release claim is gated on it per
> `accuracy-decision-rule.md`.

## Per-criterion concordance

- [ ] 18. For each evidence tag: count fires on truth-positive vs truth-negative variants
- [ ] 19. Compute per-criterion precision (fires on positive / all fires)
- [ ] 20. Identify criteria with high false-fire rate on benign variants
- [ ] 21. Report as a table in the validation document

## Consequence accuracy comparison

- [ ] 22. Implement ConsequenceComparison: run same VCF through local and API mode
- [ ] 23. Pair results by coordinate, compare consequence fields
- [ ] 24. Report: overall concordance %, discordant cases by type
- [ ] 25. Identify systematic patterns (e.g., local calls synonymous that VEP calls missense)

## Tool comparison

- [ ] 26. Install InterVar (or prepare a pre-run dataset from their web tool)
- [ ] 27. Run truth set through InterVar, parse classification output
- [ ] 28. Compare: vartriage vs InterVar vs ClinVar expert (3-way confusion matrix)
- [ ] 29. Identify unique correct/incorrect calls per tool
- [ ] 30. Document differences with explanations (different data sources, different criteria versions)

## Validation report

- [ ] 31. Implement report_generator: produce docs/validation_report.md from BenchmarkReport
- [ ] 32. Include: truth set description, pipeline config, all metrics with CIs, confusion matrix
- [ ] 33. Include: per-criterion table, known limitations, comparison summary
- [ ] 34. Generate publication-ready figures (matplotlib): ROC-like plot, bar chart per criterion

## REVEL coverage analysis

- [ ] 35. Quantify: for the truth set, % of missense variants with REVEL scores
- [ ] 36. Measure PP3 fire rate with vs without REVEL data
- [ ] 37. Measure sensitivity impact when REVEL is removed
- [ ] 38. Document recommendation in docs/configuration.md under PP3 guidance

## Release prep

- [ ] 39. Run full test suite (target: 1000+ tests after all releases)
- [ ] 40. mypy strict 3.10/3.11/3.12 clean
- [ ] 41. All documentation current
- [ ] 42. Update README with accuracy metrics summary from validation report
- [ ] 43. Changelog v1.0.0 entry
- [ ] 44. Commit, push `feature/validation-benchmarks` branch, raise PR
- [ ] 45. After merge: tag v1.0.0, publish to PyPI, create GitHub Release
