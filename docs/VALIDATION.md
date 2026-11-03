# SeekNSee verification and scientific regression (v0.2.5)

> **Documentation draft; not a test-execution report.** Existing repository inventory contained only README placeholders in `tests/packaging/` and `tests/scientific_regression/`. Previously reported benchmark totals and assertion counts have **not been independently reproduced** by this draft.

## Objectives

The first release gate is to freeze a reproducible behavioral baseline **before** modularization or optimization. Tests should distinguish four kinds of evidence:

1. **Unit/edge-case tests:** FASTA validation, cleavage rules, coordinates, missed-cleavage enumeration, rule policy, scoring and specificity.
2. **Scientific regression:** fixed small FASTA fixtures and expected candidate/trace outputs generated from the original v0.2.5 source, checked against a refactored implementation.
3. **Integration/UI/export:** Shiny startup, upload and settings, workbook generation, file integrity, and Excel limits.
4. **Manuscript reproduction:** complete reference FASTA and experimental datasets, exact processing scripts, configurations, benchmark results, and platform logs.

These levels are **not interchangeable**: passing unit tests cannot prove manuscript-scale counts or predictive-performance comparisons.

## Candidate accounting identities

The existing `validate_internal_identities()` checks these invariants:

- `candidate_count = retained_count + blocked_count`.
- `Rule_trace` has `candidate_count × 12` rows.
- Every candidate instance has exactly 12 rule rows.
- `peptide_score = 12 + sum(rule_score)`.

Run these on small controlled data and, when feasible, on representative high-volume batches. Source functions also expose `run_validation_case()`, `evaluate_test_peptide()` and `compare_digestion()`.

## Minimal test execution

The companion starter test is `tests/scientific_regression/test_baseline.R`, intended for execution **from the repository root**:

```bash
Rscript tests/scientific_regression/test_baseline.R
```

It uses the application's existing dependencies (`shiny`, `dplyr`, `stringr`, `writexl`) plus base-R assertions, without refactoring any production code. The script loads `app/app.R` into a private environment without starting an interactive Shiny server and checks version identifiers, enzyme cleavage examples, missed cleavages, rule-state behavior, FASTA parsing and candidate-trace identities. A local R runtime and dependency installation are required. **The draft was not executed in its preparation environment.**

## Required release additions

- Turn minimal tests into a maintained suite (`testthat` is appropriate), including **all nine enzyme configurations**, 0–2 missed cleavages, peptide length boundaries, all 12 rules and ignore/apply/block states, protein-specificity policies, duplicate/isoform IDs, and Excel long-cell behavior.
- Snapshot v0.2.5 reference outputs from a fixed fixture, retain expected results under version control, and compare scientific columns without normalization that could conceal changes.
- Record **SHA-256 input checksums**, R/package versions, algorithm identifier, Git commit SHA, test commands, outcomes and dates.
- Confirm the previously reported **198 automated assertions** from the original scripts/logs or amend the manuscript to reflect what can actually be reproduced.
- Check the actual Excel row limit (**1,048,576 rows per sheet**) before export. Current v0.2.5 validates cell length and abbreviates long identifier lists, but does not check sheet row counts. Decide on explicit failure versus a lossless alternative export for scale.
- Verify Windows and Ubuntu packaging with documented clean-machine test runs; macOS packaging should be described as validated only after it is tested.
- Recover and pin all benchmark FASTAs, input filtering, datasets, score-label definitions, PeptideRanger version/settings, metrics scripts, plots and outputs before claiming numerical reproduction of the manuscript.

## Pass/fail record template

| Field | Value |
| --- | --- |
| Git commit | _record SHA_ |
| SeekNSee version | 0.2.5 |
| Algorithm version | 2026-09-enzyme-aware-missed-cleavages-p1prime-1 |
| Input FASTA SHA-256 | _record digest_ |
| R version / OS | _record_ |
| Dependency lockfile | _record hash_ |
| Test command | _record_ |
| Tests / assertions executed | _record actual counts_ |
| Failures / skipped | _record_ |
| Maintainer approval | _record_ |

## Acceptance criterion

Refactoring must preserve candidate instances, coordinates, rule applicability, scores, retention/blocking, and trace order **for the same data, configuration and algorithm identifier**, except for expressly reviewed, versioned changes. Do not merge or publish claims of equivalence before runtime checks complete.
