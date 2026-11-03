# Reproducibility and scientific regression — SeekNSee 0.2.5

**Status:** implementation-characterization tests passed on one Linux workstation. CI is proposed but has not yet executed on GitHub. This document does not claim empirical peptide detectability, cross-platform validation or recreation of all 198 assertions mentioned in the manuscript.

## Original baseline

- Source commit: `fff3fc7f2daf714e75dbbe4b6eec8f428e833b1e` (`main` at start of publication preparation).
- Baseline source SHA-256: `3cd007d023e6e0ddd59908a74accb156b709c8693d3dd9e7567694db0b907c27` (`app/app.R`, before comment-only documentation).
- Algorithm identifier: `2026-09-enzyme-aware-missed-cleavages-p1prime-1`.
- Software version: `0.2.5`.
- Frozen regression snapshot: `tests/scientific_regression/v0.2.5_full_outputs.rds`.
- Snapshot SHA-256: `689c7b71785f0a42782ec5801f29898799dfac66ccce76b406e1f60faf44c83b`.

The frozen snapshot represents four small, synthetic scenarios: reference Trypsin/1 missed cleavage, external Trypsin/P/2 missed cleavages under full-sequence exception policy, external AspN/0 missed cleavages under strict absence, and external Trypsin/0 missed cleavages with customized rules.

The four scenarios collectively produce 39 candidate instances and 468 candidate-rule records (12 per instance). The comparison includes every field and row of `candidates`, `rule_trace`, `target_summary` and all seven pre-export workbook tables, except the explicitly normalized runtime-specific `R version` and `Run date/time` entries in `Run_summary`.

## Observed local test environment

Recorded on 2026-11-03 from the author's Linux workstation:

| Component | Version |
| --- | --- |
| R | 4.5.0 |
| shiny | 1.14.0 |
| dplyr | 1.2.1 |
| stringr | 1.6.0 |
| writexl | 1.5.4 |

Additional transitive dependency versions are **not** yet locked. The initial GitHub Actions workflow installs current releases of these four packages, so its environment may differ from the local test environment. The next reproducibility milestone is reviewing and committing an `renv.lock` generated and restored in a clean environment.

## Run all currently available local checks

From the repository root:

```bash
Rscript -e 'invisible(parse(file="app/app.R")); cat("PASS syntax\n")'
Rscript tests/scientific_regression/test_baseline.R
Rscript tests/scientific_regression/test_edge_cases.R
Rscript tests/scientific_regression/compare_full_outputs.R compare
sha256sum tests/scientific_regression/v0.2.5_full_outputs.rds
```

**Do not run** `compare_full_outputs.R capture` on the refactored version. A captured baseline is meaningful only when generated from the original implementation and deliberately reviewed before freezing. Do not replace the snapshot solely to make a failing test pass.

For **comment-only** edits to `app/app.R`, an additional check can prove that executable R expressions are unchanged:

```bash
git show fff3fc7:app/app.R > /tmp/seeknsee_original_app.R
Rscript -e 'stopifnot(identical(parse("/tmp/seeknsee_original_app.R", keep.source=FALSE), parse("app/app.R", keep.source=FALSE))); cat("PASS executable expressions unchanged\n")'
```

This expression-level comparison is intended only for comment-only changes, **not** for later functional modularization. After modularization, equivalent outputs and independent targeted tests will be the acceptance criteria.

## Verified, not verified and follow-ups

- **Local tests verified:** 23 starter assertions, 49 edge-case assertions, and 4 full-output regression scenarios (passed on the local R installation).
- **Not independently verified:** 198 assertions and benchmark totals reported by the paper; experimental peptide detectability; large-proteome memory/time; end-to-end Shiny interactions; Windows/macOS installers.
- **Export limitations:** Excel preflight checks cell-character limits and may abbreviate identifier lists with omission counts. Excel's maximum worksheet row count (1,048,576) is not explicitly validated by the current application, and no lossless multi-file export exists yet for very large rule traces.
- **CI:** `.github/workflows/scientific-validation.yml` is an initial Linux test gate, not a cross-platform certification. Review its first remote execution logs before describing it as passing.

## Release gate

A public scientific release should be tagged only when the maintainer approves the code, dependency lockfile, license, citation metadata, installation instructions, CI, examples and correspondence with the version used in the manuscript. Keep `README_FIRST.md` and other temporary phase instructions out of release commits.
