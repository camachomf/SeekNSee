# Environment reproducibility — SeekNSee v0.2.5

## Recorded local environment

The following versions were reported from the local Linux regression environment:

| Component | Version |
| --- | --- |
| R | 4.5.0 |
| shiny | 1.14.0 |
| dplyr | 1.2.1 |
| stringr | 1.6.0 |
| writexl | 1.5.4 |

These values alone do not pin transitive dependencies. A validated `renv.lock` must be generated on a machine where the project is functional and reviewed before publication.

## Step 1: Record dependencies without changing the active project library

From the repository root:

```bash
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org"); print(renv::dependencies(path=".", progress=FALSE))'
```

Inspect discovered packages; dynamically loaded dependencies or platform-specific packaging dependencies may require separate handling.

## Step 2: Generate a *candidate* lockfile

`renv::snapshot()` can snapshot a project that has **not yet been activated with renv**, using the available library paths. This initially avoids changing `~/.Rprofile`, the project `.Rprofile`, or R's active library setup. It **does write `renv.lock`**. This is a controlled preparation step, not a verified release artifact.

```bash
Rscript -e 'renv::snapshot(project=getwd(), packages=c("shiny","dplyr","stringr","writexl"), prompt=FALSE)'
```

`packages=` includes recursive package dependencies; inspect the resulting lockfile for provenance, repositories, and versions. Do not manually type a synthetic lockfile from the version table above. The direct packages must match the validated environment; inspect transitive packages and the actual `R` version in the generated lockfile.

Confirm that the application remains unchanged and rerun the tests:

```bash
git diff --check
git diff -- app/app.R
Rscript tests/scientific_regression/test_baseline.R
Rscript tests/scientific_regression/test_edge_cases.R
Rscript tests/scientific_regression/compare_full_outputs.R compare
```

## Step 3: Restore in a disposable clone before approving

Do **not** test `renv::restore()` against a shared/global R library. Create a separate clean clone or isolated project directory, initialize project-local renv infrastructure there, and restore the proposed lockfile. The restoration may download and compile packages and therefore requires network access and platform toolchains.

Only after the restore and the three test commands pass in that isolated environment should the maintainers approve adding `renv.lock` and project activation files to the release branch.

## Step 4: CI integration

Once `renv.lock` and the project activation strategy are reviewed, update `.github/workflows/scientific-validation.yml` to restore locked dependencies using `r-lib/actions/setup-renv@v2` rather than installing unconstrained current packages. Test the workflow on the branch **before** marking the pull request ready.

## Caveats

- `renv` manages R package versions, not the R executable, Linux system libraries, or browser behavior.
- The Linux environment and bundled Windows/Ubuntu installers require distinct acceptance evidence.
- The four available regression scenarios cover 39 candidate instances and 468 rule rows, not every input combination or proteome-scale performance.
- SeekNSee's heuristic scores are not calibrated peptide-detectability probabilities.
