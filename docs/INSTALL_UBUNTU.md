# Install SeekNSee from source — Ubuntu x86_64

**Status:** Source-code instructions for **SeekNSee v0.2.5**, prepared for local verification. These instructions are distinct from the already reported, bundled Ubuntu distribution. The README identifies Ubuntu 24.04 LTS x86_64 as a validated packaging target; this repository does not yet include its build scripts or release archive.

## Prerequisites

- Ubuntu 24.04 LTS x86_64 (the reported Linux packaging target).
- **R 4.5.0** to reproduce the environment used for the local regression runs; other R versions are not established as equivalent by those runs.
- Git, an internet connection for package installation, and a browser.
- R build tools/system libraries may be needed when compiled dependencies lack available binaries.

Install Git and development prerequisites (these commands do **not** guarantee installation of R 4.5.0):

```bash
sudo apt-get update
sudo apt-get install -y git build-essential libcurl4-openssl-dev libssl-dev libxml2-dev
```

Install an appropriate R release using the official [CRAN Ubuntu instructions](https://cran.r-project.org/bin/linux/ubuntu/) and confirm:

```bash
Rscript -e 'cat(R.version.string, "\n")'
```

## Clone and install R dependencies

```bash
git clone https://github.com/camachomf/SeekNSee.git
cd SeekNSee
Rscript -e 'install.packages(c("shiny", "dplyr", "stringr", "writexl"), repos="https://cloud.r-project.org")'
```

The clone command only works after the project is made public, or for an account authorized to access it. Until the project lockfile has been generated and restore-tested, package installation through CRAN is **not fully reproducible**; versions may change. After `renv.lock` is approved, prefer the procedure in [ENVIRONMENT.md](ENVIRONMENT.md).

## Start the application

From the repository root:

```bash
Rscript -e 'shiny::runApp("app", launch.browser=TRUE)'
```

The Shiny application is located in `app/app.R`, not in the repository root. Use `Ctrl+C` to stop a terminal-launched instance; the current local-session implementation can also stop when its last browser session disconnects.

## Verify the source implementation

```bash
Rscript -e 'invisible(parse("app/app.R")); cat("PASS syntax\n")'
Rscript tests/scientific_regression/test_baseline.R
Rscript tests/scientific_regression/test_edge_cases.R
Rscript tests/scientific_regression/compare_full_outputs.R compare
```

The regression comparison requires the reference snapshot `tests/scientific_regression/v0.2.5_full_outputs.rds`, which should be tracked and kept unchanged. **Do not run `capture` during validation.** These tests characterize sequence-based software behavior, not experimental peptide detectability.

## Supported input and outputs

The application reads reference and target protein FASTA/FAA files and exports a seven-sheet XLSX report. Large datasets require particular care: Excel has a limit of 1,048,576 rows per worksheet, and the current implementation does not explicitly check that row bound. See [TROUBLESHOOTING.md](TROUBLESHOOTING.md).
