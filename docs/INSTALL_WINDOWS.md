# Install SeekNSee from source — Windows x64

**Status:** Draft source-installation instructions; they still require hands-on testing on Windows. The project README reports a previously built, manually accepted Windows 10/11 x64 distribution with a bundled R runtime and package library. Its validated release archive and build scripts are not yet available in this repository; the instructions below describe a **separate, source-based** approach.

## Prerequisites

- Windows 10/11 x64, a browser, and Git (or the GitHub ZIP once the repository is public).
- R 4.5.0 to mirror the reported local scientific test environment; use the official [CRAN Windows downloads](https://cran.r-project.org/bin/windows/base/).
- Rtools matching the installed R version may be required if packages must be compiled from source; check the official [Rtools documentation](https://cran.r-project.org/bin/windows/Rtools/).

Open PowerShell, clone the repository (once accessible), and enter its root:

```powershell
git clone https://github.com/camachomf/SeekNSee.git
Set-Location SeekNSee
```

Install required packages from R or from a PowerShell session where `Rscript.exe` is on `PATH`:

```powershell
Rscript -e 'install.packages(c("shiny", "dplyr", "stringr", "writexl"), repos="https://cloud.r-project.org")'
```

Start the app from the repository root:

```powershell
Rscript -e 'shiny::runApp("app", launch.browser=TRUE)'
```

Run local scientific checks:

```powershell
Rscript tests/scientific_regression/test_baseline.R
Rscript tests/scientific_regression/test_edge_cases.R
Rscript tests/scientific_regression/compare_full_outputs.R compare
```

**Caution:** These commands have not yet been executed on Windows as part of the present repository-preparation work. The README's claims for the bundled Windows build must be supported by its original build artifacts and acceptance records before the corresponding release is published.

The future version-pinned installation procedure is documented in [ENVIRONMENT.md](ENVIRONMENT.md).
