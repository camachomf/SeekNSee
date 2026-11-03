# Troubleshooting SeekNSee v0.2.5

These checks apply primarily to execution of the source-based R/Shiny application.

## `app.R` not found

The application lives in `app/app.R`, **not** `./app.R`. From the repository root run:

```bash
Rscript -e 'shiny::runApp("app", launch.browser=TRUE)'
```

## R package is missing

Confirm R and direct package versions:

```bash
Rscript -e 'cat(R.version.string, "\n"); for (p in c("shiny","dplyr","stringr","writexl")) cat(p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "MISSING", "\n")'
```

If a project lockfile has been approved, restore that environment using `renv` rather than installing arbitrary current package versions. Otherwise install missing packages from CRAN.

## FASTA upload fails

The parser expects `>` headers followed by nonempty amino-acid sequences. It removes whitespace, normalizes case, removes trailing `*`, and rejects internal stop symbols or unsupported sequence characters. Validate input provenance separately; parsing success does not prove a sequence is biologically correct.

## Peptide not retained

Check (1) digestion enzyme and missed-cleavage setting, (2) candidate length (4–25 amino acids before rule evaluation), (3) configured blocking rules, (4) reference-target provenance and specificity policy, and (5) rule-by-rule trace in the Excel workbook. A blocked peptide may still have a numeric score.

## The application exits after closing the browser

The existing Shiny server tracks active sessions and calls `stopApp()` when the final session ends. This behavior is suited to the present local application but may need revision before deployment as a hosted multiuser service.

## Excel export fails or shows abbreviated reference IDs

Excel limits a cell to 32,767 characters. SeekNSee adapts certain semicolon-separated identifier lists and reports counts of omitted identifiers; this does **not** mean the omitted identifiers were excluded from scoring. Oversized unsupported text fields can raise a preflight error.

**Important:** Excel also limits each worksheet to 1,048,576 rows. The v0.2.5 implementation does not explicitly preflight that row limit; an XLSX export with millions of rule trace rows must not be assumed to be complete or supported without additional testing.

## Scientific comparison fails

Run the functional checks first, then the comparison:

```bash
Rscript tests/scientific_regression/test_baseline.R
Rscript tests/scientific_regression/test_edge_cases.R
Rscript tests/scientific_regression/compare_full_outputs.R compare
```

A failed regression should be investigated against the original v0.2.5 baseline. **Never overwrite the snapshot using `capture` simply to force the comparison to pass.**
