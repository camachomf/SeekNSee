#!/usr/bin/env Rscript
# SeekNSee v0.2.5: complete scientific-output snapshot and comparison.
# Usage from repository root:
#   Rscript tests/scientific_regression/compare_full_outputs.R capture
#   Rscript tests/scientific_regression/compare_full_outputs.R compare
# This tool deliberately does not change app/app.R or write Excel files.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("capture", "compare")) {
  stop("Usage: Rscript tests/scientific_regression/compare_full_outputs.R {capture|compare}", call. = FALSE)
}
mode <- args[[1L]]
if (!file.exists("app/app.R")) {
  stop("Run from the SeekNSee repository root (app/app.R not found).", call. = FALSE)
}

required_packages <- c("shiny", "dplyr", "stringr", "writexl")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)

app <- new.env(parent = globalenv())
sys.source("app/app.R", envir = app, keep.source = FALSE)

# Stable *file basenames* are important: they appear in the Run_summary.
fixture_dir <- file.path(tempdir(), "seeknsee_complete_regression")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
reference_path <- file.path(fixture_dir, "reference.fasta")
external_path <- file.path(fixture_dir, "external_targets.fasta")
writeLines(c(
  ">sp|P00001|SYNTH_A GN=GENE1", "MPEPTIDEKAAAAAAAKGGGGR",
  ">sp|P00002|SYNTH_B GN=GENE2", "MNQVAAAAAAAKTTTK",
  ">sp|P00003|SYNTH_C GN=GENE3", "AKPRAKPEPTIDERKSTVVVVVVVK",
  ">sp|P00001-2|SYNTH_A_ISOFORM GN=GENE1", "MPEPTIDEKAAAAAAAKGGGGR"
), reference_path, useBytes = TRUE)
writeLines(c(
  ">EXTERNAL_A", "MPEPTIDEKAAAAAAAKGGGGR",
  ">EXTERNAL_B", "MPEPTIDERAKPEPTIDERKSTVVVVVVVK",
  ">EXTERNAL_C", "MNQWCMNQGDKTTTK"
), external_path, useBytes = TRUE)
reference <- app$read_fasta(reference_path)

custom_rules <- app$default_rules()
custom_rules$unique_proteome <- "apply"
custom_rules$length <- "apply"
custom_rules$isobaric <- "block"
custom_rules$charge <- "ignore"
app$validate_rule_configuration(custom_rules)

# Each run exercises a different input policy, enzyme or rule configuration.
# All cases use the same small, deterministic FASTA input bytes.
scenarios <- list(
  reference_trypsin_m1 = list(
    target_source = "reference_proteome", target_id = reference$id[[1]],
    enzyme = "Trypsin", missed_cleavages = 1L, rules = app$default_rules()
  ),
  external_trypsinp_m2_allow = list(
    target_source = "external_fasta", target_fasta = external_path,
    target_scope = "all", reference_policy = "allow_exact_target_match",
    enzyme = "Trypsin/P", missed_cleavages = 2L, rules = app$default_rules()
  ),
  external_aspn_m0_strict = list(
    target_source = "external_fasta", target_fasta = external_path,
    target_scope = "all", reference_policy = "strict_absence",
    enzyme = "AspN", missed_cleavages = 0L, rules = app$default_rules()
  ),
  external_trypsin_m0_custom = list(
    target_source = "external_fasta", target_fasta = external_path,
    target_scope = "all", reference_policy = "strict_absence",
    enzyme = "Trypsin", missed_cleavages = 0L, rules = custom_rules
  )
)

normalize_workbook <- function(workbook) {
  stopifnot(is.list(workbook), is.data.frame(workbook$Run_summary))
  # Remove ONLY run-specific presentation metadata, not scientific results.
  volatile_fields <- c("Run date/time", "R version")
  workbook$Run_summary <- workbook$Run_summary[
    !workbook$Run_summary$field %in% volatile_fields, , drop = FALSE
  ]
  workbook
}

run_scenario <- function(scenario) {
  validation <- do.call(
    app$run_validation_case,
    c(list(reference_fasta = reference_path), scenario)
  )
  identities <- app$validate_internal_identities(validation)
  if (anyNA(identities$passed) || !all(identities$passed)) {
    stop("An internal accounting identity failed.", call. = FALSE)
  }
  stopifnot(identical(names(validation$workbook), c(
    "Run_summary", "Rule_configuration", "Target_proteins",
    "Ranked_peptides", "Rule_trace", "Blocked_peptides", "Rule_definitions"
  )))
  stopifnot(nrow(validation$result$rule_trace) ==
              nrow(validation$result$candidates) * length(app$RULE_NAMES))
  list(
    # All candidate columns, all candidate rows, and original row order.
    candidates = validation$result$candidates,
    # All 12 rule rows per candidate, with original order and all fields.
    rule_trace = validation$result$rule_trace,
    target_summary = validation$result$target_summary,
    # All seven complete workbook sheets; no Excel truncation or row limits.
    workbook = normalize_workbook(validation$workbook),
    accounting_identities = identities,
    # Stable run inputs for auditing; no elapsed seconds or timestamp.
    run_configuration = validation$context[c(
      "reference_filename", "target_filename", "target_source",
      "reference_policy", "enzyme", "missed_cleavages", "rules"
    )]
  )
}

run_all <- function() {
  out <- lapply(names(scenarios), function(name) {
    message("Running scenario: ", name)
    run_scenario(scenarios[[name]])
  })
  names(out) <- names(scenarios)
  out
}

baseline_file <- "tests/scientific_regression/v0.2.5_full_outputs.rds"
output <- list(
  software_version = app$SEEKNSEE_VERSION,
  algorithm_version = app$ALGORITHM_VERSION,
  fixture_reference = readLines(reference_path, warn = FALSE),
  fixture_external = readLines(external_path, warn = FALSE),
  scenarios = run_all()
)

if (identical(mode, "capture")) {
  if (file.exists(baseline_file)) {
    stop("Baseline already exists. Refusing to overwrite: ", baseline_file,
         "\nCompare against it, or archive it before creating a deliberate new baseline.",
         call. = FALSE)
  }
  saveRDS(output, baseline_file, version = 3, compress = "xz")
  message("CAPTURED full baseline: ", baseline_file)
  message("Size: ", file.info(baseline_file)$size, " bytes")
  message("No application source files were changed.")
  quit(save = "no", status = 0L)
}

if (!file.exists(baseline_file)) {
  stop("No baseline found. First run with 'capture' on the unmodified v0.2.5.", call. = FALSE)
}
expected <- readRDS(baseline_file)

# Compare all content exactly, including rows, columns, types, attributes and order.
# Do not discard candidates or summarize trace records.
failures <- character()
compare_section <- function(label, before, after) {
  if (identical(before, after)) {
    cat("PASS ", label, "\n", sep = "")
  } else {
    explanation <- all.equal(before, after, tolerance = 0,
                             check.attributes = TRUE)
    explanation <- if (isTRUE(explanation)) "objects differ" else paste(head(explanation, 6), collapse = " | ")
    message("FAIL ", label, ": ", explanation)
    failures <<- c(failures, label)
  }
}
compare_section("software_version", expected$software_version, output$software_version)
compare_section("algorithm_version", expected$algorithm_version, output$algorithm_version)
compare_section("fixture_reference", expected$fixture_reference, output$fixture_reference)
compare_section("fixture_external", expected$fixture_external, output$fixture_external)
compare_section("scenario_names", names(expected$scenarios), names(output$scenarios))

for (scenario_name in intersect(names(expected$scenarios), names(output$scenarios))) {
  old <- expected$scenarios[[scenario_name]]
  new <- output$scenarios[[scenario_name]]
  for (field in c("candidates", "rule_trace", "target_summary",
                  "accounting_identities", "run_configuration")) {
    compare_section(paste(scenario_name, field, sep = "/"), old[[field]], new[[field]])
  }
  compare_section(paste0(scenario_name, "/workbook_sheets"),
                  names(old$workbook), names(new$workbook))
  for (sheet_name in intersect(names(old$workbook), names(new$workbook))) {
    compare_section(paste0(scenario_name, "/workbook/", sheet_name),
                    old$workbook[[sheet_name]], new$workbook[[sheet_name]])
  }
  cat(sprintf("  Counts for %s: %d candidates, %d rule rows\n",
              scenario_name, nrow(new$candidates), nrow(new$rule_trace)))
}
if (length(failures)) {
  stop(sprintf("COMPLETE REGRESSION FAILED: %d differences detected", length(failures)),
       call. = FALSE)
}
cat("PASS COMPLETE REGRESSION: every captured scientific table matches exactly.\n")
