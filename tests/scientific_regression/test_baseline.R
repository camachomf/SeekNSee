#!/usr/bin/env Rscript
# Minimal external baseline checks for original SeekNSee v0.2.5.
# Run at repo root: Rscript tests/scientific_regression/test_baseline.R
# This is a starter suite, NOT proof of the manuscript's 198 assertions.

stopifnot(file.exists("app/app.R"))
required_packages <- c("shiny", "dplyr", "stringr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R dependencies: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

# sys.source evaluates definitions/UI but does not print/launch the Shiny app.
# Use a dedicated environment to avoid polluting this test script.
seeknsee <- new.env(parent = globalenv())
sys.source("app/app.R", envir = seeknsee, keep.source = FALSE)

check <- function(name, expression) {
  tryCatch({
    stopifnot(isTRUE(expression))
    cat("PASS ", name, "\n", sep = "")
    invisible(TRUE)
  }, error = function(err) {
    stop("FAIL ", name, ": ", conditionMessage(err), call. = FALSE)
  })
}

check("software version", identical(seeknsee$SEEKNSEE_VERSION, "0.2.5"))
check("algorithm identifier", identical(seeknsee$ALGORITHM_VERSION, "2026-09-enzyme-aware-missed-cleavages-p1prime-1"))
check("nine enzyme choices", length(seeknsee$ENZYMES) == 9L)
check("twelve rules", length(seeknsee$RULE_NAMES) == 12L)

# Deterministic cleavage examples, independent from workbook/export code.
check("Trypsin disallows K|P", identical(seeknsee$digest_protein("AKPRK", "Trypsin")$peptide, c("AKPR", "K")))
check("Trypsin/P permits K|P", identical(seeknsee$digest_protein("AKPRK", "Trypsin/P")$peptide, c("AK", "PR", "K")))
check("AspN cleavage precedes D", identical(seeknsee$digest_protein("MADPEP", "AspN")$peptide, c("MA", "DPEP")))

missed <- seeknsee$digest_protein("AKRAK", "Trypsin", missed_cleavages = 2L)
check("0..2 missed cleavages enumerated", identical(sort(unique(missed$missed_cleavages)), 0:2))
check("missed-cleavage coordinate bounds", all(missed$start >= 1L & missed$end <= 5L & missed$length == missed$end - missed$start + 1L))
check("missed-cleavage 2 joins three fragments", any(missed$peptide == "AKRAK" & missed$missed_cleavages == 2L))

rules <- seeknsee$default_rules()
check("two default blocking rules", identical(names(rules)[vapply(rules, identical, logical(1), "block")], c("unique_proteome", "length")))
check("rule configuration accepted", isTRUE(seeknsee$validate_rule_configuration(rules)))

# Distinguish score from eligibility: terminal peptide has one apply penalty.
case_ok <- seeknsee$evaluate_peptide(
  pep = "AAAAAAAK", rules = rules,
  reference_proteome = seeknsee$empty_proteome(),
  target_id = "synthetic", target_source = "external_fasta",
  target_sequence = "AAAAAAAK", enzyme = "Trypsin"
)
check("8-aa terminal peptide retained", identical(case_ok$status, "retained"))
check("terminal penalty = -1", unname(case_ok$trace["terminal"]) == -1)
check("terminal case score = 11", identical(as.numeric(case_ok$score), 11))

case_block <- seeknsee$evaluate_peptide(
  pep = "AAAAAAK", rules = rules,
  reference_proteome = seeknsee$empty_proteome(),
  target_id = "synthetic", target_source = "external_fasta",
  target_sequence = "AAAAAAK", enzyme = "Trypsin"
)
check("7-aa peptide blocked by length", identical(case_block$status, "blocked") && "length" %in% case_block$blocked_by)
check("blocked rule does not directly lower score", unname(case_block$trace["length"]) == 0)

# Build and validate a tiny synthetic run without touching repository files.
ref_file <- tempfile(fileext = ".fasta")
writeLines(c(
  ">sp|P00001|SYNTHETIC_A GN=GENE1",
  "MPEPTIDEKAAAAAAAKGGGGR",
  ">sp|P00002|SYNTHETIC_B GN=GENE2",
  "MNQVAAAAAAAKTTTK"
), ref_file)
reference <- seeknsee$read_fasta(ref_file)
check("FASTA parsing retains two entries", nrow(reference) == 2L)
check("FASTA normalized UniProt accession", identical(reference$accession, c("P00001", "P00002")))

validation <- seeknsee$run_validation_case(
  reference_fasta = ref_file,
  target_source = "reference_proteome",
  target_id = reference$id[1],
  enzyme = "Trypsin", missed_cleavages = 1L,
  rules = rules
)
expected_sheets <- c("Run_summary", "Rule_configuration", "Target_proteins",
                     "Ranked_peptides", "Rule_trace", "Blocked_peptides", "Rule_definitions")
check("all seven workbook sheets exist", identical(names(validation$workbook), expected_sheets))
check("test FASTA yields candidate instances", nrow(validation$result$candidates) > 0L)
check("twelve trace rows per candidate", nrow(validation$result$rule_trace) == 12L * nrow(validation$result$candidates))
identities <- seeknsee$validate_internal_identities(validation)
check("all internal accounting identities", all(identities$passed))
unlink(ref_file)
cat("\nBaseline checks completed. These are starter assertions; further tests and manual review are required.\n")
