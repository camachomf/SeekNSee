#!/usr/bin/env Rscript
# SeekNSee v0.2.5 — edge-case characterization before modular refactoring.
# Run from the repository root:
#   Rscript tests/scientific_regression/test_edge_cases.R
# This test script does not modify app/app.R or the captured regression RDS.
# It uses synthetic protein sequences; it is NOT experimental validation.

if (!file.exists('app/app.R')) {
  stop('Run from the SeekNSee repository root.', call. = FALSE)
}

required <- c('shiny', 'dplyr', 'stringr', 'writexl')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Missing R packages: ', paste(missing, collapse = ', '), call. = FALSE)

app <- new.env(parent = globalenv())
sys.source('app/app.R', envir = app, keep.source = FALSE)

n_pass <- 0L
check <- function(label, result) {
  value <- tryCatch(force(result), error = function(e) {
    stop('FAIL ', label, ': ', conditionMessage(e), call. = FALSE)
  })
  if (!identical(value, TRUE)) {
    stop('FAIL ', label, ': expected TRUE, received ', paste(value, collapse = ','), call. = FALSE)
  }
  n_pass <<- n_pass + 1L
  cat('PASS ', label, '\n', sep = '')
  invisible(TRUE)
}
throws <- function(expr, contains = NULL) {
  message <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) conditionMessage(e))
  !is.null(message) && (is.null(contains) || grepl(contains, message, fixed = TRUE))
}

# Fixtures live outside the repository; do not add test runs to Git.
fixture_dir <- tempfile('seeknsee_edge_cases_')
dir.create(fixture_dir)
fixture <- function(filename, contents) {
  p <- file.path(fixture_dir, filename)
  writeLines(contents, con = p, useBytes = TRUE)
  p
}

cat('\n== FASTA normalization and error handling ==\n')
duplicate_file <- fixture('duplicates.fasta', c('>DUP', 'AAAA', '>DUP', 'AAAA*'))
duplicates <- app$read_fasta(duplicate_file)
check('duplicate IDs are disambiguated', identical(duplicates$id, c('DUP', 'DUP__duplicate_1')))
check('original duplicate IDs are retained', identical(duplicates$original_id, c('DUP', 'DUP')))
check('terminal stop markers are removed', identical(duplicates$sequence, c('AAAA', 'AAAA')))
check('malformed FASTA without header rejected', throws(app$read_fasta(fixture('no_header.fasta', 'AAAA')), 'No FASTA headers'))
check('FASTA sequence preceding header rejected', throws(app$read_fasta(fixture('before_header.fasta', c('AAAA', '>A', 'AAAA'))), 'before its first header'))
check('internal stop codon rejected', throws(app$read_fasta(fixture('internal_stop.fasta', c('>A', 'AAA*AAA'))), 'Internal stop'))
check('empty FASTA protein rejected', throws(app$read_fasta(fixture('empty_entry.fasta', '>A')), 'Empty sequence'))
check('unsupported FASTA residue rejected', throws(app$read_fasta(fixture('bad_residue.fasta', c('>A', 'AAAA#'))), 'Unsupported character'))

cat('\n== Reference occurrence and specificity policies ==\n')
reference_file <- fixture('specificity_reference.fasta', c(
  '>sp|P00001|CAN GN=GENE1', 'MPEPTIDEKAAAAAPEPTIDEK',
  '>sp|P00001-2|ISO GN=GENE1', 'MPEPTIDEK',
  '>sp|P00002|OTHER GN=GENE2', 'QQQQQQQQ'
))
reference <- app$read_fasta(reference_file)
cache <- new.env(parent = emptyenv())
matches <- app$find_reference_matches('PEPTIDEK', reference, match_cache = cache)
check('two reference entries contain the peptide', as.integer(matches$entry_count) == 2L)
check('three total peptide occurrences are recorded', as.integer(matches$occurrence_count) == 3L)
check('two matching entries share one accession group', identical(matches$bases, 'P00001'))
check('reference hit IDs include canonical and isoform', identical(matches$ids, reference$id[1:2]))
check('peptide-match cache stores the result', exists('PEPTIDEK', envir = cache, inherits = FALSE))
check('cached result is identical to first search', identical(app$find_reference_matches('PEPTIDEK', reference, cache), matches))
check('absent reference substring yields zero occurrences', app$find_reference_matches('RRRRRRRR', reference)$occurrence_count == 0L)
check('full-length exact match includes the isoform entry', identical(app$exact_reference_ids('MPEPTIDEK', reference), reference$id[2]))

rules <- app$default_rules()
evaluate <- function(ref = reference, target_source = 'external_fasta',
                     target_id = 'EXTERNAL', target_sequence = 'MPEPTIDEK',
                     reference_policy = 'strict_absence') {
  app$evaluate_peptide(
    pep = 'PEPTIDEK', rules = rules, reference_proteome = ref,
    target_id = target_id, target_source = target_source,
    target_sequence = target_sequence, enzyme = 'Trypsin',
    start = 2L, end = 9L, reference_policy = reference_policy
  )
}
strict <- evaluate()
check('external strict absence blocks shared peptide', identical(strict$status, 'blocked') && 'unique_proteome' %in% strict$blocked_by)
check('accession-group rule passes for canonical plus isoform', identical(unname(strict$pass['unique_protein']), 'pass'))
allow_partial <- evaluate(reference_policy = 'allow_exact_target_match')
check('allowed exact full-sequence entry does not excuse other matches',
      identical(allow_partial$status, 'blocked') && identical(allow_partial$disallowed_reference_ids, reference$id[1]))
allow_only_exact <- evaluate(ref = reference[c(2L, 3L), , drop = FALSE], reference_policy = 'allow_exact_target_match')
check('external exact-target policy permits only the identical reference entry', identical(allow_only_exact$status, 'retained'))
reference_target <- evaluate(target_source = 'reference_proteome',
                             target_id = reference$id[1],
                             target_sequence = reference$sequence[1])
check('reference-proteome target is blocked by match in separate isoform entry',
      identical(reference_target$status, 'blocked') && identical(reference_target$disallowed_reference_ids, reference$id[2]))
absent <- evaluate(ref = reference[3, , drop = FALSE])
check('strict absence passes peptide absent from reference', identical(absent$status, 'retained'))

cat('\n== Conditional rules and scoring semantics ==\n')
check('Trypsin excludes K|P boundaries', !app$valid_tryptic_cut('AKPAAAA', 2L))
check('Trypsin accepts internal K|A boundary', app$valid_tryptic_cut('AKAAAA', 2L))
proline_risk <- app$assess_proline_p1prime('AKPAAAAA', 3L, 'Trypsin/P')
check('Trypsin/P P1-prime proline is applicable and fails',
      isTRUE(proline_risk$applicable) && identical(proline_risk$pass, FALSE))
check('P1-prime rule is inapplicable to standard Trypsin',
      !app$assess_proline_p1prime('AKPAAAAA', 3L, 'Trypsin')$applicable)
check('P1-prime rule is inapplicable to protein N-terminus',
      !app$assess_proline_p1prime('PAKAAAAA', 1L, 'Trypsin/P')$applicable)
acidic <- app$acidic_window_at_cut('AADEKAAAAA', 5L)
check('acidic window detects D/E at positions 3 and 4',
      isTRUE(acidic$risk) && identical(acidic$positions, c(3L, 4L)))
check('acidic rule is inapplicable outside Trypsin',
      !app$assess_acidic_near_tryptic_cleavage('AADEKAAAAA', 1L, 5L, 'AspN')$applicable)

all_ignore <- setNames(as.list(rep('ignore', length(app$RULE_NAMES))), app$RULE_NAMES)
ignored <- app$evaluate_peptide('AAAAAAAK', all_ignore, app$empty_proteome(),
  'SYN', 'external_fasta', 'AAAAAAAK', 'Trypsin')
check('ignored failing rules do not change score or eligibility',
      identical(ignored$status, 'retained') && identical(as.numeric(ignored$score), 12))
block_terminal <- rules
block_terminal$terminal <- 'block'
blocked_terminal <- app$evaluate_peptide('AAAAAAAK', block_terminal, app$empty_proteome(),
  'SYN', 'external_fasta', 'AAAAAAAK', 'Trypsin')
check('blocking terminal rule changes status without penalty',
      identical(blocked_terminal$status, 'blocked') &&
      identical(as.numeric(blocked_terminal$score), 12) &&
      unname(blocked_terminal$trace['terminal']) == 0)
applicable <- app$evaluate_peptide('AAAAAAAK', rules, app$empty_proteome(),
  'SYN', 'external_fasta', 'AAAAAAAK', 'AspN')
check('dibasic, acidic and P1-prime rules marked not applicable under AspN',
      identical(unname(applicable$pass[c('dibasic', 'acidic', 'proline_p1prime')]),
                rep('not_applicable', 3L)))
mc <- app$evaluate_peptide('AKRAK', rules, app$empty_proteome(),
  'SYN', 'external_fasta', 'AKRAK', 'Trypsin', missed_cleavages = 2L)
check('two missed cleavages contribute minus two', unname(mc$trace['miscleavage']) == -2)
check('candidate score equals 12 plus all trace contributions',
      identical(as.numeric(mc$score), as.numeric(12 + sum(mc$trace))))
check('glycosylation sequon is flagged',
      identical(unname(app$evaluate_peptide('NATAAAAK', rules, app$empty_proteome(),
        'SYN', 'external_fasta', 'NATAAAAK', 'Trypsin')$pass['glyco']), 'fail'))
check('dibasic motif is flagged for Trypsin',
      identical(unname(app$evaluate_peptide('AAAKKAAA', rules, app$empty_proteome(),
        'SYN', 'external_fasta', 'AAAKKAAA', 'Trypsin')$pass['dibasic']), 'fail'))
check('I/L residues are flagged by isobaric rule',
      identical(unname(app$evaluate_peptide('PEPTIDEK', rules, app$empty_proteome(),
        'SYN', 'external_fasta', 'PEPTIDEK', 'Trypsin')$pass['isobaric']), 'fail'))
charged <- app$evaluate_peptide('AKRKAAAA', rules, app$empty_proteome(),
  'SYN', 'external_fasta', 'AKRKAAAA', 'Trypsin')
check('heuristic charge rule flags z_hat > 2',
      identical(charged$charge_hat, 4L) && identical(unname(charged$pass['charge']), 'fail'))

cat('\n== Excel adaptation and preflight ==\n')
long_ids <- paste(sprintf('P%05d', 1:45), collapse = ';')
workbook <- list(
  Run_summary = data.frame(field = 'Run', value = 'synthetic', stringsAsFactors = FALSE),
  Ranked_peptides = data.frame(
    candidate_key = 'X1', target_id = 'SYN', peptide = 'AAAAAAAK',
    matching_reference_ids = long_ids, stringsAsFactors = FALSE
  )
)
adapted <- app$adapt_workbook_for_excel(workbook, safe_limit = 130L)
check('long identifier list triggers explicit Excel adaptation', isTRUE(adapted$truncation_used))
check('Excel adaptation reports omitted identifier occurrences',
      adapted$omitted_identifier_occurrences > 0L && nrow(adapted$audit) == 1L)
check('Excel identifier audit records total = exported + omitted',
      adapted$audit$total_ids[1] == adapted$audit$exported_ids[1] + adapted$audit$omitted_ids[1])
check('Excel adaptation preserves unrelated scientific fields',
      identical(adapted$workbook$Ranked_peptides$peptide, workbook$Ranked_peptides$peptide))
check('Excel adaptation documents abbreviation in Run_summary',
      'Identifier occurrences omitted from abbreviated cells' %in% adapted$workbook$Run_summary$field)
check('adapted identifier cell stays below requested safe limit',
      nchar(adapted$workbook$Ranked_peptides$matching_reference_ids[1]) <= 130L)

oversize <- list(Test = data.frame(candidate_key = 'X1',
                                    other_text = strrep('X', 32768L),
                                    stringsAsFactors = FALSE))
issues <- app$validate_excel_workbook(oversize)
check('Excel preflight locates cells over 32,767 characters',
      nrow(issues) == 1L && issues$characters[1] == 32768L)
check('non-listable oversized Excel cell raises an error',
      throws(app$adapt_workbook_for_excel(oversize), '32,767-character limit'))

output_path <- file.path(fixture_dir, 'export', 'small_report.xlsx')
exported_path <- app$export_seeknsee_workbook(workbook, output_path)
check('Excel export creates nonempty XLSX file',
      file.exists(exported_path) && file.info(exported_path)$size > 0)

cat('\n== Empty target edge case ==\n')
short_file <- fixture('short.fasta', c('>SHORT', 'AAA'))
short_run <- app$run_validation_case(short_file, target_source = 'reference_proteome',
  target_id = 'SHORT', enzyme = 'Trypsin')
check('short target generates zero 4–25-aa candidate instances',
      nrow(short_run$workbook$Ranked_peptides) == 0L &&
      nrow(short_run$workbook$Blocked_peptides) == 0L)
check('empty candidate run preserves accounting identities',
      all(app$validate_internal_identities(short_run)$passed))

unlink(fixture_dir, recursive = TRUE)
cat(sprintf('\nPASS EDGE CASES: %d checks passed.\n', n_pass))
cat('NOTE: These checks characterize v0.2.5 behavior, not peptide detectability.\n')
cat('NOT TESTED: Excel worksheet row limit; large-proteome memory/time; UI end-to-end; multi-OS packaging.\n')
