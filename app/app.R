library(shiny)
library(dplyr)
library(stringr)
library(writexl)

options(shiny.maxRequestSize = 500 * 1024^2)

SEEKNSEE_VERSION <- "0.2.5"
ALGORITHM_VERSION <- "2026-09-enzyme-aware-missed-cleavages-p1prime-1"

ENZYMES <- c(
  "Trypsin", "Trypsin/P", "LysC", "ArgC", "AspN",
  "Chymotrypsin", "Elastase", "GluC", "ProAlanase"
)

RULE_CHOICES <- c("ignore", "apply", "block")
RULE_NAMES <- c(
  "unique_proteome", "unique_protein", "length", "dibasic", "glyco",
  "chemical", "acidic", "terminal", "miscleavage", "proline_p1prime",
  "isobaric", "charge"
)

RULE_DEFINITIONS <- data.frame(
  rule = RULE_NAMES,
  label = c(
    "Specific against reference proteome", "Unique accession group", "Length",
    "Dibasic (trypsin)", "Glycosylation motif", "Chemical risk",
    "D/E near tryptic cleavage", "Terminal fragment", "Missed cleavages",
    "Proline at P1-prime (Trypsin/P)", "Isobaric I/L", "Heuristic charge > 2"
  ),
  definition = c(
    "For a reference target, fails when the peptide occurs in a non-target reference entry. For an external target under strict absence, fails on any reference match. Under allow exact target match, sequence-identical reference entries are permitted.",
    "Passes when all matching reference entries belong to no more than one normalized accession group.",
    "Passes for peptides from 8 through 25 amino acids.",
    "Fails when KK, KR, RK, or RR is present. This rule is specific to Trypsin and Trypsin/P and is not applicable to other enzymes.",
    "Fails when an N-X-S/T motif is present and X is not proline.",
    "Fails when the peptide contains M, W, C, NG, QG, DP, or starts with Q/E because of N-terminal pyroglutamate risk.",
    "Fails when D or E occurs within the P3-P3' window of a valid tryptic cleavage site delimiting the peptide. K/R-P bonds are not treated as tryptic cleavage sites. This rule is currently applicable only to Trypsin.",
    "Fails when the peptide corresponds to the N-terminal or C-terminal fragment of the target protein.",
    "Passes with zero missed cleavages and fails with one or two. When configured as apply, the penalty is -1 per missed cleavage, with no additional internal-site penalty.",
    "For Trypsin/P candidates created at an internal K/R cleavage boundary, fails when the first peptide residue is proline, corresponding to K|P or R|P. It is not applicable to protein N-terminal candidates or other enzymes.",
    "Fails when I or L is present.",
    "Uses z_hat = 1 + nK + nR + nH and fails when z_hat is greater than 2."
  ),
  stringsAsFactors = FALSE
)

empty_proteome <- function() {
  data.frame(
    id = character(), original_id = character(), accession = character(),
    base_id = character(), isoform = character(), gene = character(),
    name = character(), sequence = character(), sequence_length = integer(),
    stringsAsFactors = FALSE
  )
}

extract_accession <- function(protein_id) {
  if (is.null(protein_id) || length(protein_id) == 0 || !nzchar(protein_id[1])) {
    return("unknown_protein")
  }
  token <- as.character(protein_id[1])
  parts <- strsplit(token, "\\|")[[1]]
  accession <- if (length(parts) >= 3 && parts[1] %in% c("sp", "tr") && nzchar(parts[2])) {
    parts[2]
  } else {
    token
  }
  accession <- gsub("[^A-Za-z0-9_-]", "_", accession)
  ifelse(nzchar(accession), accession, "unknown_protein")
}

parse_identifier <- function(id) {
  parts <- strsplit(id, "\\|")[[1]]
  is_uniprot <- length(parts) >= 3 && parts[1] %in% c("sp", "tr") && nzchar(parts[2])
  accession <- if (is_uniprot) parts[2] else id
  isoform <- if (is_uniprot && grepl("-[0-9]+$", accession)) {
    sub("^.*-([0-9]+)$", "\\1", accession)
  } else {
    ""
  }
  base_id <- if (is_uniprot) sub("-[0-9]+$", "", accession) else accession
  list(accession = accession, base_id = base_id, isoform = isoform)
}

make_export_filename <- function(target_label) {
  paste0(
    extract_accession(target_label),
    "_peptide_selection_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".xlsx"
  )
}

read_fasta <- function(
    path,
    progress_callback = NULL
) {
  lines <- readLines(
    path,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  lines <- sub(
    "^\\ufeff",
    "",
    lines
  )
  
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  
  headers <- which(
    startsWith(
      lines,
      ">"
    )
  )
  
  if (length(headers) == 0) {
    stop(
      "No FASTA headers were found. Every entry must begin with '>'.",
      call. = FALSE
    )
  }
  
  if (headers[1] != 1) {
    stop(
      "The FASTA contains sequence text before its first header.",
      call. = FALSE
    )
  }
  
  if (!is.null(progress_callback)) {
    progress_callback(
      completed = 0L,
      total = length(headers)
    )
  }
  
  records <- vector(
    "list",
    length(headers)
  )
  
  for (i in seq_along(headers)) {
    start <- headers[i]
    
    end <- if (i == length(headers)) {
      length(lines)
    } else {
      headers[i + 1] - 1
    }
    
    header <- sub(
      "^>",
      "",
      lines[start]
    )
    
    seq_lines <- if (end > start) {
      lines[(start + 1):end]
    } else {
      character()
    }
    
    raw_sequence <- paste(
      seq_lines,
      collapse = ""
    )
    
    raw_sequence <- toupper(
      gsub(
        "[[:space:]]",
        "",
        raw_sequence
      )
    )
    
    if (!nzchar(raw_sequence)) {
      stop(
        paste0(
          "Empty sequence for FASTA entry: ",
          header
        ),
        call. = FALSE
      )
    }
    
    if (
      grepl(
        "[^ACDEFGHIKLMNPQRSTVWYXBZJUO*]",
        raw_sequence
      )
    ) {
      invalid_characters <- unique(
        strsplit(
          gsub(
            "[ACDEFGHIKLMNPQRSTVWYXBZJUO*]",
            "",
            raw_sequence
          ),
          ""
        )[[1]]
      )
      
      stop(
        paste0(
          "Unsupported character in sequence: ",
          header,
          ". Invalid character(s): ",
          paste(
            invalid_characters,
            collapse = ", "
          )
        ),
        call. = FALSE
      )
    }
    
    if (
      grepl(
        "\\*.",
        raw_sequence
      )
    ) {
      stop(
        paste0(
          "Internal stop character '*' in sequence: ",
          header
        ),
        call. = FALSE
      )
    }
    
    sequence <- sub(
      "\\*+$",
      "",
      raw_sequence
    )
    
    if (!nzchar(sequence)) {
      stop(
        paste0(
          "Empty sequence after removing terminal stop character: ",
          header
        ),
        call. = FALSE
      )
    }
    
    original_id <- strsplit(
      header,
      "[[:space:]]+"
    )[[1]][1]
    
    parsed <- parse_identifier(
      original_id
    )
    
    gene_match <- str_match(
      header,
      "(?:^|[[:space:]])GN=([A-Za-z0-9_.:-]+)"
    )[, 2]
    
    records[[i]] <- data.frame(
      original_id = original_id,
      accession = parsed$accession,
      base_id = parsed$base_id,
      isoform = parsed$isoform,
      gene = ifelse(
        is.na(gene_match),
        "",
        gene_match
      ),
      name = header,
      sequence = sequence,
      sequence_length = nchar(sequence),
      stringsAsFactors = FALSE
    )
    
    if (
      !is.null(progress_callback) &&
      (
        i == 1L ||
        i %% 100L == 0L ||
        i == length(headers)
      )
    ) {
      progress_callback(
        completed = i,
        total = length(headers)
      )
    }
  }
  
  out <- bind_rows(records)
  
  out$id <- make.unique(
    out$original_id,
    sep = "__duplicate_"
  )
  
  out %>%
    select(
      id,
      everything()
    )
}

find_regex_positions <- function(sequence, pattern) {
  hits <- gregexpr(pattern, sequence, perl = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) integer(0) else as.integer(hits)
}

cleavage_positions <- function(sequence, enzyme) {
  n <- nchar(sequence)

  cuts <- switch(
    enzyme,
    "Trypsin" = {
      locs <- find_regex_positions(sequence, "[KR]")
      locs[
        locs == n |
          str_sub(sequence, locs + 1L, locs + 1L) != "P"
      ]
    },
    "Trypsin/P" = find_regex_positions(sequence, "[KR]"),
    "LysC" = find_regex_positions(sequence, "K"),
    "ArgC" = find_regex_positions(sequence, "R"),
    "AspN" = find_regex_positions(sequence, "D") - 1L,
    "GluC" = find_regex_positions(sequence, "E"),
    "Chymotrypsin" = {
      locs <- find_regex_positions(sequence, "[FWY]")
      locs[
        locs == n |
          str_sub(sequence, locs + 1L, locs + 1L) != "P"
      ]
    },
    "Elastase" = find_regex_positions(sequence, "[AVLISG]"),
    "ProAlanase" = find_regex_positions(sequence, "[PA]"),
    integer(0)
  )

  sort(unique(cuts[cuts > 0L & cuts < n]))
}

digest_protein <- function(
    sequence,
    enzyme,
    missed_cleavages = 0L
) {
  n <- nchar(sequence)

  if (
    length(missed_cleavages) != 1L ||
    is.na(missed_cleavages) ||
    !is.numeric(missed_cleavages) ||
    missed_cleavages != as.integer(missed_cleavages) ||
    missed_cleavages < 0L ||
    missed_cleavages > 2L
  ) {
    stop(
      "missed_cleavages must be a single integer from 0 to 2.",
      call. = FALSE
    )
  }

  missed_cleavages <- as.integer(missed_cleavages)
  cuts <- cleavage_positions(sequence, enzyme)
  boundaries <- c(0L, cuts, n)
  fragment_starts <- boundaries[-length(boundaries)] + 1L
  fragment_ends <- boundaries[-1L]

  rows <- list()
  row_index <- 1L

  for (start_index in seq_along(fragment_starts)) {
    maximum_end_index <- min(
      length(fragment_starts),
      start_index + missed_cleavages
    )

    for (end_index in seq.int(start_index, maximum_end_index)) {
      peptide_start <- fragment_starts[start_index]
      peptide_end <- fragment_ends[end_index]

      rows[[row_index]] <- data.frame(
        peptide = substr(sequence, peptide_start, peptide_end),
        start = as.integer(peptide_start),
        end = as.integer(peptide_end),
        length = as.integer(peptide_end - peptide_start + 1L),
        missed_cleavages = as.integer(end_index - start_index),
        stringsAsFactors = FALSE
      )

      row_index <- row_index + 1L
    }
  }

  bind_rows(rows)
}

has_internal_cleavage_site <- function(pep, enzyme) {
  n <- nchar(pep)
  if (n <= 1L) return(FALSE)

  if (enzyme %in% c("Trypsin", "Trypsin/P")) {
    locs <- find_regex_positions(pep, "[KR]")

    if (identical(enzyme, "Trypsin")) {
      return(any(
        locs < n &
          str_sub(pep, locs + 1L, locs + 1L) != "P"
      ))
    }

    return(any(locs < n))
  }

  if (identical(enzyme, "AspN")) {
    return(any(find_regex_positions(pep, "D") > 1L))
  }

  if (identical(enzyme, "Chymotrypsin")) {
    locs <- find_regex_positions(pep, "[FWY]")
    return(any(
      locs < n &
        str_sub(pep, locs + 1L, locs + 1L) != "P"
    ))
  }

  pattern <- switch(
    enzyme,
    "LysC" = "K",
    "ArgC" = "R",
    "GluC" = "E",
    "Elastase" = "[AVLISG]",
    "ProAlanase" = "[PA]",
    NULL
  )

  if (is.null(pattern)) return(FALSE)
  any(find_regex_positions(pep, pattern) < n)
}

count_fixed_occurrences <- function(pattern, sequence) {
  hits <- gregexpr(pattern, sequence, fixed = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) 0L else length(hits)
}

find_reference_matches <- function(peptide, reference_proteome, match_cache = NULL) {
  if (!is.null(match_cache) && exists(peptide, envir = match_cache, inherits = FALSE)) {
    return(get(peptide, envir = match_cache, inherits = FALSE))
  }
  
  empty_result <- list(
    entry_count = 0L,
    occurrence_count = 0L,
    ids = character(),
    genes = character(),
    bases = character()
  )
  
  if (is.null(reference_proteome) || nrow(reference_proteome) == 0) {
    if (!is.null(match_cache)) assign(peptide, empty_result, envir = match_cache)
    return(empty_result)
  }
  
  hit <- grepl(peptide, reference_proteome$sequence, fixed = TRUE)
  
  if (!any(hit)) {
    if (!is.null(match_cache)) assign(peptide, empty_result, envir = match_cache)
    return(empty_result)
  }
  
  occurrence_counts <- vapply(
    reference_proteome$sequence[hit],
    function(sequence) count_fixed_occurrences(peptide, sequence),
    integer(1)
  )
  
  result <- list(
    entry_count = sum(hit),
    occurrence_count = sum(occurrence_counts),
    ids = reference_proteome$id[hit],
    genes = unique(reference_proteome$gene[hit & nzchar(reference_proteome$gene)]),
    bases = unique(reference_proteome$base_id[hit])
  )
  
  if (!is.null(match_cache)) assign(peptide, result, envir = match_cache)
  result
}

valid_tryptic_cut <- function(sequence, cut_after) {
  n <- nchar(sequence)
  if (length(cut_after) != 1L || is.na(cut_after) || cut_after < 1L || cut_after >= n) {
    return(FALSE)
  }
  p1 <- substr(sequence, cut_after, cut_after)
  p1_prime <- substr(sequence, cut_after + 1L, cut_after + 1L)
  p1 %in% c("K", "R") && p1_prime != "P"
}

acidic_window_at_cut <- function(sequence, cut_after) {
  if (!valid_tryptic_cut(sequence, cut_after)) {
    return(list(valid_site = FALSE, risk = FALSE, positions = integer(), residues = character()))
  }
  n <- nchar(sequence)
  positions <- c(cut_after - 2L, cut_after - 1L, cut_after + 1L, cut_after + 2L, cut_after + 3L)
  positions <- positions[positions >= 1L & positions <= n]
  residues <- substring(sequence, positions, positions)
  acidic <- residues %in% c("D", "E")
  list(
    valid_site = TRUE,
    risk = any(acidic),
    positions = as.integer(positions[acidic]),
    residues = residues[acidic]
  )
}

assess_acidic_near_tryptic_cleavage <- function(target_sequence, start, end, enzyme) {
  if (!identical(enzyme, "Trypsin")) {
    return(list(
      applicable = FALSE, pass = NA, evaluated_site_count = 0L,
      affected_boundaries = character(), affected_positions = integer(),
      affected_residues = character()
    ))
  }
  n <- nchar(target_sequence)
  boundaries <- character()
  cuts <- integer()
  if (start > 1L) {
    boundaries <- c(boundaries, "N")
    cuts <- c(cuts, start - 1L)
  }
  if (end < n) {
    boundaries <- c(boundaries, "C")
    cuts <- c(cuts, end)
  }
  details <- lapply(cuts, function(cut) acidic_window_at_cut(target_sequence, cut))
  valid <- if (length(details)) vapply(details, function(x) isTRUE(x$valid_site), logical(1)) else logical()
  details <- details[valid]
  boundaries <- boundaries[valid]
  risks <- if (length(details)) vapply(details, function(x) isTRUE(x$risk), logical(1)) else logical()
  list(
    applicable = TRUE,
    pass = !any(risks),
    evaluated_site_count = as.integer(length(details)),
    affected_boundaries = boundaries[risks],
    affected_positions = as.integer(unlist(lapply(details[risks], `[[`, "positions"), use.names = FALSE)),
    affected_residues = unlist(lapply(details[risks], `[[`, "residues"), use.names = FALSE)
  )
}

assess_proline_p1prime <- function(
    target_sequence,
    start,
    enzyme
) {
  if (!identical(enzyme, "Trypsin/P") || start <= 1L) {
    return(list(
      applicable = FALSE,
      pass = NA,
      preceding_residue = "",
      p1prime_residue = if (start >= 1L) substr(target_sequence, start, start) else ""
    ))
  }

  preceding_residue <- substr(
    target_sequence,
    start - 1L,
    start - 1L
  )

  p1prime_residue <- substr(
    target_sequence,
    start,
    start
  )

  applicable <- preceding_residue %in% c("K", "R")

  list(
    applicable = applicable,
    pass = if (applicable) p1prime_residue != "P" else NA,
    preceding_residue = preceding_residue,
    p1prime_residue = p1prime_residue
  )
}

exact_reference_ids <- function(target_sequence, reference_proteome) {
  if (is.null(reference_proteome) || nrow(reference_proteome) == 0) {
    return(character())
  }
  reference_proteome$id[reference_proteome$sequence == target_sequence]
}

evaluate_peptide <- function(
    pep,
    rules,
    reference_proteome,
    target_id,
    target_source,
    target_sequence,
    enzyme,
    start = 1L,
    end = nchar(target_sequence),
    target_sequence_length = nchar(target_sequence),
    missed_cleavages = 0L,
    reference_policy = "strict_absence",
    match_cache = NULL,
    precomputed_allowed_ids = NULL
) {
  missed_cleavages <- as.integer(missed_cleavages)

  if (
    length(missed_cleavages) != 1L ||
    is.na(missed_cleavages) ||
    missed_cleavages < 0L ||
    missed_cleavages > 2L
  ) {
    stop(
      "missed_cleavages must be a single integer from 0 to 2.",
      call. = FALSE
    )
  }

  match_info <- find_reference_matches(pep, reference_proteome, match_cache)

  if (!is.null(precomputed_allowed_ids)) {
    allowed_ids <- precomputed_allowed_ids
  } else if (identical(target_source, "reference_proteome")) {
    allowed_ids <- target_id
  } else if (identical(reference_policy, "allow_exact_target_match")) {
    allowed_ids <- exact_reference_ids(target_sequence, reference_proteome)
  } else {
    allowed_ids <- character()
  }

  disallowed_ids <- setdiff(match_info$ids, allowed_ids)
  reference_specific <- length(disallowed_ids) == 0L
  accession_group_unique <- length(unique(match_info$bases)) <= 1L
  charge_hat <- 1L + str_count(pep, "K") + str_count(pep, "R") + str_count(pep, "H")

  acidic_assessment <- assess_acidic_near_tryptic_cleavage(
    target_sequence,
    start,
    end,
    enzyme
  )

  proline_p1prime_assessment <- assess_proline_p1prime(
    target_sequence,
    start,
    enzyme
  )

  checks <- c(
    unique_proteome = reference_specific,
    unique_protein = accession_group_unique,
    length = nchar(pep) >= 8L && nchar(pep) <= 25L,
    dibasic = if (enzyme %in% c("Trypsin", "Trypsin/P")) {
      !grepl("KK|KR|RK|RR", pep)
    } else {
      NA
    },
    glyco = !grepl("N[^P][ST]", pep),
    chemical = !grepl("M|W|C|NG|QG|DP|^[QE]", pep),
    acidic = if (isTRUE(acidic_assessment$applicable)) acidic_assessment$pass else NA,
    terminal = !(start == 1L || end == target_sequence_length),
    miscleavage = missed_cleavages == 0L,
    proline_p1prime = if (isTRUE(proline_p1prime_assessment$applicable)) {
      proline_p1prime_assessment$pass
    } else {
      NA
    },
    isobaric = !grepl("[IL]", pep),
    charge = charge_hat <= 2L
  )

  applicable <- !is.na(checks)
  failed <- applicable & !checks
  states <- vapply(names(checks), function(rule) rules[[rule]], character(1))
  blocked_by <- names(checks)[failed & states == "block"]

  penalties <- setNames(numeric(length(checks)), names(checks))
  penalties[failed & states == "apply"] <- 1

  if (
    isTRUE(failed[["miscleavage"]]) &&
    identical(states[["miscleavage"]], "apply")
  ) {
    penalties[["miscleavage"]] <- missed_cleavages
  }

  trace <- -penalties
  score <- length(checks) + sum(trace)
  pass_fail <- ifelse(
    !applicable,
    "not_applicable",
    ifelse(checks, "pass", "fail")
  )
  names(pass_fail) <- names(checks)

  list(
    valid = length(blocked_by) == 0L,
    status = if (length(blocked_by) == 0L) "retained" else "blocked",
    score = as.numeric(score),
    trace = trace,
    pass = pass_fail,
    applicable = setNames(applicable, names(checks)),
    blocked_by = blocked_by,
    charge_hat = as.integer(charge_hat),
    missed_cleavages = missed_cleavages,
    acidic_assessment = acidic_assessment,
    proline_p1prime_assessment = proline_p1prime_assessment,
    match_info = match_info,
    disallowed_reference_ids = disallowed_ids,
    allowed_reference_ids = allowed_ids
  )
}

build_rules <- function(input) {
  setNames(
    lapply(RULE_NAMES, function(r) input[[paste0("rule_", r)]]),
    RULE_NAMES
  )
}

make_target_candidates <- function(
    targets,
    source_label,
    enzyme,
    rules,
    reference_proteome,
    reference_policy,
    missed_cleavages = 0L,
    match_cache = NULL
) {
  candidate_rows <- list()
  trace_rows <- list()
  raw_counts <- integer(nrow(targets))
  candidate_counts <- integer(nrow(targets))

  for (i in seq_len(nrow(targets))) {
    target <- targets[i, , drop = FALSE]

    digestion <- digest_protein(
      target$sequence,
      enzyme,
      missed_cleavages = missed_cleavages
    )

    raw_counts[i] <- nrow(digestion)
    candidates <- digestion %>% filter(length >= 4L, length <= 25L)
    candidate_counts[i] <- nrow(candidates)

    if (nrow(candidates) == 0L) next

    allowed_ids_for_target <- if (identical(source_label, "reference_proteome")) {
      target$id
    } else if (identical(reference_policy, "allow_exact_target_match")) {
      exact_reference_ids(target$sequence, reference_proteome)
    } else {
      character()
    }

    for (j in seq_len(nrow(candidates))) {
      cand <- candidates[j, , drop = FALSE]

      ev <- evaluate_peptide(
        pep = cand$peptide,
        rules = rules,
        reference_proteome = reference_proteome,
        target_id = target$id,
        target_source = source_label,
        target_sequence = target$sequence,
        enzyme = enzyme,
        start = cand$start,
        end = cand$end,
        target_sequence_length = target$sequence_length,
        missed_cleavages = cand$missed_cleavages,
        reference_policy = reference_policy,
        match_cache = match_cache,
        precomputed_allowed_ids = allowed_ids_for_target
      )

      candidate_key <- paste(
        target$id,
        cand$start,
        cand$end,
        cand$missed_cleavages,
        cand$peptide,
        sep = "|"
      )

      candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
        candidate_key = candidate_key,
        target_id = target$id,
        target_original_id = target$original_id,
        target_accession = target$accession,
        target_gene = target$gene,
        target_source = source_label,
        peptide = cand$peptide,
        start = cand$start,
        end = cand$end,
        length = cand$length,
        enzyme = enzyme,
        missed_cleavages = cand$missed_cleavages,
        candidate_status = ev$status,
        peptide_score = ev$score,
        blocked_by = paste(ev$blocked_by, collapse = ";"),
        charge_hat = ev$charge_hat,
        proline_p1prime_applicable = ev$proline_p1prime_assessment$applicable,
        proline_p1prime_preceding_residue = ev$proline_p1prime_assessment$preceding_residue,
        proline_p1prime_residue = ev$proline_p1prime_assessment$p1prime_residue,
        acidic_evaluated_site_count = ev$acidic_assessment$evaluated_site_count,
        acidic_affected_boundaries = paste(ev$acidic_assessment$affected_boundaries, collapse = ";"),
        acidic_affected_positions = paste(ev$acidic_assessment$affected_positions, collapse = ";"),
        acidic_affected_residues = paste(ev$acidic_assessment$affected_residues, collapse = ";"),
        reference_entry_count = ev$match_info$entry_count,
        reference_occurrence_count = ev$match_info$occurrence_count,
        matching_reference_ids = paste(ev$match_info$ids, collapse = ";"),
        disallowed_reference_ids = paste(ev$disallowed_reference_ids, collapse = ";"),
        stringsAsFactors = FALSE
      )

      for (rule in RULE_NAMES) {
        trace_rows[[length(trace_rows) + 1L]] <- data.frame(
          candidate_key = candidate_key,
          target_id = target$id,
          peptide = cand$peptide,
          start = cand$start,
          end = cand$end,
          missed_cleavages = cand$missed_cleavages,
          rule = rule,
          state = rules[[rule]],
          applicable = as.logical(ev$applicable[rule]),
          pass_fail = unname(ev$pass[rule]),
          rule_score = as.numeric(ev$trace[rule]),
          candidate_status = ev$status,
          peptide_score = ev$score,
          blocks_candidate = as.logical(rule %in% ev$blocked_by),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  candidates_df <- if (length(candidate_rows)) bind_rows(candidate_rows) else data.frame()
  trace_df <- if (length(trace_rows)) bind_rows(trace_rows) else data.frame()

  target_summary <- targets %>%
    transmute(
      target_id = id,
      target_original_id = original_id,
      target_accession = accession,
      target_gene = gene,
      target_source = source_label,
      sequence_length = sequence_length,
      maximum_missed_cleavages = as.integer(missed_cleavages),
      exact_reference_match_ids = vapply(
        sequence,
        function(x) paste(exact_reference_ids(x, reference_proteome), collapse = ";"),
        character(1)
      ),
      raw_fragments = raw_counts,
      candidates_4_25 = candidate_counts
    )

  if (nrow(candidates_df)) {
    retained_counts <- candidates_df %>%
      filter(candidate_status == "retained") %>%
      count(target_id, name = "retained")

    blocked_counts <- candidates_df %>%
      filter(candidate_status == "blocked") %>%
      count(target_id, name = "blocked")

    target_summary <- target_summary %>%
      left_join(retained_counts, by = "target_id") %>%
      left_join(blocked_counts, by = "target_id")
  }

  if (!"retained" %in% names(target_summary)) target_summary$retained <- 0L
  if (!"blocked" %in% names(target_summary)) target_summary$blocked <- 0L
  target_summary$retained[is.na(target_summary$retained)] <- 0L
  target_summary$blocked[is.na(target_summary$blocked)] <- 0L

  list(
    candidates = candidates_df,
    rule_trace = trace_df,
    target_summary = target_summary
  )
}


##############################
# VALIDATION AND SHARED EXPORT HELPERS
##############################

default_rules <- function() {
  setNames(
    as.list(c(
      "block", # unique_proteome
      "apply", # unique_protein
      "block", # length
      "apply", # dibasic
      "apply", # glyco
      "apply", # chemical
      "apply", # acidic
      "apply", # terminal
      "apply", # miscleavage
      "apply", # proline_p1prime
      "apply", # isobaric
      "apply"  # charge
    )),
    RULE_NAMES
  )
}

validate_rule_configuration <- function(rules) {
  if (!is.list(rules)) {
    stop("rules must be a named list.", call. = FALSE)
  }
  
  missing_rules <- setdiff(RULE_NAMES, names(rules))
  extra_rules <- setdiff(names(rules), RULE_NAMES)
  
  if (length(missing_rules) > 0) {
    stop(
      "Missing rules: ",
      paste(missing_rules, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (length(extra_rules) > 0) {
    stop(
      "Unknown rules: ",
      paste(extra_rules, collapse = ", "),
      call. = FALSE
    )
  }
  
  states <- unlist(rules[RULE_NAMES], use.names = FALSE)
  invalid_states <- setdiff(states, RULE_CHOICES)
  
  if (length(invalid_states) > 0) {
    stop(
      "Invalid rule states: ",
      paste(unique(invalid_states), collapse = ", "),
      ". Allowed values: ",
      paste(RULE_CHOICES, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

empty_candidate_table <- function() {
  data.frame(
    candidate_key = character(),
    target_id = character(),
    target_original_id = character(),
    target_accession = character(),
    target_gene = character(),
    target_source = character(),
    peptide = character(),
    start = integer(),
    end = integer(),
    length = integer(),
    enzyme = character(),
    missed_cleavages = integer(),
    candidate_status = character(),
    peptide_score = numeric(),
    blocked_by = character(),
    charge_hat = integer(),
    proline_p1prime_applicable = logical(),
    proline_p1prime_preceding_residue = character(),
    proline_p1prime_residue = character(),
    acidic_evaluated_site_count = integer(),
    acidic_affected_boundaries = character(),
    acidic_affected_positions = character(),
    acidic_affected_residues = character(),
    reference_entry_count = integer(),
    reference_occurrence_count = integer(),
    matching_reference_ids = character(),
    disallowed_reference_ids = character(),
    stringsAsFactors = FALSE
  )
}

empty_rule_trace_table <- function() {
  data.frame(
    candidate_key = character(),
    target_id = character(),
    peptide = character(),
    start = integer(),
    end = integer(),
    missed_cleavages = integer(),
    rule = character(),
    state = character(),
    applicable = logical(),
    pass_fail = character(),
    rule_score = numeric(),
    candidate_status = character(),
    peptide_score = numeric(),
    blocks_candidate = logical(),
    stringsAsFactors = FALSE
  )
}

build_seeknsee_workbook <- function(result, context, reference_proteome) {
  candidates <- result$candidates
  if (is.null(candidates) || ncol(candidates) == 0) {
    candidates <- empty_candidate_table()
  }
  
  rule_trace <- result$rule_trace
  if (is.null(rule_trace) || ncol(rule_trace) == 0) {
    rule_trace <- empty_rule_trace_table()
  }
  
  ranked <- candidates %>%
    filter(candidate_status == "retained") %>%
    arrange(desc(peptide_score), target_id, peptide, start)
  
  blocked <- candidates %>%
    filter(candidate_status == "blocked") %>%
    arrange(target_id, peptide, start)
  
  rules_df <- data.frame(
    rule = RULE_NAMES,
    state = unlist(context$rules[RULE_NAMES], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  retained_count <- sum(candidates$candidate_status == "retained")
  blocked_count <- sum(candidates$candidate_status == "blocked")
  
  summary_df <- data.frame(
    field = c(
      "SeekNSee version",
      "Algorithm version",
      "R version",
      "Run date/time",
      "Reference FASTA",
      "Reference proteins",
      "Target FASTA",
      "Target source",
      "Reference specificity policy",
      "Enzyme",
      "Maximum missed cleavages",
      "Targets processed",
      "Raw fragments",
      "Candidates 4-25 aa",
      "Retained candidates",
      "Blocked candidates"
    ),
    value = c(
      SEEKNSEE_VERSION,
      ALGORITHM_VERSION,
      R.version.string,
      context$run_datetime,
      context$reference_filename,
      nrow(reference_proteome),
      context$target_filename,
      context$target_source,
      context$reference_policy,
      context$enzyme,
      context$missed_cleavages,
      nrow(result$target_summary),
      sum(result$target_summary$raw_fragments),
      sum(result$target_summary$candidates_4_25),
      retained_count,
      blocked_count
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    Run_summary = summary_df,
    Rule_configuration = rules_df,
    Target_proteins = result$target_summary,
    Ranked_peptides = ranked,
    Rule_trace = rule_trace,
    Blocked_peptides = blocked,
    Rule_definitions = RULE_DEFINITIONS
  )
}

EXCEL_MAX_CELL_CHARACTERS <- 32767L
EXCEL_SAFE_CELL_CHARACTERS <- 30000L
EXCEL_LIST_COLUMNS <- c(
  "matching_reference_ids",
  "disallowed_reference_ids",
  "exact_reference_match_ids"
)

truncate_delimited_identifier_list <- function(
    value,
    safe_limit = EXCEL_SAFE_CELL_CHARACTERS,
    separator = ";"
) {
  value <- if (is.na(value)) "" else as.character(value)
  
  if (!nzchar(value) || nchar(value, type = "chars") <= safe_limit) {
    total_ids <- if (nzchar(value)) {
      length(strsplit(value, separator, fixed = TRUE)[[1]])
    } else {
      0L
    }
    
    return(list(
      value = value,
      total_ids = as.integer(total_ids),
      exported_ids = as.integer(total_ids),
      omitted_ids = 0L,
      truncated = FALSE
    ))
  }
  
  identifiers <- strsplit(value, separator, fixed = TRUE)[[1]]
  total_ids <- length(identifiers)
  kept <- character()
  
  for (identifier in identifiers) {
    proposed_kept <- c(kept, identifier)
    omitted_if_kept <- total_ids - length(proposed_kept)
    suffix <- if (omitted_if_kept > 0L) {
      paste0(
        separator,
        "[",
        omitted_if_kept,
        " additional reference IDs omitted from this Excel cell]"
      )
    } else {
      ""
    }
    proposed_value <- paste0(
      paste(proposed_kept, collapse = separator),
      suffix
    )
    
    if (nchar(proposed_value, type = "chars") > safe_limit) {
      break
    }
    
    kept <- proposed_kept
  }
  
  if (length(kept) == 0L) {
    stop(
      "An individual reference identifier exceeds the safe Excel cell limit.",
      call. = FALSE
    )
  }
  
  omitted_ids <- total_ids - length(kept)
  suffix <- paste0(
    separator,
    "[",
    omitted_ids,
    " additional reference IDs omitted from this Excel cell]"
  )
  output_value <- paste0(
    paste(kept, collapse = separator),
    suffix
  )
  
  list(
    value = output_value,
    total_ids = as.integer(total_ids),
    exported_ids = as.integer(length(kept)),
    omitted_ids = as.integer(omitted_ids),
    truncated = TRUE
  )
}

adapt_sheet_identifier_lists_for_excel <- function(
    sheet,
    sheet_name,
    safe_limit = EXCEL_SAFE_CELL_CHARACTERS
) {
  if (!is.data.frame(sheet) || nrow(sheet) == 0L) {
    return(list(sheet = sheet, audit = data.frame()))
  }
  
  columns_to_adapt <- intersect(EXCEL_LIST_COLUMNS, names(sheet))
  if (length(columns_to_adapt) == 0L) {
    return(list(sheet = sheet, audit = data.frame()))
  }
  
  audit_parts <- list()
  
  for (column_name in columns_to_adapt) {
    adapted <- lapply(
      sheet[[column_name]],
      truncate_delimited_identifier_list,
      safe_limit = safe_limit
    )
    
    truncated_flags <- vapply(
      adapted,
      function(x) x$truncated,
      logical(1)
    )
    
    affected_rows <- which(truncated_flags)
    
    if (length(affected_rows) > 0L) {
      sheet[[column_name]] <- vapply(
        adapted,
        function(x) x$value,
        character(1)
      )
      
      sheet[[paste0(column_name, "_total")]] <- vapply(
        adapted,
        function(x) x$total_ids,
        integer(1)
      )
      
      sheet[[paste0(column_name, "_exported")]] <- vapply(
        adapted,
        function(x) x$exported_ids,
        integer(1)
      )
      
      sheet[[paste0(column_name, "_omitted")]] <- vapply(
        adapted,
        function(x) x$omitted_ids,
        integer(1)
      )
      
      sheet[[paste0(column_name, "_truncated")]] <- truncated_flags
      contextual_value <- function(name) {
        if (name %in% names(sheet)) {
          as.character(sheet[[name]][affected_rows])
        } else {
          rep("", length(affected_rows))
        }
      }
      
      audit_parts[[length(audit_parts) + 1L]] <- data.frame(
        sheet = rep(sheet_name, length(affected_rows)),
        column = rep(column_name, length(affected_rows)),
        row = as.integer(affected_rows),
        candidate_key = contextual_value("candidate_key"),
        target_id = contextual_value("target_id"),
        peptide = contextual_value("peptide"),
        total_ids = sheet[[paste0(column_name, "_total")]][affected_rows],
        exported_ids = sheet[[paste0(column_name, "_exported")]][affected_rows],
        omitted_ids = sheet[[paste0(column_name, "_omitted")]][affected_rows],
        stringsAsFactors = FALSE
      )
    }
  }
  
  audit <- if (length(audit_parts) > 0L) {
    bind_rows(audit_parts)
  } else {
    data.frame()
  }
  
  list(sheet = sheet, audit = audit)
}

adapt_workbook_for_excel <- function(
    workbook,
    safe_limit = EXCEL_SAFE_CELL_CHARACTERS
) {
  if (!is.list(workbook) || is.null(names(workbook))) {
    stop(
      "workbook must be a named list of data frames.",
      call. = FALSE
    )
  }
  
  adapted_workbook <- workbook
  audit_parts <- list()
  
  for (sheet_name in names(adapted_workbook)) {
    adapted <- adapt_sheet_identifier_lists_for_excel(
      sheet = adapted_workbook[[sheet_name]],
      sheet_name = sheet_name,
      safe_limit = safe_limit
    )
    
    adapted_workbook[[sheet_name]] <- adapted$sheet
    
    if (nrow(adapted$audit) > 0L) {
      audit_parts[[length(audit_parts) + 1L]] <- adapted$audit
    }
  }
  
  audit <- if (length(audit_parts) > 0L) {
    bind_rows(audit_parts) %>%
      arrange(sheet, column, row)
  } else {
    data.frame(
      sheet = character(),
      column = character(),
      row = integer(),
      candidate_key = character(),
      target_id = character(),
      peptide = character(),
      total_ids = integer(),
      exported_ids = integer(),
      omitted_ids = integer(),
      stringsAsFactors = FALSE
    )
  }
  
  truncation_used <- nrow(audit) > 0L
  total_omitted <- if (truncation_used) sum(audit$omitted_ids) else 0L
  
  if (truncation_used && "Run_summary" %in% names(adapted_workbook)) {
    export_summary <- data.frame(
      field = c(
        "Excel identifier-list adaptation used",
        "Excel maximum characters per cell",
        "SeekNSee safe characters per cell",
        "Identifier-list cells abbreviated",
        "Identifier occurrences omitted from abbreviated cells"
      ),
      value = c(
        if (truncation_used) "TRUE" else "FALSE",
        EXCEL_MAX_CELL_CHARACTERS,
        safe_limit,
        nrow(audit),
        total_omitted
      ),
      stringsAsFactors = FALSE
    )
    
    adapted_workbook$Run_summary <- bind_rows(
      adapted_workbook$Run_summary,
      export_summary
    )
  }
  
  remaining_issues <- validate_excel_workbook(
    adapted_workbook,
    limit = EXCEL_MAX_CELL_CHARACTERS
  )
  
  if (nrow(remaining_issues) > 0L) {
    stop(
      "Excel export adaptation could not reduce every cell below Excel's 32,767-character limit.",
      call. = FALSE
    )
  }
  
  list(
    workbook = adapted_workbook,
    audit = audit,
    truncation_used = truncation_used,
    abbreviated_cells = as.integer(nrow(audit)),
    omitted_identifier_occurrences = as.integer(total_omitted)
  )
}

empty_excel_preflight_issues <- function() {
  data.frame(
    sheet = character(),
    column = character(),
    row = integer(),
    characters = integer(),
    excel_limit = integer(),
    candidate_key = character(),
    target_id = character(),
    peptide = character(),
    stringsAsFactors = FALSE
  )
}

validate_excel_workbook <- function(
    workbook,
    limit = EXCEL_MAX_CELL_CHARACTERS
) {
  issue_parts <- list()
  
  for (sheet_name in names(workbook)) {
    sheet <- workbook[[sheet_name]]
    
    if (!is.data.frame(sheet) || nrow(sheet) == 0L || ncol(sheet) == 0L) {
      next
    }
    
    character_columns <- names(sheet)[
      vapply(sheet, is.character, logical(1))
    ]
    
    for (column_name in character_columns) {
      lengths <- nchar(
        sheet[[column_name]],
        type = "chars",
        allowNA = TRUE,
        keepNA = TRUE
      )
      affected_rows <- which(!is.na(lengths) & lengths > limit)
      
      if (length(affected_rows) == 0L) {
        next
      }
      
      contextual_value <- function(name) {
        if (name %in% names(sheet)) {
          as.character(sheet[[name]][affected_rows])
        } else {
          rep("", length(affected_rows))
        }
      }
      
      issue_parts[[length(issue_parts) + 1L]] <- data.frame(
        sheet = rep(sheet_name, length(affected_rows)),
        column = rep(column_name, length(affected_rows)),
        row = as.integer(affected_rows),
        characters = as.integer(lengths[affected_rows]),
        excel_limit = rep(as.integer(limit), length(affected_rows)),
        candidate_key = contextual_value("candidate_key"),
        target_id = contextual_value("target_id"),
        peptide = contextual_value("peptide"),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(issue_parts) == 0L) {
    return(empty_excel_preflight_issues())
  }
  
  bind_rows(issue_parts) %>%
    arrange(desc(characters), sheet, column, row)
}

export_seeknsee_workbook <- function(workbook, output_file) {
  output_directory <- dirname(output_file)
  if (!identical(output_directory, ".") && !dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  
  adapted_export <- adapt_workbook_for_excel(workbook)
  write_xlsx(adapted_export$workbook, path = output_file)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

run_validation_case <- function(
    reference_fasta,
    target_source = c("reference_proteome", "external_fasta"),
    target_id = NULL,
    target_fasta = NULL,
    target_scope = c("selected", "all"),
    selected_external_id = NULL,
    reference_policy = c("strict_absence", "allow_exact_target_match"),
    enzyme = "Trypsin",
    missed_cleavages = 0L,
    rules = default_rules(),
    output_file = NULL,
    reference_data = NULL,
    match_cache = NULL
) {
  target_source <- match.arg(target_source)
  target_scope <- match.arg(target_scope)
  reference_policy <- match.arg(reference_policy)
  validate_rule_configuration(rules)

  if (!enzyme %in% ENZYMES) {
    stop(
      "Unknown enzyme: ", enzyme,
      ". Allowed values: ", paste(ENZYMES, collapse = ", "),
      call. = FALSE
    )
  }

  missed_cleavages <- as.integer(missed_cleavages)
  if (
    length(missed_cleavages) != 1L ||
    is.na(missed_cleavages) ||
    missed_cleavages < 0L ||
    missed_cleavages > 2L
  ) {
    stop("missed_cleavages must be a single integer from 0 to 2.", call. = FALSE)
  }

  if (is.null(reference_data)) reference_data <- read_fasta(reference_fasta)
  if (is.null(match_cache)) match_cache <- new.env(parent = emptyenv())

  if (identical(target_source, "reference_proteome")) {
    if (is.null(target_id) || !nzchar(target_id)) {
      stop("target_id is required for a reference-proteome target.", call. = FALSE)
    }

    targets <- reference_data %>% filter(id == target_id)
    if (nrow(targets) == 0L) stop("Target not found in reference proteome: ", target_id, call. = FALSE)
    if (nrow(targets) > 1L) stop("Target ID is not unique: ", target_id, call. = FALSE)

    policy <- "allow_selected_reference_target"
    target_filename <- basename(reference_fasta)
  } else {
    if (is.null(target_fasta) || !nzchar(target_fasta)) {
      stop("target_fasta is required for external targets.", call. = FALSE)
    }

    external_data <- read_fasta(target_fasta)

    if (identical(target_scope, "all")) {
      targets <- external_data
    } else {
      if (is.null(selected_external_id) || !nzchar(selected_external_id)) {
        stop("selected_external_id is required when target_scope = 'selected'.", call. = FALSE)
      }

      targets <- external_data %>% filter(id == selected_external_id)
      if (nrow(targets) == 0L) stop("External target not found: ", selected_external_id, call. = FALSE)
      if (nrow(targets) > 1L) stop("External target ID is not unique: ", selected_external_id, call. = FALSE)
    }

    policy <- reference_policy
    target_filename <- basename(target_fasta)
  }

  start_time <- Sys.time()

  result <- make_target_candidates(
    targets = targets,
    source_label = target_source,
    enzyme = enzyme,
    rules = rules,
    reference_proteome = reference_data,
    reference_policy = policy,
    missed_cleavages = missed_cleavages,
    match_cache = match_cache
  )

  end_time <- Sys.time()

  context <- list(
    run_datetime = format(start_time, "%Y-%m-%d %H:%M:%S %z"),
    reference_filename = basename(reference_fasta),
    target_filename = target_filename,
    target_source = target_source,
    reference_policy = policy,
    enzyme = enzyme,
    missed_cleavages = missed_cleavages,
    rules = rules,
    elapsed_seconds = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  workbook <- build_seeknsee_workbook(
    result = result,
    context = context,
    reference_proteome = reference_data
  )

  exported_file <- NULL
  if (!is.null(output_file)) exported_file <- export_seeknsee_workbook(workbook, output_file)

  structure(
    list(
      result = result,
      context = context,
      workbook = workbook,
      reference_proteome = reference_data,
      match_cache = match_cache,
      exported_file = exported_file
    ),
    class = "seeknsee_validation"
  )
}

summarize_validation_result <- function(validation) {
  workbook <- validation$workbook
  ranked <- workbook$Ranked_peptides
  trace <- workbook$Rule_trace
  blocked <- workbook$Blocked_peptides
  targets <- workbook$Target_proteins
  
  score_summary <- if (nrow(ranked) > 0) {
    data.frame(
      minimum_score = min(ranked$peptide_score),
      maximum_score = max(ranked$peptide_score),
      median_score = median(ranked$peptide_score),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      minimum_score = NA_real_,
      maximum_score = NA_real_,
      median_score = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  score_distribution <- if (nrow(ranked) > 0) {
    ranked %>%
      count(peptide_score, name = "candidate_count") %>%
      arrange(desc(peptide_score))
  } else {
    data.frame(peptide_score = numeric(), candidate_count = integer())
  }
  
  rule_failures <- if (nrow(trace) > 0) {
    trace %>%
      group_by(rule) %>%
      summarise(
        tested = sum(applicable),
        not_applicable = sum(!applicable),
        failed = sum(applicable & pass_fail == "fail"),
        failure_percentage = ifelse(tested > 0, 100 * failed / tested, NA_real_),
        .groups = "drop"
      ) %>%
      arrange(match(rule, RULE_NAMES))
  } else {
    data.frame(
      rule = character(),
      tested = integer(),
      not_applicable = integer(),
      failed = integer(),
      failure_percentage = numeric()
    )
  }
  
  top_five <- if (nrow(ranked) > 0) {
    ranked %>%
      arrange(desc(peptide_score), target_id, peptide, start) %>%
      slice_head(n = 5) %>%
      select(target_id, peptide, start, end, peptide_score)
  } else {
    data.frame()
  }
  
  block_counts <- setNames(integer(length(RULE_NAMES)), RULE_NAMES)
  if (nrow(trace) > 0) {
    observed <- trace %>%
      filter(blocks_candidate) %>%
      count(rule, name = "blocked_candidates")
    block_counts[observed$rule] <- observed$blocked_candidates
  }
  
  multiple_blocks <- if (nrow(blocked) > 0) {
    sum(str_count(blocked$blocked_by, fixed(";")) >= 1)
  } else {
    0L
  }
  
  list(
    targets = targets,
    score_summary = score_summary,
    score_distribution = score_distribution,
    rule_failures = rule_failures,
    top_five = top_five,
    blocked_by_rule = data.frame(
      rule = names(block_counts),
      blocked_candidates = as.integer(block_counts),
      stringsAsFactors = FALSE
    ),
    blocked_total = nrow(blocked),
    multiple_blocks = multiple_blocks,
    rule_trace_rows = nrow(trace),
    elapsed_seconds = validation$context$elapsed_seconds
  )
}

validate_internal_identities <- function(validation) {
  workbook <- validation$workbook
  ranked <- workbook$Ranked_peptides
  blocked <- workbook$Blocked_peptides
  trace <- workbook$Rule_trace
  targets <- workbook$Target_proteins

  total_candidates <- sum(targets$candidates_4_25)
  retained_count <- nrow(ranked)
  blocked_count <- nrow(blocked)

  identity_candidates <- total_candidates == retained_count + blocked_count
  identity_trace <- nrow(trace) == total_candidates * length(RULE_NAMES)

  if (nrow(trace) > 0L) {
    score_check <- trace %>%
      group_by(candidate_key) %>%
      summarise(
        expected_score = length(RULE_NAMES) + sum(rule_score),
        observed_score = first(peptide_score),
        score_correct = expected_score == observed_score,
        rule_rows = n(),
        .groups = "drop"
      )

    identity_scores <- all(score_check$score_correct)
    identity_candidate_trace <- all(score_check$rule_rows == length(RULE_NAMES))
  } else {
    identity_scores <- TRUE
    identity_candidate_trace <- total_candidates == 0L
  }

  data.frame(
    identity = c(
      "candidates = retained + blocked",
      paste0("rule trace = candidates x ", length(RULE_NAMES)),
      paste0("each candidate has ", length(RULE_NAMES), " rule rows"),
      "score = rule count + sum of rule_score"
    ),
    passed = c(
      identity_candidates,
      identity_trace,
      identity_candidate_trace,
      identity_scores
    ),
    stringsAsFactors = FALSE
  )
}

evaluate_test_peptide <- function(
    peptide,
    rules = default_rules(),
    enzyme = "Trypsin",
    reference_proteome = empty_proteome(),
    target_id = "synthetic_target",
    target_source = "external_fasta",
    target_sequence = peptide,
    start = 1L,
    end = nchar(target_sequence),
    target_sequence_length = nchar(target_sequence),
    missed_cleavages = 0L,
    reference_policy = "strict_absence"
) {
  validate_rule_configuration(rules)

  evaluation <- evaluate_peptide(
    pep = peptide,
    rules = rules,
    reference_proteome = reference_proteome,
    target_id = target_id,
    target_source = target_source,
    target_sequence = target_sequence,
    enzyme = enzyme,
    start = start,
    end = end,
    target_sequence_length = target_sequence_length,
    missed_cleavages = missed_cleavages,
    reference_policy = reference_policy
  )

  data.frame(
    peptide = peptide,
    missed_cleavages = as.integer(missed_cleavages),
    rule = RULE_NAMES,
    state = unlist(rules[RULE_NAMES], use.names = FALSE),
    applicable = as.logical(evaluation$applicable[RULE_NAMES]),
    pass_fail = unname(evaluation$pass[RULE_NAMES]),
    rule_score = as.numeric(evaluation$trace[RULE_NAMES]),
    candidate_status = evaluation$status,
    peptide_score = evaluation$score,
    blocks_candidate = RULE_NAMES %in% evaluation$blocked_by,
    charge_hat = evaluation$charge_hat,
    stringsAsFactors = FALSE
  )
}

compare_digestion <- function(sequence, enzyme, expected_peptides) {
  if (!enzyme %in% ENZYMES) {
    stop("Unknown enzyme: ", enzyme, call. = FALSE)
  }
  
  observed <- digest_protein(sequence, enzyme)
  observed_peptides <- observed$peptide
  
  data.frame(
    enzyme = enzyme,
    sequence = sequence,
    expected = paste(expected_peptides, collapse = ";"),
    observed = paste(observed_peptides, collapse = ";"),
    expected_count = length(expected_peptides),
    observed_count = length(observed_peptides),
    concordance = identical(observed_peptides, expected_peptides),
    stringsAsFactors = FALSE
  )
}

ui <- fluidPage(
  tags$head(
    tags$title("SeekNSee"),
    tags$style(HTML("
      .sidebar-flex {
        display: flex;
        flex-direction: column;
        height: calc(100vh - 140px);
      }
      .sidebar-body {
        flex: 1 1 auto;
        overflow-y: auto;
        padding-bottom: 10px;
      }
      .sidebar-footer {
        flex: 0 0 auto;
        padding-top: 10px;
        border-top: 1px solid #e5e5e5;
      }
      .fullwidth-btn > .btn,
      .fullwidth-btn > a.btn {
        width: 100%;
      }
      .info-page {
        max-width: 950px;
        line-height: 1.6;
        font-size: 15px;
      }
      .info-page h3 {
        margin-top: 0;
        font-weight: bold;
      }
      .info-page h4 {
        margin-top: 22px;
        font-weight: bold;
      }
      .info-page ul,
      .info-page ol {
        margin-bottom: 14px;
      }
      .info-page code {
        background-color: #f5f5f5;
        padding: 2px 4px;
        border-radius: 3px;
        color: #c7254e;
      }
      .status-box {
        padding: 8px;
        margin: 8px 0 14px 0;
        background: #f5f5f5;
        border-left: 4px solid #337ab7;
      }
      .peptide-column {
        min-width: 135px;
      }
    "))
  ),
  
  fluidRow(
    column(
      12,
      div(
        style = "position:relative;height:80px;margin-bottom:10px;",
        img(
          src = "logo.png",
          height = "72px",
          style = "position:absolute;left:10px;top:4px;"
        ),
        div(
          "SeekNSee",
          style = "text-align:center;font-size:26px;font-weight:bold;line-height:80px;"
        )
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      div(
        class = "sidebar-flex",
        div(
          class = "sidebar-body",
          tabsetPanel(
            tabPanel(
              "Selection",
              fileInput("reference_fasta", "Reference proteome FASTA/FAA"),
              radioButtons(
                "target_source",
                "Target source",
                choices = c(
                  "Reference proteome" = "reference_proteome",
                  "External target FASTA" = "external_fasta"
                ),
                selected = "reference_proteome"
              ),
              conditionalPanel(
                "input.target_source == 'reference_proteome'",
                textInput("reference_search", "Search target protein"),
                uiOutput("reference_target_select")
              ),
              conditionalPanel(
                "input.target_source == 'external_fasta'",
                fileInput("target_fasta", "External target FASTA/FAA"),
                radioButtons(
                  "target_scope",
                  "Process",
                  choices = c(
                    "Selected target" = "selected",
                    "All uploaded targets" = "all"
                  ),
                  selected = "selected"
                ),
                uiOutput("external_target_select"),
                selectInput(
                  "reference_policy",
                  "Reference specificity policy",
                  choices = c(
                    "Strict absence" = "strict_absence",
                    "Allow sequence-identical reference entries" = "allow_exact_target_match"
                  ),
                  selected = "strict_absence"
                )
              ),
              selectInput("enzyme", "Enzyme", choices = ENZYMES),
              selectInput(
                "missed_cleavages",
                "Maximum missed cleavages",
                choices = setNames(0:2, 0:2),
                selected = 0
              ),
              div(
                class = "fullwidth-btn",
                actionButton("run_analysis", "Verify", class = "btn-primary")
              )
            ),
            
            tabPanel(
              "Rule settings",
              selectInput("rule_unique_proteome", "Specific against reference proteome", RULE_CHOICES, "block"),
              selectInput("rule_unique_protein", "Unique accession group", RULE_CHOICES, "apply"),
              selectInput("rule_length", "Length", RULE_CHOICES, "block"),
              selectInput("rule_dibasic", "Dibasic (Trypsin and Trypsin/P)", RULE_CHOICES, "apply"),
              selectInput("rule_glyco", "Glycosylation", RULE_CHOICES, "apply"),
              selectInput("rule_chemical", "Chemical risk", RULE_CHOICES, "apply"),
              selectInput("rule_acidic", "D/E near tryptic cleavage", RULE_CHOICES, "apply"),
              selectInput("rule_terminal", "Terminal fragment", RULE_CHOICES, "apply"),
              selectInput("rule_miscleavage", "Missed cleavages", RULE_CHOICES, "apply"),
              selectInput(
                "rule_proline_p1prime",
                "Proline at P1-prime (Trypsin/P)",
                RULE_CHOICES,
                "apply"
              ),
              selectInput("rule_isobaric", "Isobaric I/L", RULE_CHOICES, "apply"),
              selectInput("rule_charge", "Heuristic charge > 2", RULE_CHOICES, "apply"),
              div(
                class = "fullwidth-btn",
                actionButton("update_rules", "Update", class = "btn-primary")
              )
            )
          )
        ),
        div(
          class = "sidebar-footer",
          div(
            class = "fullwidth-btn",
            downloadButton(
              "download_table",
              "Download Excel",
              class = "btn-primary",
              onclick = "Shiny.setInputValue('excel_download_clicked', Date.now(), {priority: 'event'});"
            )
          )
        )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        selected = "About",
        
        tabPanel(
          "About",
          br(),
          div(
            class = "info-page",
            h3("About SeekNSee"),
            p(
              "SeekNSee is a portable Shiny-based application for rule-driven ",
              "exploration, filtering, and prioritization of peptide candidates ",
              "derived from in silico enzymatic digestion."
            ),
            p(
              "The application is designed to help users identify, inspect, rank, ",
              "and export digestion-derived peptide candidates using a transparent ",
              "and configurable rule system."
            ),
            p(
              "SeekNSee separates the reference proteome from the target proteins. ",
              "A target may be selected directly from the uploaded reference proteome ",
              "or loaded independently from an external FASTA or FAA file."
            ),
            
            h4("What SeekNSee does"),
            tags$ul(
              tags$li("Loads a reference proteome from a FASTA or FAA file."),
              tags$li("Allows protein searches by protein ID, accession, gene name, or FASTA header text."),
              tags$li("Lets the user select a target protein from the uploaded reference proteome."),
              tags$li("Accepts a separate FASTA or FAA file containing one or multiple target proteins."),
              tags$li("Supports external targets that are absent from the reference proteome."),
              tags$li("Supports exploratory analysis of recombinant proteins, novel isoforms, sequence variants, synthetic constructs, and fusion proteins."),
              tags$li("Processes either one selected external target or all proteins in the uploaded target FASTA."),
              tags$li("Simulates enzymatic digestion using the selected enzyme, including Trypsin/P."),
              tags$li("Generates candidates with zero, one, or two missed cleavages."),
              tags$li("Generates peptide candidates between 4 and 25 amino acids."),
              tags$li("Records the start position, end position, and length of each candidate."),
              tags$li("Evaluates each candidate using 12 configurable selection rules."),
              tags$li("Allows every rule to be configured as ignore, apply, or block."),
              tags$li("Evaluates all rules even when a candidate fails one or more blocking rules."),
              tags$li("Ranks retained candidates according to their score."),
              tags$li("Records blocked candidates and all rules responsible for exclusion."),
              tags$li("Searches for candidate occurrences in the uploaded reference proteome."),
              tags$li("Exports a traceable Excel report containing run metadata, targets, retained and blocked candidates, rule settings, and rule-level results.")
            ),
            
            h4("Target protein sources"),
            p("SeekNSee supports two target-protein workflows:"),
            tags$ol(
              tags$li(
                tags$strong("Reference proteome: "),
                "the target is searched and selected directly from the uploaded reference proteome."
              ),
              tags$li(
                tags$strong("External target FASTA: "),
                "one or multiple targets are loaded from a separate FASTA or FAA file. External targets do not need to be present in the reference proteome."
              )
            ),
            
            h4("Reference-proteome specificity"),
            p(
              "The reference proteome is used to determine whether a candidate peptide ",
              "is shared with other protein entries. External targets provide two policies:"
            ),
            tags$ul(
              tags$li(
                tags$strong("Strict absence: "),
                "the peptide passes the specificity check only when it is absent from the entire reference proteome. This option is recommended for novel isoforms, fusion sequences, recombinant constructs, and synthetic proteins."
              ),
              tags$li(
                tags$strong("Allow sequence-identical reference entries: "),
                "a reference entry whose full sequence is identical to the external target is permitted. Matches to other entries can still cause the specificity rule to fail."
              )
            ),
            
            h4("What it is designed for"),
            p(
              "SeekNSee is intended for exploratory peptide candidate selection and ",
              "transparent decision support. It helps users inspect possible peptides ",
              "and prioritize them according to practical sequence-based criteria."
            ),
            p(
              "The configurable checks include specificity against the reference proteome, ",
              "normalized accession grouping, sequence length, dibasic motifs, potential ",
              "N-linked glycosylation motifs, refined chemical-risk motifs, D/E near valid tryptic cleavage sites, ",
              "terminal protein fragments, internal enzymatic cleavage sites, isobaric I/L ",
              "residues, and a heuristic charge estimate."
            ),
            
            h4("Quick workflow"),
            tags$ol(
              tags$li("Upload a reference proteome in FASTA or FAA format."),
              tags$li("Choose whether the target comes from the reference proteome or an external target FASTA."),
              tags$li("Search and select a reference target, or upload and select external targets."),
              tags$li("For an external target, choose the reference-specificity policy."),
              tags$li("Choose the digestion enzyme."),
              tags$li("Choose the maximum number of missed cleavages: 0, 1, or 2."),
              tags$li("Configure the rules as ignore, apply, or block."),
              tags$li("Click Verify to digest the targets and evaluate the candidates."),
              tags$li("Inspect retained peptides and their score indicators in the Peptides tab."),
              tags$li("Optionally modify the rules and click Update."),
              tags$li("Download the Excel report for candidate-level and rule-level traceability.")
            ),
            
            h4("Scoring and candidate status"),
            p(
              "Every candidate starts with a fixed maximum score of 12. A failed standard rule ",
              "configured as apply normally subtracts one point. For missed cleavages, ",
              "zero, one, and two missed cleavages contribute 0, -1, and -2 points, respectively. ",
              "A failed ignore rule is recorded but does not change the score."
            ),
            p(
              "A failed block rule marks the candidate as blocked. Blocked candidates are ",
              "excluded from the retained ranking, but all 12 rules are evaluated and ",
              "recorded in the Excel report. Candidate score and candidate status are ",
              "stored separately."
            ),
            
            h4("Important interpretation"),
            p(
              tags$strong("SeekNSee is an exploratory screening and decision-support tool.")
            ),
            p(
              "The application does not predict experimental peptide detectability and ",
              "does not replace empirical validation. The score is a transparent, ",
              "sequence-based prioritization framework rather than a validated physical, ",
              "chemical, or instrument-specific prediction model."
            ),
            p(
              "Results should be interpreted together with experimental design, sample ",
              "preparation, digestion conditions, analytical platform, acquisition method, ",
              "biological context, and any additional evidence available to the user."
            ),
            
            h4("Version information"),
            p(tags$strong("SeekNSee version: "), SEEKNSEE_VERSION),
            p(tags$strong("Algorithm version: "), ALGORITHM_VERSION)
          )
        ),
        
        tabPanel(
          "Rules",
          br(),
          div(
            class = "info-page",
            h3("Rules and Scoring"),
            p(
              "SeekNSee evaluates every generated peptide candidate using a fixed set ",
              "of 12 checks. Each check can be configured as ",
              tags$code("ignore"), ", ", tags$code("apply"), ", or ",
              tags$code("block"), "."
            ),
            
            h4("Rule states"),
            tags$ul(
              tags$li(
                tags$strong("ignore: "),
                "the result is recorded but does not affect filtering or scoring."
              ),
              tags$li(
                tags$strong("apply: "),
                "the candidate remains eligible for ranking. A failed apply rule normally subtracts 1 point; the missed-cleavage rule subtracts 1 or 2 points for one or two missed cleavages."
              ),
              tags$li(
                tags$strong("block: "),
                "the candidate is classified as blocked when the rule fails. All remaining rules are still evaluated and recorded."
              ),
              tags$li(
                tags$strong("not applicable: "),
                "the rule is not interpreted in that context and does not affect score or candidate status. D/E near tryptic cleavage applies only to Trypsin; dibasic applies to Trypsin and Trypsin/P; proline at P1-prime applies only to Trypsin/P and internal K/R cleavage boundaries."
              )
            ),
            
            h4("Checks currently evaluated"),
            tags$ul(
              tags$li(
                tags$strong("Specific against reference proteome: "),
                "evaluates whether the peptide is specific relative to the uploaded reference proteome. Interpretation depends on the target source and specificity policy."
              ),
              tags$li(
                tags$strong("Unique accession group: "),
                "passes when all matching reference entries belong to no more than one normalized accession group. This is not a complete isoform-specificity test."
              ),
              tags$li(
                tags$strong("Length: "),
                "passes when peptide length is between 8 and 25 amino acids, inclusive."
              ),
              tags$li(
                tags$strong("Dibasic: "),
                "flags peptides containing ", tags$code("KK"), ", ",
                tags$code("KR"), ", ", tags$code("RK"), ", or ",
                tags$code("RR"), "."
              ),
              tags$li(
                tags$strong("Glycosylation: "),
                "flags potential N-linked glycosylation motifs using ",
                tags$code("N-X-S/T"), ", where X is not proline."
              ),
              tags$li(
                tags$strong("Chemical risk: "),
                "flags peptides containing M, W, C, NG, QG, DP, or beginning with Q/E because of N-terminal pyroglutamate risk."
              ),
              tags$li(
                tags$strong("D/E near tryptic cleavage: "),
                "flags D or E within the P3-P3' window of a valid tryptic cleavage site delimiting the peptide. K/R-P is excluded. This rule is not applicable to non-trypsin enzymes."
              ),
              tags$li(
                tags$strong("Terminal fragment: "),
                "flags candidates corresponding to the N-terminal or C-terminal fragment of the target protein."
              ),
              tags$li(
                tags$strong("Missed cleavages: "),
                "records the number of missed cleavages selected during digestion. With apply, one or two missed cleavages subtract one or two points, respectively."
              ),
              tags$li(
                tags$strong("Proline at P1-prime: "),
                "for Trypsin/P, flags candidates created at internal ",
                tags$code("K|P"), " or ", tags$code("R|P"),
                " boundaries. The inverse motifs P|K and P|R are not P1-prime proline."
              ),
              tags$li(
                tags$strong("Isobaric I/L: "),
                "flags peptides containing isoleucine ", tags$code("I"),
                " or leucine ", tags$code("L"), "."
              ),
              tags$li(
                tags$strong("Heuristic charge > 2: "),
                "uses ", tags$code("z_hat = 1 + nK + nR + nH"),
                " and flags candidates with an estimated value greater than 2."
              )
            ),
            
            h4("Reference-proteome specificity"),
            tags$ul(
              tags$li(
                tags$strong("Target selected from the reference proteome: "),
                "the peptide may occur in the selected entry, but fails specificity when it is also found in another reference entry."
              ),
              tags$li(
                tags$strong("External target with strict absence: "),
                "the peptide passes only when absent from the entire reference proteome."
              ),
              tags$li(
                tags$strong("External target with an identical reference entry allowed: "),
                "a reference entry with a complete sequence identical to the external target is permitted. Matches to other entries still cause failure."
              )
            ),
            p(
              tags$strong("Recommended use: "),
              "strict absence is generally the most appropriate policy when searching for peptides that distinguish a novel isoform, recombinant construct, fusion protein, synthetic sequence, or sequence variant."
            ),
            
            h4("Cleavage logic by enzyme"),
            tags$ul(
              tags$li(tags$strong("Trypsin: "), "cleaves after K or R, except when followed by P."),
              tags$li(tags$strong("Trypsin/P: "), "cleaves after K or R, including K-P and R-P bonds."),
              tags$li(tags$strong("LysC: "), "internal K sites."),
              tags$li(tags$strong("ArgC: "), "internal R sites."),
              tags$li(tags$strong("AspN: "), "internal D positions because AspN cleaves before D."),
              tags$li(tags$strong("GluC: "), "internal E sites."),
              tags$li(tags$strong("Chymotrypsin: "), "internal F, W, or Y sites, except when followed by P."),
              tags$li(tags$strong("Elastase: "), "internal A, V, L, I, S, or G sites."),
              tags$li(tags$strong("ProAlanase: "), "internal P or A sites.")
            ),
            
            h4("Scoring model"),
            p(
              "Every candidate starts with a fixed maximum score of ",
              tags$strong("12"), "."
            ),
            tags$ul(
              tags$li("A failed ", tags$code("apply"), " rule usually subtracts 1 point."),
              tags$li("The missed-cleavage rule subtracts 1 point for one missed cleavage and 2 points for two."),
              tags$li("A failed ", tags$code("ignore"), " rule is recorded but does not change the score."),
              tags$li("A failed ", tags$code("block"), " rule changes candidate status to ", tags$code("blocked"), ".")
            ),
            p(
              style = "font-family:monospace;font-size:16px;margin-left:25px;",
              tags$strong("peptide score = 12 + sum of rule_score")
            ),
            p(
              "The score starts at 12 and is calculated as 12 plus the sum of rule_score values. ",
              "Standard failed apply rules contribute -1; the missed-cleavage rule contributes 0, -1, or -2."
            ),
            p(
              "Candidate score and candidate status are stored separately. A blocked ",
              "candidate keeps a numeric score but is excluded from the retained ranking."
            ),
            p(
              tags$strong("Important: "),
              "the maximum remains 12 even when rules are configured as ",
              tags$code("ignore"), ". Comparisons between runs must also consider rule configuration."
            ),
            
            h4("Complete rule evaluation"),
            p(
              "All 12 rules are evaluated before final status is assigned. Evaluation ",
              "does not stop at the first failed blocking rule. Multiple reasons are ",
              "recorded in the ", tags$code("blocked_by"), " field."
            ),
            p(
              style = "font-family:monospace;font-size:15px;margin-left:25px;",
              "blocked_by = length;unique_proteome"
            ),
            
            h4("Excel export"),
            tags$ul(
              tags$li(tags$strong("Run_summary: "), "run, input, version, enzyme, policy, and candidate totals."),
              tags$li(tags$strong("Rule_configuration: "), "the state of every rule."),
              tags$li(tags$strong("Target_proteins: "), "one row per processed target and its candidate counts."),
              tags$li(tags$strong("Ranked_peptides: "), "one row per retained candidate."),
              tags$li(tags$strong("Rule_trace: "), "one row per candidate-rule combination."),
              tags$li(tags$strong("Blocked_peptides: "), "one row per blocked candidate and its exclusion reasons."),
              tags$li(tags$strong("Rule_definitions: "), "definitions associated with the current algorithm version.")
            ),
            
            h4("Interpretation"),
            p(tags$strong("SeekNSee is an exploratory screening and decision-support tool.")),
            p(
              "The rules and scores form a transparent sequence-based prioritization ",
              "framework. They are not a validated prediction of experimental peptide ",
              "detectability, ionization, chromatographic behavior, or instrument performance."
            )
          )
        ),
        
        tabPanel(
          "Peptides",
          br(),
          uiOutput("run_status"),
          uiOutput("peptide_display")
        )
      )
    )
  )
)

active_sessions <- 0L

server <- function(input, output, session) {
  active_sessions <<- active_sessions + 1L
  session$onSessionEnded(function() {
    active_sessions <<- max(0L, active_sessions - 1L)
    if (active_sessions == 0L) stopApp()
  })
  
  reference_proteome <- reactiveVal(NULL)
  reference_filename <- reactiveVal("")
  external_targets <- reactiveVal(NULL)
  target_filename <- reactiveVal("")
  analysis_results <- reactiveVal(NULL)
  analysis_context <- reactiveVal(NULL)
  reference_match_cache <- new.env(parent = emptyenv())
  
  observeEvent(input$reference_fasta, {
    req(input$reference_fasta$datapath)
    
    progress <- Progress$new(
      session = session,
      min = 0,
      max = 1
    )
    
    progress_closed <- FALSE
    
    close_progress <- function() {
      if (!progress_closed) {
        progress$close()
        progress_closed <<- TRUE
      }
      
      invisible(NULL)
    }
    
    progress$set(
      value = 0,
      message = "Processing reference proteome...",
      detail = "Reading and validating the FASTA file."
    )
    
    on.exit(
      close_progress(),
      add = TRUE
    )
    
    tryCatch({
      parsed_reference <- read_fasta(
        path = input$reference_fasta$datapath,
        
        progress_callback = function(
    completed,
    total
        ) {
          proportion <- if (total > 0L) {
            completed / total
          } else {
            0
          }
          
          progress$set(
            value = proportion,
            message = "Processing reference proteome...",
            detail = paste0(
              "Reading and validating protein ",
              format(
                completed,
                big.mark = ",",
                scientific = FALSE
              ),
              " of ",
              format(
                total,
                big.mark = ",",
                scientific = FALSE
              ),
              "."
            )
          )
        }
      )
      
      progress$set(
        value = 0.98,
        message = "Finalizing reference proteome...",
        detail = "Preparing protein search and clearing previous results."
      )
      
      reference_proteome(
        parsed_reference
      )
      
      reference_filename(
        input$reference_fasta$name
      )
      
      cache_keys <- ls(
        envir = reference_match_cache,
        all.names = TRUE
      )
      
      if (length(cache_keys) > 0L) {
        rm(
          list = cache_keys,
          envir = reference_match_cache
        )
      }
      
      analysis_results(NULL)
      analysis_context(NULL)
      
      progress$set(
        value = 1,
        message = "Reference proteome ready.",
        detail = paste0(
          format(
            nrow(parsed_reference),
            big.mark = ",",
            scientific = FALSE
          ),
          " proteins loaded."
        )
      )
      
      close_progress()
      
      showNotification(
        paste0(
          "Reference proteome loaded: ",
          format(
            nrow(parsed_reference),
            big.mark = ",",
            scientific = FALSE
          ),
          " proteins."
        ),
        type = "message"
      )
    }, error = function(e) {
      close_progress()
      
      showNotification(
        conditionMessage(e),
        type = "error",
        duration = NULL
      )
    })
  })
  
  observeEvent(input$target_fasta, {
    req(input$target_fasta$datapath)
    tryCatch({
      external_targets(read_fasta(input$target_fasta$datapath))
      target_filename(input$target_fasta$name)
      analysis_results(NULL)
      showNotification("External target FASTA loaded.", type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })
  
  output$reference_target_select <- renderUI({
    df <- reference_proteome()
    if (is.null(df)) return(NULL)
    
    q <- trimws(if (is.null(input$reference_search)) "" else input$reference_search)
    if (nchar(q) < 3) return(helpText("Enter at least three characters."))
    
    hit <- grepl(q, df$id, ignore.case = TRUE) |
      grepl(q, df$gene, ignore.case = TRUE) |
      grepl(q, df$name, ignore.case = TRUE)
    matches <- df[hit, , drop = FALSE]
    
    if (!nrow(matches)) return(div("No matches found"))
    
    labels <- ifelse(
      nzchar(matches$gene),
      paste(matches$id, matches$gene, sep = " | "),
      matches$id
    )
    
    selectInput(
      "selected_reference_target",
      "Protein",
      choices = setNames(matches$id, labels)
    )
  })
  
  output$external_target_select <- renderUI({
    df <- external_targets()
    if (is.null(df) || identical(input$target_scope, "all")) return(NULL)
    
    labels <- ifelse(
      nzchar(df$gene),
      paste(df$id, df$gene, sep = " | "),
      df$id
    )
    
    selectInput(
      "selected_external_target",
      "Target protein",
      choices = setNames(df$id, labels)
    )
  })
  
  selected_targets <- reactive({
    if (identical(input$target_source, "reference_proteome")) {
      req(reference_proteome(), input$selected_reference_target)
      reference_proteome() %>% filter(id == input$selected_reference_target)
    } else {
      req(external_targets())
      if (identical(input$target_scope, "all")) {
        external_targets()
      } else {
        req(input$selected_external_target)
        external_targets() %>% filter(id == input$selected_external_target)
      }
    }
  })
recompute <- function() {
    req(reference_proteome())
    targets <- selected_targets()
    source_label <- input$target_source

    policy <- if (source_label == "reference_proteome") {
      "allow_selected_reference_target"
    } else {
      input$reference_policy
    }

    rules <- build_rules(input)
    missed_cleavages <- suppressWarnings(as.integer(input$missed_cleavages))

    if (
      length(missed_cleavages) != 1L ||
      is.na(missed_cleavages) ||
      missed_cleavages < 0L ||
      missed_cleavages > 2L
    ) {
      stop(
        "Maximum missed cleavages must be 0, 1, or 2.",
        call. = FALSE
      )
    }

    withProgress(message = "Digesting and evaluating targets...", value = 0.2, {
      result <- make_target_candidates(
        targets = targets,
        source_label = source_label,
        enzyme = input$enzyme,
        rules = rules,
        reference_proteome = reference_proteome(),
        reference_policy = policy,
        missed_cleavages = missed_cleavages,
        match_cache = reference_match_cache
      )
      incProgress(0.8)
    })

    analysis_results(result)
    analysis_context(list(
      run_datetime = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      reference_filename = reference_filename(),
      target_filename = if (source_label == "external_fasta") {
        target_filename()
      } else {
        reference_filename()
      },
      target_source = source_label,
      reference_policy = policy,
      enzyme = input$enzyme,
      missed_cleavages = missed_cleavages,
      rules = rules
    ))
  }
  
  observeEvent(input$run_analysis, {
    tryCatch(
      recompute(),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      }
    )
  })
  
  observeEvent(input$update_rules, {
    req(analysis_results())
    tryCatch(
      recompute(),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      }
    )
  })
  
  output$run_status <- renderUI({
    result <- analysis_results()
    if (is.null(result)) {
      return(div(
        class = "status-box",
        "Load a reference proteome, select the target source, and click Verify."
      ))
    }
    
    candidates <- result$candidates
    retained <- if (nrow(candidates)) {
      sum(candidates$candidate_status == "retained")
    } else {
      0
    }
    blocked <- if (nrow(candidates)) {
      sum(candidates$candidate_status == "blocked")
    } else {
      0
    }
    
    div(
      class = "status-box",
      paste0(
        nrow(result$target_summary), " target(s); ",
        nrow(candidates), " candidates; ",
        retained, " retained; ",
        blocked, " blocked."
      )
    )
  })
  
  output$peptide_display <- renderUI({
    result <- analysis_results()
    req(result)
    
    candidates <- result$candidates
    if (is.null(candidates) || nrow(candidates) == 0) {
      return(div(class = "status-box", "No peptide candidates were generated."))
    }
    
    retained <- candidates %>%
      filter(candidate_status == "retained") %>%
      arrange(target_id, desc(peptide_score), peptide, start)
    
    if (nrow(retained) == 0) {
      return(div(
        class = "status-box",
        "No retained peptide candidates were found with the current rule settings."
      ))
    }
    
    target_groups <- split(retained, retained$target_id, drop = TRUE)
    
    tagList(lapply(names(target_groups), function(target_name) {
      target_data <- target_groups[[target_name]]
      scores <- as.numeric(
        target_data$peptide_score
      )
      
      score_minimum <- 0
      score_maximum <- length(RULE_NAMES)
      
      normalized_scores <- (
        scores - score_minimum
      ) / (
        score_maximum - score_minimum
      )
      
      normalized_scores[
        !is.finite(normalized_scores)
      ] <- 0
      
      normalized_scores <- pmax(
        0,
        pmin(
          1,
          normalized_scores
        )
      )
      
      items_per_column <- 10L
      peptide_chunks <- split(
        seq_len(nrow(target_data)),
        ceiling(seq_len(nrow(target_data)) / items_per_column)
      )
      
      target_heading <- NULL
      if (length(target_groups) > 1) {
        target_heading <- div(
          style = paste0(
            "font-size:16px;font-weight:bold;margin-top:14px;",
            "margin-bottom:12px;padding-bottom:5px;",
            "border-bottom:1px solid #dddddd;"
          ),
          target_name
        )
      }
      
      tagList(
        target_heading,
        fluidRow(lapply(peptide_chunks, function(chunk_indices) {
          column_width <- if (length(peptide_chunks) <= 2) {
            3
          } else if (length(peptide_chunks) <= 3) {
            4
          } else {
            3
          }
          
          column(
            width = column_width,
            class = "peptide-column",
            lapply(chunk_indices, function(i) {
              peptide <- target_data$peptide[i]
              score <- target_data$peptide_score[i]
              start_position <- target_data$start[i]
              end_position <- target_data$end[i]
              reference_count <- target_data$reference_entry_count[i]
              marker_position <- max(0, min(100, normalized_scores[i] * 100))
              
              tooltip <- paste0(
                "Score: ", score,
                " | Position: ", start_position, "-", end_position,
                " | Missed cleavages: ", target_data$missed_cleavages[i],
                " | Reference entries: ", reference_count
              )
              
              div(
                style = "margin-bottom:14px;font-family:monospace;max-width:230px;",
                title = tooltip,
                div(
                  style = paste0(
                    "font-size:15px;",
                    "white-space:nowrap;",
                    "overflow:hidden;",
                    "text-overflow:ellipsis;",
                    "margin-bottom:3px;"
                  ),
                  
                  peptide,
                  
                  span(
                    paste0(
                      "  [",
                      score,
                      "/",
                      length(RULE_NAMES),
                      "]"
                    ),
                    
                    style = paste0(
                      "font-family:Arial,sans-serif;",
                      "font-size:12px;",
                      "color:#666666;",
                      "margin-left:5px;"
                    )
                  )
                ),
                div(
                  style = paste0(
                    "width:100px;height:6px;",
                    "background:linear-gradient(to right,",
                    "#ff1f1f 0%,#ff8c00 25%,#ffff00 50%,",
                    "#7ecb20 75%,#008000 100%);",
                    "position:relative;"
                  ),
                  div(
                    style = paste0(
                      "position:absolute;left:", marker_position, "%;top:0;",
                      "width:3px;height:6px;background:#000000;",
                      "transform:translateX(-1px);"
                    )
                  )
                )
              )
            })
          )
        }))
      )
    }))
  })
  
  observeEvent(input$excel_download_clicked, {
    result <- analysis_results()
    context <- analysis_context()
    
    if (is.null(result) || is.null(context)) {
      showNotification(
        "Run an analysis before downloading the Excel report.",
        type = "warning",
        duration = NULL
      )
      return(invisible(NULL))
    }
    
    tryCatch({
      workbook <- build_seeknsee_workbook(
        result = result,
        context = context,
        reference_proteome = reference_proteome()
      )
      
      adapted_export <- adapt_workbook_for_excel(workbook)
      
      if (isTRUE(adapted_export$truncation_used)) {
        showNotification(
          paste0(
            "Excel report adapted successfully. ",
            adapted_export$abbreviated_cells,
            if (adapted_export$abbreviated_cells == 1L) {
              " reference-ID cell was abbreviated; "
            } else {
              " reference-ID cells were abbreviated; "
            },
            format(
              adapted_export$omitted_identifier_occurrences,
              big.mark = ",",
              scientific = FALSE
            ),
            " identifier occurrences were omitted only from the displayed Excel lists. ",
            "All identifiers were used in the scientific evaluation, and no rule, score, status, or match count was changed."
          ),
          type = "warning",
          duration = NULL
        )
      } else {
        showNotification(
          "Excel report validated. Download starting...",
          type = "message",
          duration = 4
        )
      }
    }, error = function(e) {
      showNotification(
        paste0(
          "Excel export validation failed: ",
          conditionMessage(e)
        ),
        type = "error",
        duration = NULL
      )
    })
  })
  
  output$download_table <- downloadHandler(
    filename = function() {
      result <- analysis_results()
      label <- if (!is.null(result) && nrow(result$target_summary) == 1) {
        result$target_summary$target_id[1]
      } else {
        "target_panel"
      }
      make_export_filename(label)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    content = function(file) {
      result <- analysis_results()
      context <- analysis_context()
      req(result, context)
      
      workbook <- build_seeknsee_workbook(
        result = result,
        context = context,
        reference_proteome = reference_proteome()
      )
      
      adapted_export <- adapt_workbook_for_excel(workbook)
      
      write_xlsx(
        adapted_export$workbook,
        path = file
      )
      
    }
  )
  
}

shinyApp(ui, server)
