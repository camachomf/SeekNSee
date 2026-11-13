###############################################################################
# Analysis:
# DIA peptide score benchmarking (fractionated dataset)
#
# Purpose:
# Process DIA-NN output from fractionated acquisitions, generate
# summary tables, integrate  scores, perform
# score binning, and calculate metric associations used in
# downstream analyses.
#
# Inputs:
#   scores.parquet (SeekNSee R script output)
#   diann_report.parquet (DIA-NN report)
#
# Outputs:
#   fractionated_score_binned_features_raw.parquet
#   fractionated_score_binned_features.parquet
#   fractionated_correlation_results.tsv
#
###############################################################################

###############################################################################
# Load libraries
###############################################################################

library(arrow)
library(dplyr)
library(stringr)

###############################################################################
# Import data
###############################################################################

scores <- read_parquet("scores.parquet")

report <- read_parquet("diann_report.parquet") %>%
    filter(!grepl("crap", Protein.Group, ignore.case = TRUE))

###############################################################################
# Parse replicate and fraction information
###############################################################################

report <- report %>%
    mutate(Replicate = str_extract(Run, "(?<=rep)\\d+") |> as.integer(),
        Fraction = str_extract(Run, "(?<=Frac)\\d+") |> as.integer())

###############################################################################
# Generate precursor consensus table
# Retain only precursors detected in all triplicates
# Average across Fraction
###############################################################################

precursor_consensus <- report_replicates %>%
    group_by(Precursor.Id, Fraction) %>%
    filter(n_distinct(Replicate) == 3) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
        across(where(is.character), first), .groups = "drop")

###############################################################################
# Select most abundant fraction per precursor
###############################################################################

report_best_fraction <- precursor_consensus %>%
    group_by(Precursor.Id) %>%
    slice_max(order_by = Precursor.Quantity, n = 1, with_ties = FALSE) %>%
    ungroup()

###############################################################################
# Collapse charge states
###############################################################################

report_collapsed <- report_best_fraction %>%
    group_by(Stripped.Sequence) %>%
    slice_max(order_by = Precursor.Quantity, n = 1, with_ties = FALSE) %>%
    ungroup()

###############################################################################
# Select DIA features
###############################################################################

analytical_features <- report_collapsed %>%
    mutate(Peak.width = RT.Stop - RT.Start) %>%
    select(
        Stripped.Sequence,
        Precursor.Mz,
        Precursor.Charge,
        RT.Start,
        Precursor.Quantity,
        Precursor.Normalised,
        Ms1.Area,
        Ms1.Normalised,
        Ms1.Total.Signal.Before,
        FWHM,
        Peak.width,
        Empirical.Quality,
        Quantity.Quality,
        Evidence,
        Q.Value,
        PEP)

###############################################################################
# Log-transform quantitative metrics
###############################################################################

log_transformed_features <- analytical_features %>%
    mutate(
        Precursor.Quantity      = log10(Precursor.Quantity + 1),
        Precursor.Normalised    = log10(Precursor.Normalised + 1),
        Ms1.Total.Signal.Before = log10(Ms1.Total.Signal.Before + 1),
        Ms1.Area                = log10(Ms1.Area + 1),
        Ms1.Normalised          = log10(Ms1.Normalised + 1),
        Q.Value                 = -log10(Q.Value),
        PEP                     = -log10(PEP))

###############################################################################
# Integrate scores
###############################################################################

crossref <- scores %>%
    left_join(log_transformed_features, by = c("peptide" = "Stripped.Sequence"))

crossref_dataset <- result %>%
    filter(!is.na(Peak.width))

crossref_unlog <- scores %>%
    left_join(analytical_features, by = c("peptide" = "Stripped.Sequence"))

crossref_dataset_unlog <- result_nolog %>%
    filter(!is.na(Peak.width))

###############################################################################
# Assign score bins
###############################################################################

create_score_bins <- function(data) {
    data %>%
        mutate(
            score_bin = case_when(
                peptide_score < 7 ~ "5-6",
                peptide_score == 7 ~ "7",
                peptide_score == 8 ~ "8",
                peptide_score == 9 ~ "9",
                peptide_score == 10 ~ "10",
                peptide_score == 11 ~ "11",
                peptide_score == 12 ~ "12"),
            score_bin = factor(score_bin,
                levels = c("5-6", "7", "8", "9", "10", "11", "12"))) %>% 
        rename(Length = length)}

score_binned_features <- create_score_bins(crossref_dataset)

score_binned_features_raw <- create_score_bins(crossref_dataset_unlog)

###############################################################################
# Define metrics for downstream analyses
###############################################################################

columns <- setdiff(names(score_binned_features),
    c("peptide_score", "score_bin", "peptide"))

###############################################################################
# Evaluate associations between scores and DIA features
###############################################################################

correlation_results <- lapply(columns, function(metric) {
  
  corr <- cor.test(
    score_binned_features$peptide_score,
    score_binned_features[[metric]],
    method = "spearman",
    exact = FALSE)
  
  data.frame(
    metric = metric,
    rho = unname(corr$estimate),
    p_value = corr$p.value)}) %>%
  bind_rows() %>%
  arrange(desc(abs(rho)))


###############################################################################
# Export processed datasets
###############################################################################

write_parquet(score_binned_features, "fractionated_score_binned_features.parquet")
write_parquet(score_binned_features_raw, "fractionated_score_binned_features_raw.parquet")
write_parquet(analytical_features, "fractionated_analytical_features.parquet")
write_parquet(crossref, "fractionated_crossref.parquet")
write.table(correlation_results, "fractionated_correlation_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

