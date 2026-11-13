###############################################################################
# Analysis:
# DDA peptide score benchmarking
#
# Purpose:
# Process MaxQuant evidence output, generate peptide-level summary
# tables, integrate  scores, perform score binning, and
# calculate metric associations used in downstream analyses.
#
# Inputs:
#   maxquant_evidence.txt (Maxquant evidence)
#   scores.parquet (SeekNSee R script output)
#
# Outputs:
#   score_binned_features.parquet
#   score_binned_features_raw.parquet
#   analytical_features.parquet
#   crossref.parquet
#   correlation_results.tsv
###############################################################################

###############################################################################
# Load libraries
###############################################################################

library(arrow)
library(dplyr)

###############################################################################
# Import data
###############################################################################

evidence <- read.table("maxquant_evidence.txt", header = TRUE, sep = "\t")

evidence <- evidence %>%
    filter(Decoy != "+", Potential.contaminant != "+")

scores <- read_parquet("scores.parquet")

###############################################################################
# Generate peptide-level DDA summary table
# Retain the most abundant PSM for each precursor in each run
# Retain only precursors detected in all runs
# Average precursor-level metrics across runs
# Collapse charge states by keeping the most abundant precursor
###############################################################################

n_runs <- 3

evidence_best_psm <- evidence %>%
    group_by(Raw.file, Sequence, Charge) %>%
    slice_max(order_by = Intensity, n = 1, with_ties = FALSE) %>%
    ungroup()

evidence_complete <- evidence_best_psm %>%
    group_by(Sequence, Charge) %>%
    filter(n_distinct(Raw.file) == n_runs) %>%
    ungroup()

evidence_averaged <- evidence_complete %>%
    group_by(Sequence, Charge) %>%
    summarise(across(
            where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            across(where(is.character), first),
            .groups = "drop")

evidence_collapsed <- evidence_averaged %>%
    group_by(Sequence) %>%
    slice_max(order_by = Intensity, n = 1, with_ties = FALSE) %>%
    ungroup()

###############################################################################
# Select features of interest
###############################################################################

analytical_features <- evidence_collapsed %>%
    select(
        Sequence,
        m.z,
        Retention.length,
        Number.of.scans,
        Number.of.isotopic.peaks,
        Mass,
        PEP,
        Score,
        Delta.score,
        Intensity,
        MS.MS.scan.number,
        MS.MS.count,
        MS.MS.m.z,
        Best.MS.MS,
        Retention.time,
        Retention.length)

###############################################################################
# Transform quantitative features
###############################################################################

analytical_features <- analytical_features %>%
    mutate(PEP = ifelse(PEP == 0, 1e-300, PEP))

log_transformed_features <- analytical_features %>%
    mutate(Intensity = log10(Intensity + 1), PEP = -log10(PEP))

###############################################################################
# Integrate  scores with DDA metrics
###############################################################################

crossref <- scores %>%
    left_join(log_transformed_features, by = c("peptide" = "Sequence"))

crossref_dataset <- crossref %>%
    filter(!is.na(Mass)) %>%
    rename(Length = length)

crossref_unlog <- scores %>%
    left_join(analytical_features, by = c("peptide" = "Sequence"))

crossref_dataset_unlog <- crossref_unlog %>%
    filter(!is.na(Mass)) %>%
    rename(Length = length)

###############################################################################
# Assign  score bins
###############################################################################

create_score_bins <- function(data) {
    
    data %>%
        mutate(
            score_bin = case_when(
                peptide_score < 9 ~ "5-8",
                peptide_score == 9 ~ "9",
                peptide_score == 10 ~ "10",
                peptide_score == 11 ~ "11",
                peptide_score == 12 ~ "12"
            ),
            score_bin = factor(score_bin,
                levels = c("5-8", "9", "10", "11", "12")))}

score_binned_features <- create_score_bins(crossref_dataset)

score_binned_features_raw <- create_score_bins(crossref_dataset_unlog)

###############################################################################
# Define metrics for downstream analyses
###############################################################################

columns <- setdiff(names(score_binned_features), c("peptide_score", "score_bin", "peptide"))

###############################################################################
# Evaluate associations between  scores and DDA features
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
write_parquet(score_binned_features, "DDA_score_binned_features.parquet")
write_parquet(score_binned_features_raw, "DDA_score_binned_features_raw.parquet")
write_parquet(analytical_features , "DDA_analytical_features.parquet")
write_parquet(crossref, "DDA_crossref.parquet")
write_parquet(crossref_dataset, "DDA_crossref_dataset.parquet")
write.table(correlation_results, "DDA_correlation_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
