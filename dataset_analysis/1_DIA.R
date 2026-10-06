###############################################################################
# Analysis:
# DIA peptide score benchmarking
#
# Purpose:
# Process DIA-NN output, generate summary tables,
# integrate scores and calculate metric associations used in downstream analyses.
#
# Inputs:
#   diann_report (DIA-NN parquet)
#   result_peptide_scores.parquet (Tool R script output)
#
# Outputs:
#   score_binned_features.parquet
#   score_binned_features_raw.parquet
#   analytical_features.parquet
#   crossref.parquet
#   correlation_results.tsv
#   scores.parquet
#
###############################################################################

###############################################################################
# Load libraries
###############################################################################

library(arrow)
library(dplyr)
library(tidyr)

###############################################################################
# Import DIA-NN data
###############################################################################

report <- read_parquet("diann_report.parquet") %>%
    filter(!grepl("crap", Protein.Group, ignore.case = TRUE))

candidates <- read_parquet("result_peptide_scores.parquet")

scores <- candidates %>% 
    filter(candidate_status != "blocked") %>% 
    distinct(peptide, .keep_all = TRUE)

###############################################################################
# Generate peptide-level DIA summary table
#   Filtering ID: 
#   - (1) Only found in all samples
#   - (2) Average by precursor ID
#   - (3) Keep only the most abundant per Stripped sequence (e.g. collapse charges)
###############################################################################

complete_precursors <- report %>%
    group_by(Precursor.Id) %>%
    filter(n_distinct(Run.Index) == 5) %>%
    ungroup()

averaged_precursors <- complete_precursors %>%
    group_by(Precursor.Id) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
        across(where(is.character), first),
        .groups = "drop")

dominant_sequences <- averaged_precursors %>%
    group_by(Stripped.Sequence) %>%
    slice_max(order_by = Precursor.Quantity, n = 1, with_ties = FALSE) %>%
    ungroup()

dominant_sequences <- dominant_sequences %>%
    mutate(Peak.width = RT.Stop - RT.Start)

###############################################################################
# Select DIA metrics
###############################################################################

analytical_features <- dominant_sequences %>%
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
# Integrate  scores with DIA metrics
###############################################################################

crossref <- scores  %>%
    left_join(log_transformed_features, by = c("peptide" = "Stripped.Sequence"))

crossref_dataset <- crossref %>%
    filter(!is.na(Peak.width))

crossref_unlog <- scores  %>%
    left_join(analytical_features, by = c("peptide" = "Stripped.Sequence"))

crossref_dataset_unlog <- crossref_unlog %>%
    filter(!is.na(Peak.width))

###############################################################################
# Assign  score bins
###############################################################################

create_score_bins <- function(data) {
    
    data %>%
        mutate(
            score_bin = case_when(
                peptide_score < 8  ~ "5-7",
                peptide_score == 8 ~ "8",
                peptide_score == 9 ~ "9",
                peptide_score == 10 ~ "10",
                peptide_score == 11 ~ "11",
                peptide_score == 12 ~ "12"),
            score_bin = factor(score_bin,
                levels = c("5-7", "8", "9",  "10", "11",  "12"))) %>%
                rename(Length = length)
}

score_binned_features <- create_score_bins(crossref_dataset)

score_binned_features_raw <- create_score_bins(crossref_dataset_unlog)

###############################################################################
# Define metrics for downstream analyses
###############################################################################

columns <- setdiff(names(score_binned_features), c("peptide_score","score_bin","peptide"))

quality_metrics <- c(
    "Length",
    "FWHM",
    "Peak.width",
    "Empirical.Quality",
    "Quantity.Quality",
    "Evidence",
    "Q.Value",
    "PEP")

###############################################################################
# Evaluate associations between  scores and DIA features
###############################################################################

lm_results <- lapply(columns,function(metric) {          
            fit <- lm(reformulate(metric, response = "peptide_score"),
            data = score_binned_features)
            
            fit_summary <- summary(fit)
            
            data.frame(metric = metric,
            r_squared = fit_summary$r.squared,
            adj_r_squared = fit_summary$adj.r.squared,
            slope = coef(fit)[2],
            p_value = coef(fit_summary)[2, 4])
}) %>%
    bind_rows() %>%
    arrange(desc(r_squared))


cor_results <- lapply(columns, function(metric) {
            data.frame(metric = metric,
            rho = cor(
                score_binned_features$peptide_score,
                score_binned_features[[metric]],
                method = "spearman",
                use = "complete.obs"))}) %>%
    bind_rows() %>%
    arrange(desc(abs(rho)))


correlation_results <- lm_results %>%
    left_join(cor_results, by = "metric") %>%
    arrange(desc(abs(rho)))

print(correlation_results)

###############################################################################
# Pairwise correlation analyses LC-MS
###############################################################################

comparisons <- list(
    "Peak.width vs FWHM"            = c("Peak.width", "FWHM"),
    "Ms1.Area vs Signal.Before"     = c("Ms1.Area", "Ms1.Total.Signal.Before"),
    "FWHM vs Ms1.Area"              = c("FWHM", "Ms1.Area"),
    "FWHM vs Precursor.Quantity"    = c("FWHM", "Precursor.Quantity")
)

corr_metrics_LCMS <- lapply(names(comparisons), function(name) {
    x <- comparisons[[name]][1]
    y <- comparisons[[name]][2]
       spearman_test <- cor.test(
        log_transformed_features[[x]],
        log_transformed_features[[y]],
        method = "spearman",
        exact = FALSE)
    
    rho <- unname(spearman_test$estimate)
    
    n <- sum(complete.cases(
            log_transformed_features[[x]],
            log_transformed_features[[y]]))
    
    z <- atanh(rho)
    se <- 1 / sqrt(n - 3)
    
    data.frame(
        comparison = name,
        rho = rho,
        spearman_p = spearman_test$p.value,
        ci_low = tanh(z - 1.96 * se),
        ci_high = tanh(z + 1.96 * se),
        n = n)}) %>%
    bind_rows() %>%
    arrange(desc(abs(rho)))

print(corr_metrics)

###############################################################################
# Export processed datasets
###############################################################################
write_parquet(score_binned_features, "DIA_score_binned_features.parquet")
write_parquet(score_binned_features_raw, "DIA_score_binned_features_raw.parquet")
write_parquet(analytical_features, "DIA_analytical_features.parquet")
write_parquet(crossref, "DIA_crossref.parquet")
write_parquet(crossref_dataset, "DIA_crossref_dataset.parquet")
write_parquet(scores,"scores.parquet")
write.table(corr, "DIA_correlation_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
