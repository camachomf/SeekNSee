###############################################################################
# Analysis:
# PeptideRanger benchmarking 
#
# Purpose:
# Evaluate associations between PeptideRanger RF scores and
# analytical features from acquisitions.
#
# Inputs:
#   proteomicsDB_detection.parquet
#
# Outputs:
#   pDB_score_binned_raw
#   pDB_score_binned_log
#   correlation_results.tsv
###############################################################################

###############################################################################
# Load libraries
###############################################################################

library(arrow)
library(dplyr)

###############################################################################
# Import data
# Alternatively, use integrated DDA data
###############################################################################

proteomicsDB_detection <- read_parquet("DIA_proteomicsDB_detection.parquet") %>%
    mutate(Length = nchar(sequence)) %>%
    filter(!is.na(Precursor.Mz)) %>%
    select(-detected)

###############################################################################
# Assign RF score bins
###############################################################################

pDB_score_binned <- proteomicsDB_detection %>%
    mutate(RF_bin = ntile(RF_score, 7)) %>%
    arrange(RF_score)

###############################################################################
# Log-transform quantitative metrics
###############################################################################

pDB_score_binned_log <- pDB_score_binned %>%
    mutate(
        Precursor.Quantity      = log10(Precursor.Quantity + 1),
        Precursor.Normalised    = log10(Precursor.Normalised + 1),
        Ms1.Total.Signal.Before = log10(Ms1.Total.Signal.Before + 1),
        Ms1.Area                = log10(Ms1.Area + 1),
        Ms1.Normalised          = log10(Ms1.Normalised + 1),
        Q.Value                 = -log10(Q.Value),
        PEP                     = -log10(PEP))

###############################################################################
# Define metrics for downstream analyses
###############################################################################

columns <- setdiff(names(pDB_score_binned_log), c("sequence", "RF_score", "RF_bin"))

###############################################################################
# Evaluate associations between RF scores and features
###############################################################################

lm_results <- lapply(columns, function(metric) {
        
        fit <- lm(reformulate(metric, response = "RF_score"),
            data = pDB_score_binned_log )
        
        fit_summary <- summary(fit)
        
        data.frame(
            metric = metric,
            r_squared = fit_summary$r.squared,
            adj_r_squared = fit_summary$adj.r.squared,
            slope = coef(fit)[2],
            p_value = coef(fit_summary)[2, 4])}) %>%
    bind_rows() %>%
    arrange(desc(r_squared))

cor_results <- lapply(columns, function(metric) {
        
        data.frame(metric = metric,
            rho = cor(
                pDB_score_binned_log$RF_score,
                pDB_score_binned_log[[metric]],
                method = "spearman",
                use = "complete.obs"))}) %>%
    bind_rows() %>%
    arrange(desc(abs(rho)))

correlation_results <- lm_results %>%
    left_join(cor_results, by = "metric") %>%
    arrange(desc(abs(rho)))

print(correlation_results)


###############################################################################
# Export processed datasets
###############################################################################

write_parquet(pDB_score_binned, "pDB_score_binned_raw.parquet")
write_parquet(pDB_score_binned_log, "pDB_score_binned_log.parquet")
write.table(corr, "pDB_correlation_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


