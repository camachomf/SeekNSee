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

correlation_results <- lapply(columns, function(metric) {
  
  corr <- cor.test(
    pDB_score_binned_log$RF_score,
    pDB_score_binned_log[[metric]],
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

write_parquet(pDB_score_binned, "pDB_score_binned_raw.parquet")
write_parquet(pDB_score_binned_log, "pDB_score_binned_log.parquet")
write.table(correlation_results, "pDB_correlation_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


