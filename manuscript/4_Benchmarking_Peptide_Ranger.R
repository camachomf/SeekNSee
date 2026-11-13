###############################################################################
# Analysis:
# Benchmarking against PeptideRanger
#
# Purpose:
# Generate PeptideRanger predictions, evaluate peptide detection,
# compare classification performance against , and calculate
# ROC statistics.
#
# Inputs:
#   crossref.parquet
#   analytical_features.parquet
#   proteome FASTA
#
# Outputs:
#   CPTAC_detection.parquet
#   proteomicsDB_detection.parquet
#   
#
###############################################################################

###############################################################################
# Load libraries
###############################################################################

library(PeptideRanger)
library(arrow)
library(dplyr)
library(ggplot2)
library(pROC)


###############################################################################
# Generate theoretical peptidome using PeptideRanger
###############################################################################

peptidome <- create_peptidome(
    proteome_dir =
    "20260323_uniprotkb_UP000005640_AND_revi_20416entries_canonical.fasta",
    missed_cleavages = c(0, 2),
    synth_peps = FALSE,
    aa_range = c(7, 25))

###############################################################################
# Generate PeptideRanger predictions
###############################################################################

proteomicsDB_predictions <- peptide_predictions(
    peptides = peptidome$sequence,
    missed_cleavages = peptidome$missed_cleavages,
    prediction_model = RFmodel_ProteomicsDB)

CPTAC_predictions <- peptide_predictions(
    peptides = peptidome$sequence,
    missed_cleavages = peptidome$missed_cleavages,
    prediction_model = RFmodel_CPTAC)

###############################################################################
# Import processed data
###############################################################################

crossref <- read_parquet("DIA_crossref.parquet")

DIA <- read_parquet("DIA_analytical_features.parquet")

DDA <- read_parquet("DDA_analytical_features.parquet")

###############################################################################
# Integrate predictions with DIA data
# Alternatively, integrate with DDA data
###############################################################################

ranger_proteomicsDB_crossref <- left_join(
    proteomicsDB_predictions,
    DIA,
    by = c("sequence" = "Stripped.Sequence"))

ranger_CPTAC_crossref <- left_join(
    CPTAC_predictions,
    DIA,
    by = c("sequence" = "Stripped.Sequence"))


###############################################################################
# Define peptide detection status
# Use a different detection column for DDA/maxquant evidence.txt
###############################################################################

proteomicsDB_detection <- ranger_proteomicsDB_crossref %>%
    mutate(detected = ifelse(
            is.na(Precursor.Mz),
            "Not detected",
            "Detected"))

CPTAC_detection <- ranger_CPTAC_crossref %>%
    mutate(detected = ifelse(
            is.na(Precursor.Mz),
            "Not detected",
            "Detected"))

seeknsee_detection <- crossref %>%
    mutate(detected = ifelse(
            is.na(Peak.width),
            "Not detected",
            "Detected"))

###############################################################################
# Compute ROC curves
###############################################################################

roc_seeknsee <- roc(seeknsee_detection$detected, seeknsee_detection$peptide_score)

roc_proteomicsDB <- roc(proteomicsDB_detection$detected, proteomicsDB_detection$RF_score)

roc_CPTAC <- roc(CPTAC_detection$detected, CPTAC_detection$RF_score)

###############################################################################
# Prepare ROC plotting data
###############################################################################

roc_comparison_df <- bind_rows(
    data.frame(FPR = 1 - roc_seeknsee$specificities, TPR = roc_seeknsee$sensitivities,
        Model = paste0("SeekNSee (AUC = ", round(auc(roc_seeknsee), 2), ")")),
    
    data.frame(FPR = 1 - roc_proteomicsDB$specificities, TPR = roc_proteomicsDB$sensitivities,
        Model = paste0("ProteomicsDB (AUC = ", round(auc(roc_proteomicsDB), 2), ")")),
    
    data.frame(FPR = 1 - roc_CPTAC$specificities, TPR = roc_CPTAC$sensitivities,
        Model = paste0("CPTAC (AUC = ", round(auc(roc_CPTAC), 2), ")")))

###############################################################################
# Statistical comparison of ROC curves
###############################################################################

p_seeknsee_vs_proteomicsDB <- roc.test(
    roc_seeknsee,
    roc_proteomicsDB,
    method = "delong")$p.value

p_seeknsee_vs_CPTAC <- roc.test(
  roc_seeknsee,
  roc_CPTAC,
  method = "delong")$p.value


p_seeknsee_vs_CPTAC
p_seeknsee_vs_proteomicsDB 

###############################################################################
# Plot ROC comparison
###############################################################################

p_roc_comparison <- ggplot(roc_comparison_df,
    aes(x = FPR, y = TPR, colour = Model)) +
    geom_line(linewidth = 1.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    scale_colour_manual(
        values = c("#2d4789ff", "#2a6892ff", "#3bb29dff")) +
    annotate("text", x = 0.65, y = 0.50, hjust = 0, size = 3.5, 
             label = paste0("DeLong p = ", format.pval(p_seeknsee_vs_proteomicsDB))) +
    theme_bw() +
    theme(legend.position = c(0.75, 0.25)) +
    labs(x = "False Positive Rate", y = "True Positive Rate", colour = NULL) 

print(p_roc_comparison)


###############################################################################
# Results table
###############################################################################
ci_seeknsee <- ci.auc(roc_seeknsee)
ci_proteomicsDB <- ci.auc(roc_proteomicsDB)
ci_CPTAC <- ci.auc(roc_CPTAC)


auc_table <- data.frame(
    Model = c("SeekNSee", "ProteomicsDB", "CPTAC"),
    AUC = c(
        auc(roc_seeknsee),
        auc(roc_proteomicsDB),
        auc(roc_CPTAC)),
    Lower = c(
        ci_seeknsee[1],
        ci_proteomicsDB[1],
        ci_CPTAC[1]),
    Upper = c(
        ci_seeknsee[3],
        ci_proteomicsDB[3],
        ci_CPTAC[3])) %>%
    transmute(Model, `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)",
            AUC, Lower, Upper))

###############################################################################
# Export Peptide Ranger 
###############################################################################
write_parquet(CPTAC_detection, "DIA_CPTAC_detection.parquet")
write_parquet(proteomicsDB_detection, "DIA_proteomicsDB_detection.parquet")
write.table(auc_table, "auc_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


