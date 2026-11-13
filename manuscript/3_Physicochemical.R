###############################################################################
# Analysis:
# Physicochemical evalutaion of peptides
#
# Purpose:
# Calculate peptide indexes using "Peptides" library
# and calculate metric associations with  score
#
# Input:
#  score_binned_features.parquet from 1_DIA or 2_DDA (R files)
#
# Output:
# DDA_physchem_corr.tsv
#
###############################################################################

###############################################################################
# Load libraries
###############################################################################
library(arrow)
library(dplyr)
library(Peptides)
library(tidyr)
library(ggplot2)

###############################################################################
# Import data
###############################################################################

physchem <- read_parquet("DDA_score_binned_features.parquet")

physchem <- physchem %>%
    rename(Stripped.Sequence = peptide) %>%
    select(Stripped.Sequence, peptide_score, score_bin)

###############################################################################
# Calculate physicochemical properties
###############################################################################

physchem <- physchem %>%
    mutate(
        GRAVYscore        = sapply(Stripped.Sequence, hydrophobicity),
        pI                = sapply(Stripped.Sequence, pI),
        Mw                = sapply(Stripped.Sequence, mw),
        Charge            = sapply(Stripped.Sequence, charge),
        BomanIndex        = sapply(Stripped.Sequence, boman),
        Instability       = sapply(Stripped.Sequence, instaIndex),
        AliphaticIndex    = sapply(Stripped.Sequence, aIndex),
        HydrophobicMoment = sapply(Stripped.Sequence, hmoment))

###############################################################################
# Correlation analysis
###############################################################################

physchem_metrics <- c(
    "Mw",
    "pI",
    "Charge",
    "BomanIndex",
    "Instability",
    "HydrophobicMoment",
    "GRAVYscore",
    "AliphaticIndex")


physchem_corr <- lapply(physchem_metrics, function(metric) {
  
  corr <- cor.test(
    physchem$peptide_score,
    physchem[[metric]],
    method = "spearman",
    exact = FALSE)
  
  data.frame(metric = metric,
    rho = unname(corr$estimate),
    p_value = corr$p.value)}) %>%
  bind_rows() %>%
  arrange(desc(abs(rho)))


###############################################################################
# Prepare heatmap data
###############################################################################

heatmap <- physchem %>%
    group_by(score_bin) %>%
    summarise(
        GRAVY             = median(GRAVYscore, na.rm = TRUE),
        AliphaticIndex    = median(AliphaticIndex, na.rm = TRUE),
        Mw                = median(Mw, na.rm = TRUE),
        pI                = median(pI, na.rm = TRUE),
        HydrophobicMoment = median(HydrophobicMoment, na.rm = TRUE),
        BomanIndex        = median(BomanIndex, na.rm = TRUE),
        Instability       = median(Instability, na.rm = TRUE),
        .groups = "drop")

heatmap <- heatmap %>%
    pivot_longer(-score_bin, names_to = "metric", values_to = "value")

heatmap <- heatmap %>%
    group_by(metric) %>%
    mutate(z = as.numeric(scale(value))) %>%
    ungroup()

heatmap$metric <- factor(heatmap$metric,
    levels = c(
        "Mw",
        "GRAVY",
        "AliphaticIndex",
        "pI",
        "HydrophobicMoment",
        "BomanIndex",
        "Instability"))

###############################################################################
# Generate heatmap
###############################################################################

p_heat <- ggplot(heatmap,
    aes(x = score_bin, y = metric, fill = z)) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_y_discrete(limits = rev(levels(heatmap$metric))) +
    scale_fill_gradient2(
        low = "#3d4b6bff",
        mid = "white",
        high = "#ccbb68ff",
        midpoint = 0) +
    labs(x = "Score Bin", y = NULL, fill = "Z-score") +
    theme_classic(base_size = 14)

print(p_heat)

###############################################################################
# Export 
###############################################################################

write.table(physchem_corr, "DDA_physchem_corr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
