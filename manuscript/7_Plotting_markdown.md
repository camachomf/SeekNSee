Correlation Visualizations
================

<div style="font-size: 85%;">

### Analysis

**Correlation analysis and visualization of score-associated features**

**Purpose:**  
Quantify the relationship between scores and analytical metrics,
generate correlation summaries, and visualize score-bin trends.

**Description:**  
This workflow evaluates associations between peptide scores and
identification, chromatographic, abundance, and quantitative quality
metrics. Correlations are calculated using the Spearman method and
summarized through graphical representations across score categories.

**Example dataset:**  
Deep-fractionated proteome data are provided as an example.

**Inputs** - score_binned_features.parquet -
score_binned_features_raw.parquet - crossref.parquet -
correlation_results.tsv

**Outputs** - Correlation figures - Score-category visualizations -
Summary correlation tables

</div>

``` r
library(arrow)
```

    ## Some features are not enabled in this build of Arrow. Run `arrow_info()` for more information.

    ## 
    ## Attaching package: 'arrow'

    ## The following object is masked from 'package:utils':
    ## 
    ##     timestamp

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(ggridges)
library(svglite)
library(forcats)
library(scales)
```

``` r
score_binned_features <- read_parquet("score_binned_features.parquet")

score_binned_features_raw <- read_parquet("score_binned_features_raw.parquet")

crossref <- read_parquet("crossref.parquet")

correlation_results <- read.table("correlation_results.tsv", sep = "\t", header = TRUE)
```

``` r
columns <- setdiff(names(score_binned_features), c("peptide_score", "score_bin", "peptide", "Precursor.Charge", "Ms1.Total.Signal.Before"))

quality_metrics <- c(
    "Length",
    "FWHM",
    "Peak.width",
    "Empirical.Quality",
    "Quantity.Quality",
    "Evidence",
    "Q.Value",
    "PEP")
```

``` r
p_corr_bar <- ggplot(
    correlation_results,
    aes(x = reorder(metric, rho), y = rho)) +
    geom_col(fill = "#413980ff") +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 11)) +
    labs(x = NULL,  y = expression(rho))

print(p_corr_bar)
```

![](7_Plotting_markdown_files/figure-gfm/elation-summary-1.png)<!-- -->

``` r
p_corr_dot <- ggplot(correlation_results,
    aes(x = rho, y = reorder(metric, rho))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 4, colour = "#413980ff") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 11)) +
    labs(x = expression(rho),  y = NULL)

print(p_corr_dot)
```

![](7_Plotting_markdown_files/figure-gfm/elation-summary-2.png)<!-- -->

``` r
for (metric in columns) {
        rho <- correlation_results %>%
        filter(metric == !!metric) %>%
        pull(rho)
    
    pval <- correlation_results %>%
        filter(metric == !!metric) %>%
        pull(p_value)
    
    summary_df <- score_binned_features %>%
        group_by(score_bin) %>%
        summarise(median = median(.data[[metric]], na.rm = TRUE),
            q25 = quantile(.data[[metric]], 0.25, na.rm = TRUE),
            q75 = quantile(.data[[metric]], 0.75, na.rm = TRUE  ),
            .groups = "drop")
    
    p_metric <- ggplot(summary_df,
        aes(x = score_bin, y = median, group = 1)) +
        geom_ribbon(aes(ymin = q25, ymax = q75),
            fill = "cornflowerblue",
            alpha = 0.25) +
        geom_line(linewidth = 1.5, colour = "royalblue3") +
        geom_point(size = 4, colour = "royalblue3") +
        annotate("text", x = Inf, y = Inf,
                label = paste0("ρ = ", round(rho, 3), 
                "\n", "p = ",    
                format.pval(pval, digits = 3, eps = .Machine$double.xmin)),
            hjust = 1.1, vjust = 3, size = 5) +
        theme_bw() +
        theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 11)) +
        labs(x = "Score bin", y = metric)
    
    print(p_metric)
    }
```

![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-1.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-2.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-3.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-4.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-5.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-6.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-7.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-8.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-9.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-10.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-11.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-12.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-13.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/score-bin-trends-14.png)<!-- -->

``` r
for (metric in columns) {
    
    p_distribution <- ggplot(score_binned_features,
        aes(x = .data[[metric]], y = score_bin, fill = score_bin)) +
        geom_density_ridges(alpha = 0.7, scale = 1) +
        scale_fill_viridis_d(option = "mako", direction = 1, begin = 0.1, end = 0.8) +
        labs(x = metric, y = "Score bin")
        print(p_distribution)
}
```

    ## Picking joint bandwidth of 0.587

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-1.png)<!-- -->

    ## Picking joint bandwidth of 19.1

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-2.png)<!-- -->

    ## Picking joint bandwidth of 0.18

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-3.png)<!-- -->

    ## Picking joint bandwidth of 0.147

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-4.png)<!-- -->

    ## Picking joint bandwidth of 0.187

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-5.png)<!-- -->

    ## Picking joint bandwidth of 0.126

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-6.png)<!-- -->

    ## Picking joint bandwidth of 0.172

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-7.png)<!-- -->

    ## Picking joint bandwidth of 0.00109

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-8.png)<!-- -->

    ## Picking joint bandwidth of 0.00251

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-9.png)<!-- -->

    ## Picking joint bandwidth of 0.00813

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-10.png)<!-- -->

    ## Picking joint bandwidth of 0.00581

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-11.png)<!-- -->

    ## Picking joint bandwidth of 0.139

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-12.png)<!-- -->

    ## Picking joint bandwidth of 0.0446

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-13.png)<!-- -->

    ## Picking joint bandwidth of 0.115

![](7_Plotting_markdown_files/figure-gfm/ridge-distributions-14.png)<!-- -->

``` r
for (metric in columns) {

boxplot <- ggplot(score_binned_features,
  aes(x = factor(score_bin), y = .data[[metric]], fill = factor(score_bin))) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "mako", direction = -1, begin = 0.3, end = 1) +
  theme_classic() +
  theme(legend.position = "none") + 
    theme(axis.title = element_text(size = 15),
               axis.text = element_text(size = 11)) +
  labs(x = "Score bin", y = metric)
  print(boxplot)
}
```

![](7_Plotting_markdown_files/figure-gfm/boxplots-1.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-2.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-3.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-4.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-5.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-6.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-7.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-8.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-9.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-10.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-11.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-12.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-13.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/boxplots-14.png)<!-- -->

``` r
for (metric in columns) {
    
    summary_df <- score_binned_features_raw %>%
        group_by(score_bin) %>%
        summarise(mean = mean(.data[[metric]]),
            sd = sd(.data[[metric]]),
            n = n(),
            se = sd / sqrt(n),
            .groups = "drop")
    
    barplot_se <- ggplot(summary_df,
        aes(x = score_bin, y = mean, fill = score_bin)) +
        geom_col(width = 0.8) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
        scale_fill_viridis_d(option = "mako", direction = -1, begin = 0.3, end = 1) +
        theme_classic(base_size = 14) +
        labs(title = paste(metric, ""), x = "Score bin", y = metric) +
        theme(legend.position = "none")
    print(barplot_se)
}
```

![](7_Plotting_markdown_files/figure-gfm/barplots-1.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-2.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-3.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-4.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-5.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-6.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-7.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-8.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-9.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-10.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-11.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-12.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-13.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/barplots-14.png)<!-- -->

``` r
for (metric in quality_metrics) {
    
    stacked <- score_binned_features %>%
        mutate(value_bin = cut_number(.data[[metric]], n = 4), 
               value_bin = fct_rev(value_bin)) %>%
        ggplot(aes(x = score_bin, fill = value_bin)) +
        geom_bar(position = "fill") +
        scale_y_continuous(labels = percent) +
        scale_fill_viridis_d(option = "mako", direction = 1, begin = 0.1, end = 0.8) +
        labs(title = metric, x = "Score Bin", y = "Proportion", fill = metric) +
        theme_classic()
    
    print(stacked)
}
```

![](7_Plotting_markdown_files/figure-gfm/stacked-quality-1.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-2.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-3.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-4.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-5.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-6.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-7.png)<!-- -->![](7_Plotting_markdown_files/figure-gfm/stacked-quality-8.png)<!-- -->

``` r
detection <- crossref %>%
    mutate(detected = ifelse(is.na(Peak.width), "Not detected", "Detected"))

p_detection <- ggplot(
    detection,
    aes(x = factor(peptide_score), fill = detected)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("Detected" = "#413980ff", "Not detected" = "#d7f4dfff")) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
    labs( x = "Score", y = "Stripped Sequences (%)", fill = NULL)

print(p_detection)
```

![](7_Plotting_markdown_files/figure-gfm/detection-frequency-1.png)<!-- -->
