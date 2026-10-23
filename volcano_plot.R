# ============================================================
# Parkinson's Disease RNA-seq Meta-Analysis
# Random-Effects Meta-Analysis Volcano Plot
#
# Datasets:
#   GSE133101
#   GSE136666
#   GSE168496
#
# Input:
#   Parkinson_random_effects_meta_analysis.csv
#
# Effect size:
#   randomSummary
#
# P-value:
#   randomP
#
# Significance:
#   P < 0.05
#   |randomSummary| >= 1
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(ggrepel)
library(rio)


# ------------------------------------------------------------
# 2. Pin dplyr verbs
# ------------------------------------------------------------

mutate    <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize
arrange   <- dplyr::arrange
rename    <- dplyr::rename
count     <- dplyr::count
desc      <- dplyr::desc
select    <- dplyr::select
filter    <- dplyr::filter


# ------------------------------------------------------------
# 3. Project directories
# ------------------------------------------------------------

base_dir <- "/media/nabil/ssd/Parkinson-s"

meta_dir <- file.path(
  base_dir,
  "results/tables/meta-analysis"
)

figure_dir <- file.path(
  base_dir,
  "results/figures/meta-analysis"
)


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load random-effects meta-analysis results
# ------------------------------------------------------------

meta_file <- file.path(
  meta_dir,
  "Parkinson_random_effects_meta_analysis.csv"
)


meta_result <- import(meta_file)


# ------------------------------------------------------------
# 5. Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "Gene_ID",
  "Gene_Symbol",
  "randomSummary",
  "randomP"
)


missing_columns <- setdiff(
  required_columns,
  colnames(meta_result)
)


if (length(missing_columns) > 0) {
  
  stop(
    "The following required columns are missing:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 6. Prepare volcano plot data
# ------------------------------------------------------------

meta_result <- meta_result |>
  
  select(
    Gene_ID,
    Gene_Symbol,
    randomSummary,
    randomP
  ) |>
  
  mutate(
    
    Gene_ID =
      as.character(Gene_ID),
    
    Gene_Symbol =
      as.character(Gene_Symbol),
    
    log2FC =
      as.numeric(randomSummary),
    
    P.Value =
      as.numeric(randomP)
    
  )


# ------------------------------------------------------------
# 7. Remove unusable values
# ------------------------------------------------------------

meta_result <- meta_result |>
  
  filter(
    !is.na(log2FC),
    !is.na(P.Value),
    P.Value > 0,
    P.Value <= 1
  )


# ------------------------------------------------------------
# 8. Define significance
#
# P < 0.05
# |log2FC| >= 1
#
# Positive randomSummary:
#   higher expression in Parkinson's disease
#
# Negative randomSummary:
#   lower expression in Parkinson's disease
# ------------------------------------------------------------

meta_result <- meta_result |>
  
  mutate(
    
    Significance = case_when(
      
      P.Value < 0.05 &
        log2FC >= 1 ~
        "Up",
      
      P.Value < 0.05 &
        log2FC <= -1 ~
        "Down",
      
      TRUE ~
        "Not significant"
      
    )
    
  )


# ------------------------------------------------------------
# 9. Print numbers of genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PARKINSON'S META-ANALYSIS VOLCANO PLOT\n")
cat("============================================================\n\n")


cat(
  "Total genes plotted: ",
  nrow(meta_result),
  "\n",
  sep = ""
)


cat(
  "Significant UP: ",
  sum(
    meta_result$Significance == "Up"
  ),
  "\n",
  sep = ""
)


cat(
  "Significant DOWN: ",
  sum(
    meta_result$Significance == "Down"
  ),
  "\n",
  sep = ""
)


cat(
  "Not significant: ",
  sum(
    meta_result$Significance == "Not significant"
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Select genes for labeling
#
# Only annotated, significant genes.
#
# We label the genes with the largest absolute
# meta-analysis effect sizes.
# ------------------------------------------------------------

top_genes <- meta_result |>
  
  mutate(
    Gene_Symbol =
      na_if(
        Gene_Symbol,
        ""
      )
  ) |>
  
  filter(
    !is.na(Gene_Symbol),
    P.Value < 0.05,
    abs(log2FC) >= 1
  ) |>
  
  slice_max(
    order_by = abs(log2FC),
    n = 30,
    with_ties = FALSE
  )


# ------------------------------------------------------------
# 11. Calculate x-axis limits
# ------------------------------------------------------------

x_lim <- ceiling(
  max(
    abs(meta_result$log2FC),
    na.rm = TRUE
  )
)


# Prevent an excessively narrow axis

if (x_lim < 2) {
  x_lim <- 2
}


# Choose sensible x-axis spacing

x_break <- ifelse(
  x_lim > 10,
  2,
  1
)


# ------------------------------------------------------------
# 12. Volcano plot
# ------------------------------------------------------------

volcano <- ggplot(
  
  meta_result,
  
  aes(
    x = log2FC,
    y = -log10(P.Value),
    color = Significance
  )
  
) +
  
  
  # Points
  
  geom_point(
    alpha = 0.8,
    size = 2.5
  ) +
  
  
  # Log2FC threshold lines
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  
  
  # P-value threshold
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  
  
  # Colors
  
  scale_color_manual(
    
    values = c(
      "Down" = "#2c7fb8",
      "Not significant" = "#636363",
      "Up" = "#e34a33"
    )
    
  ) +
  
  
  # Labels
  
  labs(
    
    title = "Parkinson's Disease Meta-Analysis",
    
    subtitle =
      "Random-effects model",
    
    x = "Meta-analysis log2 fold change",
    
    y = "-log10(P-value)",
    
    color = "Significance"
    
  ) +
  
  
  # Theme
  
  theme_minimal() +
  
  
  theme(
    
    plot.title =
      element_text(
        size = 18,
        face = "bold"
      ),
    
    plot.subtitle =
      element_text(
        size = 14
      ),
    
    legend.title =
      element_text(
        size = 15,
        face = "bold"
      ),
    
    legend.text =
      element_text(
        size = 14
      ),
    
    legend.position =
      "right",
    
    axis.title.x =
      element_text(
        size = 14
      ),
    
    axis.title.y =
      element_text(
        size = 14
      ),
    
    axis.text.x =
      element_text(
        colour = "black",
        size = 12
      ),
    
    axis.text.y =
      element_text(
        colour = "black",
        size = 12
      ),
    
    panel.border =
      element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    
  ) +
  
  
  # Gene labels
  
  geom_text_repel(
    
    data = top_genes,
    
    aes(
      label = Gene_Symbol
    ),
    
    colour = "black",
    
    size = 3.5,
    
    max.overlaps = 30,
    
    direction = "both",
    
    max.time = 5,
    
    force = 5,
    
    force_pull = 5,
    
    point.padding = 0.5,
    
    box.padding = 0.5,
    
    seed = 40
    
  ) +
  
  
  # Symmetric x-axis
  
  scale_x_continuous(
    
    limits =
      c(
        -x_lim,
        x_lim
      ),
    
    breaks =
      seq(
        -x_lim,
        x_lim,
        by = x_break
      )
    
  )


# ------------------------------------------------------------
# 13. Display plot
# ------------------------------------------------------------

print(volcano)


# ------------------------------------------------------------
# 14. Save plot
# ------------------------------------------------------------

volcano_file <- file.path(
  
  figure_dir,
  
  "Parkinson_random_effects_volcano_plot.png"
  
)


ggsave(
  
  volcano_file,
  
  plot = volcano,
  
  width = 12,
  
  height = 12,
  
  dpi = 600,
  
  bg = "white"
  
)


cat(
  "\nVolcano plot saved to:\n",
  volcano_file,
  "\n"
)


cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")