# ============================================================
# Parkinson's Disease RNA-seq Meta-Analysis
# Summary of Random-Effects Meta-Analysis
# ============================================================

library(tidyverse)
library(rio)

# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

base_dir <- "/media/nabil/ssd/Parkinson-s"

meta_dir <- file.path(
  base_dir,
  "results/tables/meta-analysis"
)

venn_dir <- file.path(
  base_dir,
  "results/tables/venn"
)

dir.create(
  meta_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Load random-effects meta-analysis result
# ------------------------------------------------------------

random_file <- file.path(
  meta_dir,
  "Parkinson_random_effects_meta_analysis.csv"
)

meta_result <- import(random_file)


cat("\n============================================================\n")
cat("PARKINSON'S META-ANALYSIS SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Total genes analyzed: ",
  nrow(meta_result),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Standardize columns
# ------------------------------------------------------------

meta_result <- meta_result |>
  mutate(
    Gene_ID = as.character(Gene_ID),
    Gene_Symbol = as.character(Gene_Symbol),
    randomSummary = as.numeric(randomSummary),
    randomP = as.numeric(randomP)
  )


# ------------------------------------------------------------
# 4. Define significance and regulation
#
# Criteria:
#
# randomP < 0.05
# abs(randomSummary) >= 1
#
# randomSummary > 0 = higher in PD
# randomSummary < 0 = lower in PD
# ------------------------------------------------------------

meta_result <- meta_result |>
  mutate(
    
    Regulation = case_when(
      randomSummary > 0 ~ "UP",
      randomSummary < 0 ~ "DOWN",
      TRUE ~ "No change"
    ),
    
    Significance = case_when(
      
      !is.na(randomP) &
        !is.na(randomSummary) &
        randomP < 0.05 &
        randomSummary >= 1 ~
        "Significant UP",
      
      !is.na(randomP) &
        !is.na(randomSummary) &
        randomP < 0.05 &
        randomSummary <= -1 ~
        "Significant DOWN",
      
      TRUE ~ "NS"
    )
  )


# ------------------------------------------------------------
# 5. Calculate ordinary fold-change
#
# log2FC -> FC
#
# FC = 2^log2FC
# ------------------------------------------------------------

meta_result <- meta_result |>
  mutate(
    Fold_Change = 2^randomSummary
  )


# ------------------------------------------------------------
# 6. Significant meta-DEGs
# ------------------------------------------------------------

significant_meta <- meta_result |>
  filter(
    Significance != "NS"
  ) |>
  arrange(
    randomP
  )


up_genes <- significant_meta |>
  filter(
    Significance == "Significant UP"
  )


down_genes <- significant_meta |>
  filter(
    Significance == "Significant DOWN"
  )


# ------------------------------------------------------------
# 7. Print counts
# ------------------------------------------------------------

cat(
  "Significant meta-DEGs: ",
  nrow(significant_meta),
  "\n",
  sep = ""
)

cat(
  "Significant UP: ",
  nrow(up_genes),
  "\n",
  sep = ""
)

cat(
  "Significant DOWN: ",
  nrow(down_genes),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Summary statistics table
# ------------------------------------------------------------

total_genes <- nrow(meta_result)

gene_stats <- tibble(
  
  Category = c(
    "Total genes analyzed",
    "Significant meta-DEGs",
    "Significant UP",
    "Significant DOWN"
  ),
  
  Count = c(
    total_genes,
    nrow(significant_meta),
    nrow(up_genes),
    nrow(down_genes)
  ),
  
  Percentage_of_all_genes = c(
    
    100,
    
    nrow(significant_meta) /
      total_genes * 100,
    
    nrow(up_genes) /
      total_genes * 100,
    
    nrow(down_genes) /
      total_genes * 100
    
  )
)


print(gene_stats)


export(
  gene_stats,
  file.path(
    meta_dir,
    "Parkinson_meta_analysis_summary_stats.csv"
  )
)


# ------------------------------------------------------------
# 9. Significant annotated genes
# ------------------------------------------------------------

annotated_genes <- significant_meta |>
  
  filter(
    !is.na(Gene_Symbol),
    Gene_Symbol != ""
  ) |>
  
  select(
    Gene_ID,
    Gene_Symbol,
    Gene_Biotype,
    Chromosome,
    Start,
    End,
    Strand,
    randomSummary,
    randomCi.lb,
    randomCi.ub,
    randomP,
    Fold_Change,
    Regulation,
    Significance,
    Meta_Direction,
    Number_Up,
    Number_Down,
    Number_Datasets,
    Direction_Consistency,
    het_QE,
    het_QEp
  ) |>
  
  arrange(
    randomP
  )


export(
  annotated_genes,
  file.path(
    meta_dir,
    "Parkinson_significant_meta_DEGs_annotated.csv"
  )
)


# ------------------------------------------------------------
# 10. UP-regulated genes
# ------------------------------------------------------------

up_table <- annotated_genes |>
  
  filter(
    Regulation == "UP"
  ) |>
  
  arrange(
    randomP
  )


export(
  up_table,
  file.path(
    meta_dir,
    "Parkinson_meta_DEGs_UP.csv"
  )
)


# ------------------------------------------------------------
# 11. DOWN-regulated genes
# ------------------------------------------------------------

down_table <- annotated_genes |>
  
  filter(
    Regulation == "DOWN"
  ) |>
  
  arrange(
    randomP
  )


export(
  down_table,
  file.path(
    meta_dir,
    "Parkinson_meta_DEGs_DOWN.csv"
  )
)


# ------------------------------------------------------------
# 12. Direction consistency summary
# ------------------------------------------------------------

direction_summary <- significant_meta |>
  
  count(
    Direction_Consistency,
    name = "Count"
  ) |>
  
  mutate(
    Percentage =
      Count /
      sum(Count) *
      100
  )


print(direction_summary)


export(
  direction_summary,
  file.path(
    meta_dir,
    "Parkinson_direction_consistency_summary.csv"
  )
)


# ------------------------------------------------------------
# 13. Consistently directed significant genes
# ------------------------------------------------------------

consistent_genes <- significant_meta |>
  
  filter(
    Direction_Consistency %in%
      c(
        "Consistently Up",
        "Consistently Down"
      )
  ) |>
  
  select(
    Gene_ID,
    Gene_Symbol,
    Gene_Biotype,
    randomSummary,
    randomCi.lb,
    randomCi.ub,
    randomP,
    Fold_Change,
    Regulation,
    Meta_Direction,
    Number_Up,
    Number_Down,
    Number_Datasets,
    Direction_Consistency
  ) |>
  
  arrange(
    randomP
  )


export(
  consistent_genes,
  file.path(
    meta_dir,
    "Parkinson_consistent_meta_DEGs.csv"
  )
)


# ------------------------------------------------------------
# 14. Discordant significant genes
# ------------------------------------------------------------

discordant_genes <- significant_meta |>
  
  filter(
    Direction_Consistency == "Discordant"
  ) |>
  
  select(
    Gene_ID,
    Gene_Symbol,
    Gene_Biotype,
    randomSummary,
    randomCi.lb,
    randomCi.ub,
    randomP,
    Fold_Change,
    Regulation,
    Meta_Direction,
    Number_Up,
    Number_Down,
    Number_Datasets,
    Direction_Consistency
  ) |>
  
  arrange(
    randomP
  )


export(
  discordant_genes,
  file.path(
    meta_dir,
    "Parkinson_discordant_meta_DEGs.csv"
  )
)


# ------------------------------------------------------------
# 15. Final summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Total genes analyzed: ",
  total_genes,
  "\n",
  sep = ""
)

cat(
  "Significant meta-DEGs: ",
  nrow(significant_meta),
  "\n",
  sep = ""
)

cat(
  "UP-regulated: ",
  nrow(up_genes),
  "\n",
  sep = ""
)

cat(
  "DOWN-regulated: ",
  nrow(down_genes),
  "\n",
  sep = ""
)

cat(
  "Consistently directed: ",
  nrow(consistent_genes),
  "\n",
  sep = ""
)

cat(
  "Discordant: ",
  nrow(discordant_genes),
  "\n",
  sep = ""
)

cat("\nFiles created in:\n")
cat(meta_dir)
cat("\n\n")

cat("1. Parkinson_meta_analysis_summary_stats.csv\n")
cat("2. Parkinson_significant_meta_DEGs_annotated.csv\n")
cat("3. Parkinson_meta_DEGs_UP.csv\n")
cat("4. Parkinson_meta_DEGs_DOWN.csv\n")
cat("5. Parkinson_direction_consistency_summary.csv\n")
cat("6. Parkinson_consistent_meta_DEGs.csv\n")
cat("7. Parkinson_discordant_meta_DEGs.csv\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")