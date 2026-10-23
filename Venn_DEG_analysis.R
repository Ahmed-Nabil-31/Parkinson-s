# ============================================================
# COMMON DIFFERENTIALLY EXPRESSED GENES
# Parkinson's disease substantia nigra RNA-seq
#
# Datasets:
#   GSE133101
#   GSE136666
#   GSE168496
#
# Input:
#   results/tables/DESeq2/*.csv
#
# DEG criteria:
#   padj < 0.05
#   |log2FoldChange| > 1
#
# Outputs:
#   Venn diagram
#   Complete intersection summary
#   Genes shared by all datasets
#   DEG recurrence/frequency table
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(rio)


# ============================================================
# 2. PIN DPLYR VERBS
# ============================================================

mutate <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize
arrange <- dplyr::arrange
rename <- dplyr::rename
count <- dplyr::count
desc <- dplyr::desc
select <- dplyr::select
filter <- dplyr::filter


# ============================================================
# 3. DIRECTORIES
# ============================================================

base_dir <- "/media/nabil/ssd/Parkinson-s"

input_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "DESeq2"
)

figure_dir <- file.path(
  base_dir,
  "results",
  "figures",
  "venn"
)

table_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "venn"
)


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. DATASETS
# ============================================================

geo_ids <- c(
  "GSE133101",
  "GSE136666",
  "GSE168496"
)


# ============================================================
# 5. DEG THRESHOLDS
# ============================================================

sig_cutoff <- 0.05

lfc_cutoff <- 1


# ============================================================
# 6. FIND DESeq2 FILES
# ============================================================

deg_files <- file.path(
  input_dir,
  paste0(
    geo_ids,
    ".csv"
  )
)


# ============================================================
# 7. CHECK FILES EXIST
# ============================================================

missing_files <- deg_files[
  !file.exists(deg_files)
]


if (length(missing_files) > 0) {
  
  stop(
    "\nThe following DESeq2 result files were not found:\n\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\n"
  )
}


# ============================================================
# 8. NAME FILES BY DATASET
# ============================================================

names(deg_files) <- geo_ids


# ============================================================
# 9. LOAD AND FILTER DEGs
# ============================================================

gene_sets <- list()


for (geo_id in geo_ids) {
  
  cat("\n")
  cat("============================================================\n")
  cat("PROCESSING ", geo_id, "\n", sep = "")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # Load DESeq2 results
  # ----------------------------------------------------------
  
  deg <- readr::read_csv(
    deg_files[[geo_id]],
    show_col_types = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Check required columns
  # ----------------------------------------------------------
  
  required_columns <- c(
    "Gene_ID",
    "padj",
    "log2FoldChange"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(deg)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      "\nMissing required columns in ",
      geo_id,
      ":\n",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\nActual columns:\n",
      paste(
        colnames(deg),
        collapse = ", "
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Identify significant DEGs
  # ----------------------------------------------------------
  
  sig_deg <- deg |>
    
    filter(
      
      !is.na(padj),
      
      padj < sig_cutoff,
      
      !is.na(log2FoldChange),
      
      abs(log2FoldChange) >
        lfc_cutoff
      
    ) |>
    
    mutate(
      
      Gene_ID =
        as.character(Gene_ID)
      
    ) |>
    
    filter(
      
      !is.na(Gene_ID),
      
      Gene_ID != ""
      
    ) |>
    
    distinct(
      Gene_ID,
      .keep_all = TRUE
    )
  
  
  # ----------------------------------------------------------
  # Store gene set
  # ----------------------------------------------------------
  
  gene_sets[[geo_id]] <-
    unique(
      sig_deg$Gene_ID
    )
  
  
  # ----------------------------------------------------------
  # Print summary
  # ----------------------------------------------------------
  
  cat(
    "Significant DEGs: ",
    nrow(sig_deg),
    "\n",
    sep = ""
  )
  
}


# ============================================================
# 10. SAVE DEG SET SIZE SUMMARY
# ============================================================

deg_set_summary <- tibble(
  
  Dataset =
    names(gene_sets),
  
  Significant_DEGs =
    map_int(
      gene_sets,
      length
    )
  
)


write_csv(
  
  deg_set_summary,
  
  file.path(
    table_dir,
    "DEG_set_summary.csv"
  )
  
)


# ============================================================
# 11. CHECK EMPTY SETS
# ============================================================

empty_sets <- names(gene_sets)[
  map_int(
    gene_sets,
    length
  ) == 0
]


if (length(empty_sets) > 0) {
  
  warning(
    "\nNo significant DEGs found in:\n",
    paste(
      empty_sets,
      collapse = ", "
    )
  )
}


# ============================================================
# 12. REMOVE EMPTY SETS FOR VENN ANALYSIS
# ============================================================

gene_sets_plot <- gene_sets[
  map_int(
    gene_sets,
    length
  ) > 0
]


n_sets <- length(
  gene_sets_plot
)


# ============================================================
# 13. VENN DIAGRAM
# ============================================================

if (
  n_sets >= 2 &&
  n_sets <= 5
) {
  
  fill_cols <- c(
    "#99d8c9",
    "#addd8e",
    "#bcbddc",
    "#fec44f",
    "#fc9272"
  )[seq_len(n_sets)]
  
  
  cat_cols <- c(
    "#2ca25f",
    "#31a354",
    "#756bb1",
    "#e6550d",
    "#de2d26"
  )[seq_len(n_sets)]
  
  
  venn_file <- file.path(
    figure_dir,
    "PD_DEG_Venn_diagram.png"
  )
  
  
  png(
    
    venn_file,
    
    units = "in",
    
    width = 7,
    
    height = 7,
    
    res = 600
    
  )
  
  
  venn_plot <- venn.diagram(
    
    x =
      gene_sets_plot,
    
    filename =
      NULL,
    
    fill =
      fill_cols,
    
    alpha =
      0.6,
    
    cex =
      0.8,
    
    cat.cex =
      1,
    
    cat.col =
      cat_cols,
    
    margin =
      0.12,
    
    lwd =
      2,
    
    disable.logging =
      TRUE
    
  )
  
  
  grid::grid.draw(
    venn_plot
  )
  
  
  dev.off()
  
  
  cat(
    "\nVenn diagram saved:\n",
    venn_file,
    "\n"
  )
  
  
} else if (
  n_sets > 5
) {
  
  upset_plot <-
    ggVennDiagram::ggVennDiagram(
      gene_sets_plot,
      force_upset = TRUE
    )
  
  
  upset_file <- file.path(
    figure_dir,
    "PD_DEG_UpSet_plot.png"
  )
  
  
  ggsave(
    
    upset_file,
    
    upset_plot,
    
    width =
      12,
    
    height =
      7,
    
    units =
      "in",
    
    bg =
      "white",
    
    dpi =
      600
    
  )
  
  
} else {
  
  warning(
    "\nFewer than two non-empty DEG sets are available.",
    "\nVenn/UpSet plot was not generated."
  )
  
}


# ============================================================
# 14. CALCULATE INTERSECTION REGIONS
# ============================================================

if (n_sets >= 1) {
  
  venn_object <-
    ggVennDiagram::Venn(
      gene_sets
    )
  
  
  regions <-
    ggVennDiagram::process_region_data(
      venn_object
    )
  
  
  # ----------------------------------------------------------
  # Full intersection summary
  # ----------------------------------------------------------
  
  regions_out <- regions |>
    
    select(
      id,
      name,
      count
    ) |>
    
    arrange(
      desc(count)
    )
  
  
  write_csv(
    
    as.data.frame(
      regions_out
    ),
    
    file.path(
      table_dir,
      "intersection_summary.csv"
    )
    
  )
  
}


# ============================================================
# 15. GENES SHARED BY ALL THREE DATASETS
# ============================================================

common_genes <- Reduce(
  intersect,
  gene_sets_plot
)


common_df <- tibble(
  
  Gene_ID =
    as.character(
      common_genes
    )
  
)


write_csv(
  
  common_df,
  
  file.path(
    table_dir,
    "common_genes_all_datasets.csv"
  )
  
)


cat(
  "\nGenes shared by ALL datasets: ",
  length(common_genes),
  "\n",
  sep = ""
)


# ============================================================
# 16. DEG FREQUENCY TABLE
# ============================================================
#
# Shows how many datasets each gene is significantly
# differentially expressed in.
#
# This is especially useful if the three-way intersection
# is small or empty.
# ============================================================

deg_frequency <- gene_sets_plot |>
  
  unlist() |>
  
  table() |>
  
  as.data.frame() |>
  
  setNames(
    c(
      "Gene_ID",
      "n_datasets"
    )
  ) |>
  
  arrange(
    desc(n_datasets)
  )


write_csv(
  
  deg_frequency,
  
  file.path(
    table_dir,
    "deg_frequency.csv"
  )
  
)


# ============================================================
# 17. RECURRENT DEGs
# ============================================================

recurrent_genes <- deg_frequency |>
  
  filter(
    n_datasets >= 2
  )


write_csv(
  
  recurrent_genes,
  
  file.path(
    table_dir,
    "recurrent_DEGs_2_or_more_datasets.csv"
  )
  
)


cat(
  "Genes recurrent in >=2 datasets: ",
  nrow(recurrent_genes),
  "\n",
  sep = ""
)


# ============================================================
# 18. SAVE SCRIPT PARAMETERS
# ============================================================

analysis_parameters <- tibble(
  
  Parameter = c(
    "padj_cutoff",
    "absolute_log2FoldChange_cutoff",
    "datasets",
    "comparison"
  ),
  
  Value = c(
    as.character(sig_cutoff),
    as.character(lfc_cutoff),
    paste(
      geo_ids,
      collapse = ", "
    ),
    "PD vs Control"
  )
  
)


write_csv(
  
  analysis_parameters,
  
  file.path(
    table_dir,
    "Venn_analysis_parameters.csv"
  )
  
)


# ============================================================
# 19. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Venn / DEG OVERLAP ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "Datasets analyzed: ",
  paste(
    geo_ids,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Significance cutoff: padj < ",
  sig_cutoff,
  "\n",
  sep = ""
)

cat(
  "Fold-change cutoff: |log2FC| > ",
  lfc_cutoff,
  "\n",
  sep = ""
)

cat(
  "Common genes across all datasets: ",
  length(common_genes),
  "\n",
  sep = ""
)

cat(
  "Recurrent genes in >=2 datasets: ",
  nrow(recurrent_genes),
  "\n",
  sep = ""
)

cat("\nOutputs:\n")

cat(
  file.path(
    table_dir,
    "DEG_set_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    table_dir,
    "intersection_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    table_dir,
    "common_genes_all_datasets.csv"
  ),
  "\n"
)

cat(
  file.path(
    table_dir,
    "deg_frequency.csv"
  ),
  "\n"
)

cat(
  file.path(
    table_dir,
    "recurrent_DEGs_2_or_more_datasets.csv"
  ),
  "\n"
)

cat(
  file.path(
    table_dir,
    "Venn_analysis_parameters.csv"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "PD_DEG_Venn_diagram.png"
  ),
  "\n"
)

cat("============================================================\n")