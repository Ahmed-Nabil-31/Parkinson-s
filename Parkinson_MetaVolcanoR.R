# ============================================================
# Parkinson's Disease RNA-seq Meta-Analysis
#
# Datasets:
#   GSE133101
#   GSE136666
#   GSE168496
#
# Input:
#   Final DESeq2 results
#
# Main analysis:
#   MetaVolcanoR random-effects model
#
# Secondary analysis:
#   Mean fold-change combining approach
#
# Gene identifier:
#   Ensembl Gene ID
#
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(MetaVolcanoR)
library(tidyverse)
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
# 3. Define project directories
# ------------------------------------------------------------

base_dir <- "/media/nabil/ssd/Parkinson-s"

deseq_dir <- file.path(
  base_dir,
  "results/tables/DESeq2"
)

annotation_dir <- file.path(
  base_dir,
  "results/tables/annotation_audit"
)

output_table_dir <- file.path(
  base_dir,
  "results/tables/meta-analysis"
)

output_figure_dir <- file.path(
  base_dir,
  "results/figures/meta-analysis"
)


dir.create(
  output_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Define Parkinson's datasets
# ------------------------------------------------------------

geo_ids <- c(
  "GSE133101",
  "GSE136666",
  "GSE168496"
)


# ------------------------------------------------------------
# 5. Check input files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CHECKING DESeq2 INPUT FILES\n")
cat("============================================================\n\n")


for (gse in geo_ids) {
  
  f <- file.path(
    deseq_dir,
    paste0(gse, ".csv")
  )
  
  if (!file.exists(f)) {
    
    stop(
      "Missing DESeq2 file:\n",
      f
    )
    
  }
  
  cat(
    "FOUND: ",
    f,
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 6. Load existing gene annotation
#
# We do NOT query Ensembl again.
# The annotation audit files already contain:
#   GeneID
#   Gene_Symbol
#   Gene_Biotype
#   Chromosome
#   Start
#   End
#   Strand
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("LOADING EXISTING ANNOTATION\n")
cat("============================================================\n\n")


annotation_files <- file.path(
  annotation_dir,
  paste0(
    geo_ids,
    "_Ensembl_annotation.csv"
  )
)


missing_annotation <- annotation_files[
  !file.exists(annotation_files)
]


if (length(missing_annotation) > 0) {
  
  stop(
    "Missing annotation file(s):\n",
    paste(
      missing_annotation,
      collapse = "\n"
    )
  )
  
}


gene_annotation <- lapply(
  annotation_files,
  function(f) {
    
    import(f) |>
      
      select(
        any_of(
          c(
            "GeneID",
            "Gene_Symbol",
            "Gene_Biotype",
            "Chromosome",
            "Start",
            "End",
            "Strand"
          )
        )
        
      ) |>
      
      mutate(
        GeneID = sub(
          "\\..*$",
          "",
          as.character(GeneID)
        )
      )
    
  }
  
) |>
  
  bind_rows() |>
  
  distinct(
    GeneID,
    .keep_all = TRUE
  )


cat(
  "Unique annotated Ensembl IDs: ",
  nrow(gene_annotation),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Load DESeq2 results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("LOADING DESeq2 RESULTS\n")
cat("============================================================\n\n")


studies <- list()


for (gse in geo_ids) {
  
  f <- file.path(
    deseq_dir,
    paste0(gse, ".csv")
  )
  
  df <- import(f)
  
  
  # ----------------------------------------------------------
  # Check required columns
  # ----------------------------------------------------------
  
  required_columns <- c(
    "Gene_ID",
    "log2FoldChange",
    "lfcSE",
    "padj"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      "\nDataset ",
      gse,
      " is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Standardize Ensembl IDs
  # ----------------------------------------------------------
  
  df <- df |>
    
    mutate(
      
      Gene_ID = sub(
        "\\..*$",
        "",
        as.character(Gene_ID)
      ),
      
      log2FoldChange =
        as.numeric(log2FoldChange),
      
      lfcSE =
        as.numeric(lfcSE),
      
      padj =
        as.numeric(padj)
      
    )
  
  
  # ----------------------------------------------------------
  # Keep only information required by MetaVolcanoR
  # ----------------------------------------------------------
  
  df <- df |>
    
    select(
      Gene_ID,
      log2FoldChange,
      lfcSE,
      padj
    )
  
  
  # ----------------------------------------------------------
  # Remove unusable rows
  # ----------------------------------------------------------
  
  df <- df |>
    
    filter(
      
      !is.na(Gene_ID),
      
      Gene_ID != "",
      
      !is.na(log2FoldChange),
      
      !is.na(lfcSE),
      
      lfcSE > 0,
      
      !is.na(padj)
      
    )
  
  
  # ----------------------------------------------------------
  # Remove duplicate Gene IDs if any remain
  # ----------------------------------------------------------
  
  if (anyDuplicated(df$Gene_ID) > 0) {
    
    warning(
      gse,
      ": duplicated Gene_IDs detected in DESeq2 output. ",
      "Keeping the first occurrence."
    )
    
    df <- df |>
      distinct(
        Gene_ID,
        .keep_all = TRUE
      )
    
  }
  
  
  studies[[gse]] <- df
  
  
  cat(
    gse,
    ": ",
    nrow(df),
    " genes loaded\n",
    sep = ""
  )
  
}


names(studies) <- geo_ids


# ------------------------------------------------------------
# 8. Confirm common gene universe
# ------------------------------------------------------------

common_gene_universe <- Reduce(
  intersect,
  lapply(
    studies,
    function(x) x$Gene_ID
  )
)


cat(
  "\nGenes present in ALL THREE datasets: ",
  length(common_gene_universe),
  "\n",
  sep = ""
)


if (length(common_gene_universe) == 0) {
  
  stop(
    "No genes are shared across all three DESeq2 datasets."
  )
  
}


# ------------------------------------------------------------
# 9. Restrict each study to the common gene universe
#
# This makes the meta-analysis operate on the same
# Ensembl gene set across all three studies.
# ------------------------------------------------------------

studies <- lapply(
  studies,
  function(df) {
    
    df |>
      
      filter(
        Gene_ID %in% common_gene_universe
      )
    
  }
)


# ------------------------------------------------------------
# 10. Meta-analysis: Random Effects Model
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RUNNING RANDOM-EFFECTS META-ANALYSIS\n")
cat("============================================================\n\n")


meta_degs_rem <- rem_mv(
  
  diffexp =
    studies,
  
  pcriteria =
    "padj",
  
  foldchangecol =
    "log2FoldChange",
  
  genenamecol =
    "Gene_ID",
  
  geneidcol =
    NULL,
  
  collaps =
    TRUE,
  
  vcol =
    "lfcSE",
  
  cvar =
    FALSE,
  
  metathr =
    0.01,
  
  jobname =
    "Parkinson_MetaVolcano_REM",
  
  outputfolder =
    paste0(
      output_figure_dir,
      "/"
    ),
  
  draw =
    "PDF",
  
  ncores =
    1
)


# ------------------------------------------------------------
# 11. Extract random-effects results
# ------------------------------------------------------------

meta_results <- meta_degs_rem@metaresult


# ------------------------------------------------------------
# 12. Annotate meta-analysis results
# ------------------------------------------------------------

meta_results <- meta_results |>
  
  mutate(
    Gene_ID = sub(
      "\\..*$",
      "",
      as.character(Gene_ID)
    )
  ) |>
  
  left_join(
    gene_annotation,
    by = c(
      "Gene_ID" = "GeneID"
    )
  )


# ------------------------------------------------------------
# 13. Add direction of effect
#
# randomSummary is the pooled/meta-analysis effect.
# Positive = higher in Parkinson's condition
# Negative = lower in Parkinson's condition
# ------------------------------------------------------------

meta_results <- meta_results |>
  
  mutate(
    
    Meta_Direction = case_when(
      
      randomSummary > 0 ~
        "Up",
      
      randomSummary < 0 ~
        "Down",
      
      TRUE ~
        NA_character_
      
    )
    
  )


# ------------------------------------------------------------
# 14. Add consistency of direction across datasets
# ------------------------------------------------------------

direction_table <- bind_rows(
  
  lapply(
    names(studies),
    function(gse) {
      
      studies[[gse]] |>
        
        select(
          Gene_ID,
          log2FoldChange
        ) |>
        
        mutate(
          Dataset = gse
        )
      
    }
  )
  
) |>
  
  mutate(
    
    Direction = case_when(
      
      log2FoldChange > 0 ~
        "Up",
      
      log2FoldChange < 0 ~
        "Down",
      
      TRUE ~
        NA_character_
      
    )
    
  )


direction_summary <- direction_table |>
  
  group_by(Gene_ID) |>
  
  summarise(
    
    Number_Up =
      sum(Direction == "Up", na.rm = TRUE),
    
    Number_Down =
      sum(Direction == "Down", na.rm = TRUE),
    
    Number_Datasets =
      n(),
    
    Direction_Consistency =
      case_when(
        
        Number_Up == Number_Datasets ~
          "Consistently Up",
        
        Number_Down == Number_Datasets ~
          "Consistently Down",
        
        TRUE ~
          "Discordant"
        
      ),
    
    .groups = "drop"
    
  )


meta_results <- meta_results |>
  
  left_join(
    direction_summary,
    by = "Gene_ID"
  )


# ------------------------------------------------------------
# 15. Save COMPLETE random-effects result
# ------------------------------------------------------------

random_effect_file <- file.path(
  
  output_table_dir,
  
  "Parkinson_random_effects_meta_analysis.csv"
  
)


export(
  meta_results,
  random_effect_file
)


# ------------------------------------------------------------
# 16. Filter significant meta-analysis DEGs
#
# Criteria:
#
#   randomP < 0.05
#   |randomSummary| >= 1
#
# ------------------------------------------------------------

filtered_meta <- meta_results |>
  
  filter(
    
    !is.na(randomP),
    
    randomP < 0.05,
    
    !is.na(randomSummary),
    
    abs(randomSummary) >= 1
    
  ) |>
  
  arrange(
    randomP
  )


filtered_meta_file <- file.path(
  
  output_table_dir,
  
  "Parkinson_filtered_meta_DEGs.csv"
  
)


export(
  filtered_meta,
  filtered_meta_file
)


# ------------------------------------------------------------
# 17. Save consistently directed meta-DEGs
# ------------------------------------------------------------

consistent_meta <- filtered_meta |>
  
  filter(
    Direction_Consistency %in%
      c(
        "Consistently Up",
        "Consistently Down"
      )
  )


consistent_meta_file <- file.path(
  
  output_table_dir,
  
  "Parkinson_consistent_meta_DEGs.csv"
  
)


export(
  consistent_meta,
  consistent_meta_file
)


# ------------------------------------------------------------
# 18. Mean fold-change combining analysis
#
# This is a SECONDARY analysis.
# It is not the primary result.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RUNNING MEAN FOLD-CHANGE COMBINING ANALYSIS\n")
cat("============================================================\n\n")


meta_degs_comb <- combining_mv(
  
  diffexp =
    studies,
  
  pcriteria =
    "padj",
  
  foldchangecol =
    "log2FoldChange",
  
  genenamecol =
    "Gene_ID",
  
  metafc =
    "Mean",
  
  metathr =
    0.01,
  
  collaps =
    TRUE,
  
  jobname =
    "Parkinson_MetaVolcano_Combining",
  
  outputfolder =
    paste0(
      output_figure_dir,
      "/"
    ),
  
  draw =
    "PDF"
  
)


combined_meta <- meta_degs_comb@metaresult


# ------------------------------------------------------------
# 19. Annotate combining analysis
# ------------------------------------------------------------

combined_meta <- combined_meta |>
  
  mutate(
    
    Gene_ID = sub(
      "\\..*$",
      "",
      as.character(Gene_ID)
    )
    
  ) |>
  
  left_join(
    
    gene_annotation,
    
    by = c(
      "Gene_ID" = "GeneID"
    )
    
  )


# ------------------------------------------------------------
# 20. Save combining analysis
# ------------------------------------------------------------

combined_file <- file.path(
  
  output_table_dir,
  
  "Parkinson_mean_fold_change_meta_analysis.csv"
  
)


export(
  combined_meta,
  combined_file
)


# ------------------------------------------------------------
# 21. Print final summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PARKINSON'S META-ANALYSIS COMPLETED\n")
cat("============================================================\n\n")


cat(
  "Datasets included:\n"
)

for (gse in geo_ids) {
  
  cat(
    "  - ",
    gse,
    "\n",
    sep = ""
  )
  
}


cat(
  "\nCommon gene universe: ",
  length(common_gene_universe),
  "\n",
  sep = ""
)


cat(
  "\nRandom-effects meta-analysis:\n"
)

cat(
  "  Total genes analyzed: ",
  nrow(meta_results),
  "\n",
  sep = ""
)

cat(
  "  Significant meta-DEGs: ",
  nrow(filtered_meta),
  "\n",
  sep = ""
)

cat(
  "  Consistently directed meta-DEGs: ",
  nrow(consistent_meta),
  "\n",
  sep = ""
)


cat(
  "\nOutput files:\n\n"
)

cat(
  "1. ",
  random_effect_file,
  "\n"
)

cat(
  "2. ",
  filtered_meta_file,
  "\n"
)

cat(
  "3. ",
  consistent_meta_file,
  "\n"
)

cat(
  "4. ",
  combined_file,
  "\n"
)


cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")