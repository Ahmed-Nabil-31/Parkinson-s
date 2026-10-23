# ============================================================
# DESeq2 ANALYSIS
# Parkinson's disease substantia nigra RNA-seq
#
# Datasets:
#   GSE133101
#   GSE136666
#   GSE168496
#
# Input:
#   Final tissue-filtered + haplotype-collapsed count matrices
#   Final matching metadata
#
# Analysis:
#   1. Load final count matrix
#   2. Identify biological sample columns using metadata
#   3. Verify sample matching
#   4. Verify GeneID uniqueness
#   5. Construct DESeq2 object
#   6. Filter low-count genes
#   7. Estimate surrogate variables with SVA
#   8. Include SVs when detected
#   9. Run DESeq2
#  10. Export DEG results
#  11. Generate pre/post SVA PCA plots
#
# Comparison:
#   PD vs Control
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)
library(DESeq2)
library(sva)
library(rio)


# ============================================================
# 2. PIN DPLYR VERBS
# ============================================================
#
# Prevent conflicts if another package masks dplyr verbs.
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
  "final_tissue_haplotype_collapsed"
)

output_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "DESeq2"
)

pca_output_dir <- file.path(
  base_dir,
  "results",
  "figures",
  "PCA"
)

qc_output_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "DESeq2",
  "QC"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  pca_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  qc_output_dir,
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
# 5. PCA COLORS
# ============================================================

custom_colors <- c(
  "Control" = "#4575b4",
  "PD" = "#d73027"
)


# ============================================================
# 6. PROCESS EACH DATASET
# ============================================================

for (geo_id in geo_ids) {
  
  cat("\n")
  cat("============================================================\n")
  cat("PROCESSING ", geo_id, "\n", sep = "")
  cat("============================================================\n")
  
  
  # ==========================================================
  # 6.1 FILE PATHS
  # ==========================================================
  
  count_file <- file.path(
    input_dir,
    paste0(
      geo_id,
      "_Substantia_nigra_FINAL_counts.csv"
    )
  )
  
  metadata_file <- file.path(
    input_dir,
    paste0(
      geo_id,
      "_Substantia_nigra_metadata.csv"
    )
  )
  
  
  # ==========================================================
  # 6.2 CHECK FILES
  # ==========================================================
  
  if (!file.exists(count_file)) {
    
    stop(
      "\nCount file not found:\n",
      count_file
    )
  }
  
  if (!file.exists(metadata_file)) {
    
    stop(
      "\nMetadata file not found:\n",
      metadata_file
    )
  }
  
  
  # ==========================================================
  # 6.3 LOAD COUNT DATA
  # ==========================================================
  
  counts_raw <- readr::read_csv(
    count_file,
    show_col_types = FALSE
  )
  
  
  # ==========================================================
  # 6.4 LOAD METADATA
  # ==========================================================
  
  metadata <- readr::read_csv(
    metadata_file,
    show_col_types = FALSE
  )
  
  
  cat(
    "Final count dimensions: ",
    nrow(counts_raw),
    " genes x ",
    ncol(counts_raw),
    " columns\n",
    sep = ""
  )
  
  cat(
    "Metadata dimensions: ",
    nrow(metadata),
    " samples x ",
    ncol(metadata),
    " columns\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.5 CHECK GeneID
  # ==========================================================
  
  if (!"GeneID" %in% colnames(counts_raw)) {
    
    stop(
      "\nGeneID column missing in ",
      geo_id,
      "\n"
    )
  }
  
  
  # ==========================================================
  # 6.6 CHECK DUPLICATED GeneIDs
  # ==========================================================
  
  duplicated_geneids <- counts_raw |>
    
    count(
      GeneID,
      name = "n"
    ) |>
    
    filter(
      n > 1
    )
  
  
  if (nrow(duplicated_geneids) > 0) {
    
    duplicated_file <- file.path(
      
      qc_output_dir,
      
      paste0(
        geo_id,
        "_DUPLICATED_GENEIDS_BEFORE_DESEQ2.csv"
      )
      
    )
    
    write_csv(
      duplicated_geneids,
      duplicated_file
    )
    
    stop(
      "\nDuplicated GeneIDs detected in ",
      geo_id,
      ".\nCheck:\n",
      duplicated_file
    )
  }
  
  
  cat(
    "GeneID uniqueness: PASS\n"
  )
  
  
  # ==========================================================
  # 6.7 IDENTIFY METADATA SAMPLE COLUMN
  # ==========================================================
  
  possible_sample_columns <- c(
    "sample",
    "Sample",
    "SampleID",
    "sample_id"
  )
  
  sample_column <- intersect(
    possible_sample_columns,
    colnames(metadata)
  )
  
  
  if (length(sample_column) == 0) {
    
    stop(
      "\nCould not identify sample column in metadata for ",
      geo_id,
      "\nActual columns:\n",
      paste(
        colnames(metadata),
        collapse = ", "
      )
    )
  }
  
  
  sample_column <- sample_column[1]
  
  
  metadata_samples <- as.character(
    metadata[[sample_column]]
  )
  
  
  # ==========================================================
  # 6.8 IDENTIFY ACTUAL SAMPLE COLUMNS
  # ==========================================================
  #
  # Only columns present in metadata are treated as samples.
  #
  # This prevents Gene_Symbol or other annotation columns from
  # accidentally entering the count matrix.
  # ==========================================================
  
  sample_columns <- intersect(
    colnames(counts_raw),
    metadata_samples
  )
  
  
  if (length(sample_columns) != nrow(metadata)) {
    
    stop(
      "\nSample mismatch in ",
      geo_id,
      ".\n",
      "Count sample columns detected: ",
      length(sample_columns),
      "\n",
      "Metadata samples: ",
      nrow(metadata),
      "\n"
    )
  }
  
  
  # ==========================================================
  # 6.9 CHECK FOR MISSING SAMPLES
  # ==========================================================
  
  missing_from_counts <- setdiff(
    metadata_samples,
    colnames(counts_raw)
  )
  
  
  missing_from_metadata <- setdiff(
    sample_columns,
    metadata_samples
  )
  
  
  if (
    length(missing_from_counts) > 0 ||
    length(missing_from_metadata) > 0
  ) {
    
    stop(
      "\nCount/metadata sample mismatch in ",
      geo_id,
      ".\n"
    )
  }
  
  
  cat(
    "Sample matching: PASS\n"
  )
  
  cat(
    "Biological samples used: ",
    length(sample_columns),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.10 EXTRACT COUNT MATRIX
  # ==========================================================
  
  count_data <- counts_raw |>
    
    select(
      GeneID,
      all_of(sample_columns)
    ) |>
    
    column_to_rownames(
      "GeneID"
    )
  
  
  # ==========================================================
  # 6.11 CONVERT TO MATRIX
  # ==========================================================
  
  count_data <- as.matrix(
    count_data
  )
  
  
  # ==========================================================
  # 6.12 CHECK COUNT VALUES
  # ==========================================================
  
  if (!all(is.numeric(count_data))) {
    
    stop(
      "\nCount matrix contains non-numeric values in ",
      geo_id,
      ".\n"
    )
  }
  
  
  if (
    any(
      is.na(count_data)
    )
  ) {
    
    stop(
      "\nNA values detected in count matrix for ",
      geo_id,
      ".\n"
    )
  }
  
  
  if (
    any(
      count_data < 0
    )
  ) {
    
    stop(
      "\nNegative counts detected in ",
      geo_id,
      ".\n"
    )
  }
  
  
  if (
    any(
      abs(
        count_data -
        round(count_data)
      ) > 1e-8
    )
  ) {
    
    stop(
      "\nNon-integer counts detected in ",
      geo_id,
      ".\n"
    )
  }
  
  
  mode(count_data) <- "integer"
  
  
  # ==========================================================
  # 6.13 REORDER METADATA TO COUNT MATRIX
  # ==========================================================
  
  metadata <- metadata |>
    
    filter(
      .data[[sample_column]] %in%
        sample_columns
    ) |>
    
    arrange(
      match(
        .data[[sample_column]],
        sample_columns
      )
    )
  
  
  # ==========================================================
  # 6.14 VERIFY FINAL SAMPLE ORDER
  # ==========================================================
  
  if (
    !identical(
      as.character(
        metadata[[sample_column]]
      ),
      colnames(count_data)
    )
  ) {
    
    stop(
      "\nSample order mismatch after metadata alignment in ",
      geo_id,
      ".\n"
    )
  }
  
  
  # ==========================================================
  # 6.15 IDENTIFY CONDITION COLUMN
  # ==========================================================
  
  possible_condition_columns <- c(
    "condition",
    "Condition",
    "group",
    "Group"
  )
  
  condition_column <- intersect(
    possible_condition_columns,
    colnames(metadata)
  )
  
  
  if (length(condition_column) == 0) {
    
    stop(
      "\nCondition column not found in metadata for ",
      geo_id,
      ".\nActual columns:\n",
      paste(
        colnames(metadata),
        collapse = ", "
      )
    )
  }
  
  
  condition_column <- condition_column[1]
  
  
  # ==========================================================
  # 6.16 STANDARDIZE CONDITION
  # ==========================================================
  
  raw_condition <- tolower(
    trimws(
      as.character(
        metadata[[condition_column]]
      )
    )
  )
  
  
  condition <- case_when(
    
    raw_condition %in%
      c(
        "control",
        "ctrl",
        "normal"
      ) ~
      "Control",
    
    raw_condition %in%
      c(
        "pd",
        "parkinson",
        "parkinson's",
        "parkinsons",
        "disease"
      ) ~
      "PD",
    
    TRUE ~
      NA_character_
    
  )
  
  
  if (
    any(
      is.na(condition)
    )
  ) {
    
    unknown_conditions <- unique(
      raw_condition[
        is.na(condition)
      ]
    )
    
    stop(
      "\nUnrecognized condition values in ",
      geo_id,
      ":\n",
      paste(
        unknown_conditions,
        collapse = ", "
      )
    )
  }
  
  
  # ==========================================================
  # 6.17 CREATE colData
  # ==========================================================
  
  colData <- data.frame(
    
    condition = factor(
      condition,
      levels = c(
        "Control",
        "PD"
      )
    ),
    
    row.names =
      colnames(count_data)
    
  )
  
  
  cat(
    "\nCondition distribution:\n"
  )
  
  print(
    table(
      colData$condition
    )
  )
  
  
  # ==========================================================
  # 6.18 CREATE DESeq2 OBJECT
  # ==========================================================
  
  dds <- DESeqDataSetFromMatrix(
    
    countData = count_data,
    
    colData = colData,
    
    design = ~ condition
    
  )
  
  
  # ==========================================================
  # 6.19 LOW-COUNT FILTERING
  # ==========================================================
  #
  # Keep genes with >=10 counts in at least half of the
  # samples.
  # ==========================================================
  
  min_samples <- ceiling(
    nrow(metadata) / 2
  )
  
  
  keep <- rowSums(
    counts(dds) >= 10
  ) >= min_samples
  
  
  dds <- dds[
    keep,
  ]
  
  
  cat(
    "\nGenes before low-count filtering: ",
    nrow(count_data),
    "\n",
    sep = ""
  )
  
  cat(
    "Genes after low-count filtering: ",
    nrow(dds),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.20 PRE-CORRECTION VST
  # ==========================================================
  
  dds_pre <- dds
  
  design(
    dds_pre
  ) <- ~ condition
  
  
  vsd_pre <- vst(
    dds_pre,
    blind = TRUE
  )
  
  
  # ==========================================================
  # 6.21 PCA BEFORE SVA
  # ==========================================================
  
  pca_data_pre <- plotPCA(
    
    vsd_pre,
    
    intgroup = "condition",
    
    returnData = TRUE
    
  )
  
  
  percent_pre <- round(
    
    100 *
      attr(
        pca_data_pre,
        "percentVar"
      )
    
  )
  
  
  pca_pre <- ggplot(
    
    pca_data_pre,
    
    aes(
      PC1,
      PC2,
      color = condition
    )
    
  ) +
    
    geom_point(
      size = 3
    ) +
    
    xlab(
      paste0(
        "PC1: ",
        percent_pre[1],
        "% variance"
      )
    ) +
    
    ylab(
      paste0(
        "PC2: ",
        percent_pre[2],
        "% variance"
      )
    ) +
    
    ggtitle(
      paste0(
        geo_id,
        " (Pre-SVA)"
      )
    ) +
    
    scale_color_manual(
      values = custom_colors
    ) +
    
    theme_minimal() +
    
    theme(
      
      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold",
          size = 10
        ),
      
      legend.title =
        element_text(
          face = "bold",
          size = 10
        ),
      
      panel.border =
        element_rect(
          colour = "black",
          fill = NA,
          linewidth = 1
        ),
      
      legend.position = "right",
      
      aspect.ratio = 1
      
    )
  
  
  # ==========================================================
  # 6.22 SVA
  # ==========================================================
  #
  # Estimate hidden surrogate variables while protecting
  # the biological condition of interest.
  # ==========================================================
  
  mod <- model.matrix(
    
    ~ condition,
    
    data = colData(dds)
    
  )
  
  
  mod0 <- model.matrix(
    
    ~ 1,
    
    data = colData(dds)
    
  )
  
  
  svobj <- svaseq(
    
    counts(dds),
    
    mod,
    
    mod0
    
  )
  
  
  n.sv <- svobj$n.sv
  
  
  cat(
    "\nNumber of surrogate variables detected: ",
    n.sv,
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.23 ADD SVs ONLY IF THEY EXIST
  # ==========================================================
  
  if (n.sv > 0) {
    
    for (i in seq_len(n.sv)) {
      
      colData(dds)[[
        paste0(
          "SV",
          i
        )
      ]] <-
        svobj$sv[, i]
      
    }
    
    
    sv_terms <- paste0(
      
      "SV",
      
      seq_len(n.sv),
      
      collapse = " + "
      
    )
    
    
    design(dds) <- as.formula(
      
      paste(
        "~",
        sv_terms,
        "+ condition"
      )
      
    )
    
    
    cat(
      "DESeq2 design: ~ ",
      sv_terms,
      " + condition\n",
      sep = ""
    )
    
    
  } else {
    
    design(dds) <- ~ condition
    
    cat(
      "No surrogate variables detected.\n"
    )
    
    cat(
      "DESeq2 design: ~ condition\n"
    )
  }
  
  
  # ==========================================================
  # 6.24 RUN DESeq2
  # ==========================================================
  
  dds_sva <- DESeq(
    dds
  )
  
  
  # ==========================================================
  # 6.25 EXTRACT PD VS CONTROL RESULTS
  # ==========================================================
  
  res <- results(
    
    dds_sva,
    
    contrast = c(
      "condition",
      "PD",
      "Control"
    )
    
  )
  
  
  # ==========================================================
  # 6.26 CONVERT RESULTS TO DATA FRAME
  # ==========================================================
  
  res_df <- res |>
    
    as.data.frame() |>
    
    rownames_to_column(
      "Gene_ID"
    )
  
  
  # ==========================================================
  # 6.27 EXPORT DEG RESULTS
  # ==========================================================
  
  output_file <- file.path(
    
    output_dir,
    
    paste0(
      geo_id,
      ".csv"
    )
    
  )
  
  
  write_csv(
    
    res_df,
    
    output_file
    
  )
  
  
  cat(
    "\nDESeq2 results saved:\n",
    output_file,
    "\n"
  )
  
  
  # ==========================================================
  # 6.28 POST-SVA VST
  # ==========================================================
  
  vsd_post <- vst(
    
    dds_sva,
    
    blind = FALSE
    
  )
  
  
  # ==========================================================
  # 6.29 PCA AFTER SVA
  # ==========================================================
  
  pca_data_post <- plotPCA(
    
    vsd_post,
    
    intgroup = "condition",
    
    returnData = TRUE
    
  )
  
  
  percent_post <- round(
    
    100 *
      attr(
        pca_data_post,
        "percentVar"
      )
    
  )
  
  
  pca_post <- ggplot(
    
    pca_data_post,
    
    aes(
      PC1,
      PC2,
      color = condition
    )
    
  ) +
    
    geom_point(
      size = 3
    ) +
    
    xlab(
      paste0(
        "PC1: ",
        percent_post[1],
        "% variance"
      )
    ) +
    
    ylab(
      paste0(
        "PC2: ",
        percent_post[2],
        "% variance"
      )
    ) +
    
    ggtitle(
      paste0(
        geo_id,
        " (Post-SVA)"
      )
    ) +
    
    scale_color_manual(
      values = custom_colors
    ) +
    
    theme_minimal() +
    
    theme(
      
      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold",
          size = 10
        ),
      
      legend.title =
        element_text(
          face = "bold",
          size = 10
        ),
      
      panel.border =
        element_rect(
          colour = "black",
          fill = NA,
          linewidth = 1
        ),
      
      legend.position = "right",
      
      aspect.ratio = 1
      
    )
  
  
  # ==========================================================
  # 6.30 SAVE PCA PLOTS
  # ==========================================================
  
  pca_dir <- file.path(
    
    pca_output_dir,
    
    geo_id
    
  )
  
  
  dir.create(
    
    pca_dir,
    
    recursive = TRUE,
    
    showWarnings = FALSE
    
  )
  
  
  ggsave(
    
    filename = file.path(
      pca_dir,
      paste0(
        geo_id,
        "_pre.png"
      )
    ),
    
    plot = pca_pre,
    
    width = 5,
    
    height = 5,
    
    units = "in",
    
    bg = "white",
    
    dpi = 600
    
  )
  
  
  ggsave(
    
    filename = file.path(
      pca_dir,
      paste0(
        geo_id,
        "_post.png"
      )
    ),
    
    plot = pca_post,
    
    width = 5,
    
    height = 5,
    
    units = "in",
    
    bg = "white",
    
    dpi = 600
    
  )
  
  
  # ==========================================================
  # 6.31 SAVE SVA INFORMATION
  # ==========================================================
  
  sva_summary <- tibble(
    
    Dataset = geo_id,
    
    Genes_Before_Filter =
      nrow(count_data),
    
    Genes_After_Filter =
      nrow(dds),
    
    Samples =
      ncol(count_data),
    
    Control_Samples =
      sum(
        colData$condition ==
          "Control"
      ),
    
    PD_Samples =
      sum(
        colData$condition ==
          "PD"
      ),
    
    Surrogate_Variables =
      n.sv,
    
    DESeq2_Design =
      ifelse(
        n.sv > 0,
        paste0(
          "~ SV1",
          ifelse(
            n.sv > 1,
            paste0(
              " + SV2...SV",
              n.sv
            ),
            ""
          ),
          " + condition"
        ),
        "~ condition"
      )
    
  )
  
  
  write_csv(
    
    sva_summary,
    
    file.path(
      
      qc_output_dir,
      
      paste0(
        geo_id,
        "_DESeq2_SVA_SUMMARY.csv"
      )
      
    )
    
  )
  
  
  # ==========================================================
  # 6.32 FINISHED DATASET
  # ==========================================================
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("FINISHED: ", geo_id, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  
}


# ============================================================
# END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ALL DESeq2 ANALYSES COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")