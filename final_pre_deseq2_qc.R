# ============================================================
# FINAL PRE-DESeq2 QC
# Parkinson's disease RNA-seq
#
# Purpose:
# Validate the FINAL tissue-filtered + haplotype-collapsed
# count matrices immediately before DESeq2.
#
# IMPORTANT:
# Sample columns are identified by MATCHING THEM TO THE
# SAMPLE IDs PRESENT IN THE METADATA.
#
# This prevents annotation/non-sample columns from being
# incorrectly counted as biological samples.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(tidyverse)


# ============================================================
# 2. DIRECTORIES
# ============================================================

base_dir <- "/media/nabil/ssd/Parkinson-s"

final_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "final_tissue_haplotype_collapsed"
)

qc_dir <- file.path(
  final_dir,
  "pre_deseq2_qc"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. DATASETS
# ============================================================

datasets <- c(
  "GSE133101",
  "GSE136666",
  "GSE168496"
)


# Expected Substantia nigra sample numbers

expected_samples <- c(
  GSE133101 = 25,
  GSE136666 = 10,
  GSE168496 = 16
)


# ============================================================
# 4. HELPER FUNCTION
# ============================================================

find_column <- function(
    data,
    candidates,
    dataset,
    purpose
) {
  
  found <- candidates[
    candidates %in% colnames(data)
  ]
  
  if (length(found) == 0) {
    
    stop(
      "\nCould not identify ",
      purpose,
      " column in ",
      dataset,
      ".\n\nAvailable columns:\n",
      paste(
        colnames(data),
        collapse = ", "
      ),
      "\n"
    )
  }
  
  found[1]
}


# ============================================================
# 5. PROCESS DATASETS
# ============================================================

all_results <- list()


for (gse in datasets) {
  
  
  cat("\n")
  cat("============================================================\n")
  cat("FINAL PRE-DESeq2 QC: ", gse, "\n", sep = "")
  cat("============================================================\n")
  
  
  # ==========================================================
  # 5.1 FILE PATHS
  # ==========================================================
  
  count_file <- file.path(
    final_dir,
    paste0(
      gse,
      "_Substantia_nigra_FINAL_counts.csv"
    )
  )
  
  metadata_file <- file.path(
    final_dir,
    paste0(
      gse,
      "_Substantia_nigra_metadata.csv"
    )
  )
  
  
  # ==========================================================
  # 5.2 CHECK FILES
  # ==========================================================
  
  if (!file.exists(count_file)) {
    
    stop(
      "\nFinal count file not found:\n",
      count_file
    )
  }
  
  if (!file.exists(metadata_file)) {
    
    stop(
      "\nFinal metadata file not found:\n",
      metadata_file
    )
  }
  
  
  # ==========================================================
  # 5.3 READ FILES
  # ==========================================================
  
  counts <- readr::read_csv(
    count_file,
    show_col_types = FALSE
  )
  
  metadata <- readr::read_csv(
    metadata_file,
    show_col_types = FALSE
  )
  
  
  cat(
    "\nCount matrix dimensions: ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " columns\n",
    sep = ""
  )
  
  cat(
    "Metadata dimensions: ",
    nrow(metadata),
    " rows x ",
    ncol(metadata),
    " columns\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.4 CHECK GeneID
  # ==========================================================
  
  if (!"GeneID" %in% colnames(counts)) {
    
    stop(
      "\nGeneID column not found in ",
      gse,
      "\nActual columns:\n",
      paste(
        colnames(counts),
        collapse = ", "
      )
    )
  }
  
  
  # ==========================================================
  # 5.5 FIND METADATA SAMPLE COLUMN
  # ==========================================================
  
  metadata_sample_column <- find_column(
    
    metadata,
    
    c(
      "sample",
      "Sample",
      "SampleID",
      "sample_id"
    ),
    
    gse,
    
    "sample"
  )
  
  
  metadata_samples <- as.character(
    metadata[[metadata_sample_column]]
  )
  
  
  metadata_samples <- unique(
    metadata_samples
  )
  
  
  # ==========================================================
  # 5.6 IDENTIFY ACTUAL COUNT SAMPLE COLUMNS
  #
  # THIS IS THE IMPORTANT FIX.
  #
  # Only columns whose names occur in metadata are considered
  # biological samples.
  # ==========================================================
  
  actual_sample_columns <- intersect(
    
    colnames(counts),
    
    metadata_samples
    
  )
  
  
  # ==========================================================
  # 5.7 IDENTIFY NON-SAMPLE COLUMNS
  # ==========================================================
  
  non_sample_columns <- setdiff(
    
    colnames(counts),
    
    actual_sample_columns
    
  )
  
  
  cat(
    "\nActual sample columns detected: ",
    length(actual_sample_columns),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Expected samples: ",
    expected_samples[[gse]],
    "\n",
    sep = ""
  )
  
  
  cat(
    "Non-sample columns detected: ",
    length(non_sample_columns),
    "\n",
    sep = ""
  )
  
  
  if (length(non_sample_columns) > 0) {
    
    cat(
      "\nNon-sample columns:\n"
    )
    
    print(
      non_sample_columns
    )
  }
  
  
  # ==========================================================
  # 5.8 SAMPLE NUMBER CHECK
  # ==========================================================
  
  sample_number_ok <-
    
    length(actual_sample_columns) ==
    
    expected_samples[[gse]]
  
  
  cat(
    "\nSample number check: ",
    ifelse(
      sample_number_ok,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.9 COUNT / METADATA SAMPLE MATCHING
  # ==========================================================
  
  missing_from_counts <- setdiff(
    
    metadata_samples,
    
    colnames(counts)
    
  )
  
  
  missing_from_metadata <- setdiff(
    
    actual_sample_columns,
    
    metadata_samples
    
  )
  
  
  samples_match <-
    
    length(missing_from_counts) == 0 &&
    
    length(missing_from_metadata) == 0
  
  
  cat(
    "Count/metadata sample matching: ",
    ifelse(
      samples_match,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  
  if (length(missing_from_counts) > 0) {
    
    cat(
      "\nMetadata samples missing from count matrix:\n"
    )
    
    print(
      missing_from_counts
    )
  }
  
  
  if (length(missing_from_metadata) > 0) {
    
    cat(
      "\nCount samples missing from metadata:\n"
    )
    
    print(
      missing_from_metadata
    )
  }
  
  
  # ==========================================================
  # 5.10 SAMPLE ORDER CHECK
  # ==========================================================
  #
  # Metadata order is compared with the count matrix sample
  # order, using only actual sample columns.
  # ==========================================================
  
  metadata_samples_ordered <-
    as.character(
      metadata[[metadata_sample_column]]
    )
  
  
  sample_order_identical <-
    
    identical(
      actual_sample_columns,
      metadata_samples_ordered
    )
  
  
  cat(
    "Sample order identical: ",
    ifelse(
      sample_order_identical,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.11 DUPLICATED GeneIDs
  # ==========================================================
  
  duplicated_geneids <- counts |>
    
    count(
      GeneID,
      name = "n"
    ) |>
    
    filter(
      n > 1
    )
  
  
  n_duplicate_geneids <-
    nrow(
      duplicated_geneids
    )
  
  
  duplicated_geneids_pass <-
    
    n_duplicate_geneids == 0
  
  
  cat(
    "Duplicated GeneIDs: ",
    ifelse(
      duplicated_geneids_pass,
      "PASS",
      "FAIL"
    ),
    " (",
    n_duplicate_geneids,
    ")\n",
    sep = ""
  )
  
  
  if (n_duplicate_geneids > 0) {
    
    write_csv(
      
      duplicated_geneids,
      
      file.path(
        qc_dir,
        paste0(
          gse,
          "_DUPLICATED_GENEIDS.csv"
        )
      )
    )
  }
  
  
  # ==========================================================
  # 5.12 EXTRACT ONLY ACTUAL COUNT MATRIX
  # ==========================================================
  
  count_matrix <- counts |>
    
    select(
      all_of(
        actual_sample_columns
      )
    )
  
  
  # ==========================================================
  # 5.13 NUMERIC COUNTS
  # ==========================================================
  
  numeric_counts <- all(
    
    vapply(
      count_matrix,
      is.numeric,
      logical(1)
    )
    
  )
  
  
  cat(
    "Numeric count columns: ",
    ifelse(
      numeric_counts,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.14 NA VALUES
  # ==========================================================
  
  na_count_values <- any(
    
    is.na(
      as.matrix(
        count_matrix
      )
    )
    
  )
  
  
  cat(
    "NA count values: ",
    ifelse(
      !na_count_values,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.15 NEGATIVE COUNTS
  # ==========================================================
  
  if (numeric_counts) {
    
    negative_count_values <- any(
      
      as.matrix(
        count_matrix
      ) < 0,
      
      na.rm = TRUE
      
    )
    
  } else {
    
    negative_count_values <- NA
  }
  
  
  if (!is.na(negative_count_values)) {
    
    cat(
      "Negative count values: ",
      ifelse(
        !negative_count_values,
        "PASS",
        "FAIL"
      ),
      "\n",
      sep = ""
    )
    
  } else {
    
    cat(
      "Negative count values: NOT EVALUATED\n"
    )
  }
  
  
  # ==========================================================
  # 5.16 INTEGER COUNTS
  # ==========================================================
  
  if (numeric_counts) {
    
    count_values <- as.matrix(
      count_matrix
    )
    
    non_integer_values <- any(
      
      abs(
        count_values -
          round(count_values)
      ) > 1e-8,
      
      na.rm = TRUE
      
    )
    
  } else {
    
    non_integer_values <- NA
  }
  
  
  if (!is.na(non_integer_values)) {
    
    cat(
      "Integer count values: ",
      ifelse(
        !non_integer_values,
        "PASS",
        "FAIL"
      ),
      "\n",
      sep = ""
    )
    
  } else {
    
    cat(
      "Integer count values: NOT EVALUATED\n"
    )
  }
  
  
  # ==========================================================
  # 5.17 ZERO-COUNT GENES
  # ==========================================================
  
  if (numeric_counts) {
    
    gene_totals <- rowSums(
      
      as.matrix(
        count_matrix
      ),
      
      na.rm = TRUE
      
    )
    
    zero_count_genes <-
      sum(
        gene_totals == 0
      )
    
  } else {
    
    zero_count_genes <- NA
  }
  
  
  cat(
    "Zero-count genes: ",
    ifelse(
      is.na(zero_count_genes),
      "NOT EVALUATED",
      zero_count_genes
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.18 ZERO-COUNT SAMPLES
  # ==========================================================
  
  if (numeric_counts) {
    
    sample_totals <- colSums(
      
      as.matrix(
        count_matrix
      ),
      
      na.rm = TRUE
      
    )
    
    zero_count_samples <- sum(
      
      sample_totals == 0
      
    )
    
  } else {
    
    zero_count_samples <- NA
  }
  
  
  zero_sample_pass <-
    
    !is.na(zero_count_samples) &&
    
    zero_count_samples == 0
  
  
  cat(
    "Zero-count samples: ",
    ifelse(
      is.na(zero_count_samples),
      "NOT EVALUATED",
      ifelse(
        zero_sample_pass,
        "PASS (0)",
        paste0(
          "FAIL (",
          zero_count_samples,
          ")"
        )
      )
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.19 CONDITION COLUMN
  # ==========================================================
  
  condition_column <- find_column(
    
    metadata,
    
    c(
      "condition",
      "Condition",
      "group",
      "Group"
    ),
    
    gse,
    
    "condition"
  )
  
  
  condition_values <- as.character(
    
    metadata[[condition_column]]
    
  )
  
  
  condition_levels <- unique(
    condition_values
  )
  
  
  cat(
    "\nCondition column: ",
    condition_column,
    "\n",
    sep = ""
  )
  
  
  cat(
    "Condition levels:\n"
  )
  
  print(
    table(
      condition_values
    )
  )
  
  
  condition_pass <-
    
    length(
      condition_levels
    ) >= 2
  
  
  cat(
    "Condition information: ",
    ifelse(
      condition_pass,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 5.20 SAVE DETECTED SAMPLE INFORMATION
  # ==========================================================
  
  sample_audit <- tibble(
    
    Dataset = gse,
    
    Count_Column =
      colnames(counts),
    
    Is_GeneID =
      colnames(counts) == "GeneID",
    
    Is_Metadata_Sample =
      colnames(counts) %in%
      metadata_samples,
    
    Classification = case_when(
      
      colnames(counts) == "GeneID" ~
        "GENE_ID",
      
      colnames(counts) %in%
        metadata_samples ~
        "SAMPLE",
      
      TRUE ~
        "NON_SAMPLE_COLUMN"
      
    )
    
  )
  
  
  write_csv(
    
    sample_audit,
    
    file.path(
      qc_dir,
      paste0(
        gse,
        "_SAMPLE_COLUMN_AUDIT.csv"
      )
    )
    
  )
  
  
  # ==========================================================
  # 5.21 FINAL SUMMARY
  # ==========================================================
  
  result <- tibble(
    
    Dataset = gse,
    
    Final_Genes =
      nrow(counts),
    
    Total_Count_File_Columns =
      ncol(counts),
    
    Detected_Sample_Columns =
      length(actual_sample_columns),
    
    Expected_Samples =
      expected_samples[[gse]],
    
    Sample_Number =
      ifelse(
        sample_number_ok,
        "PASS",
        "FAIL"
      ),
    
    Non_Sample_Columns =
      length(non_sample_columns),
    
    Sample_Metadata_Match =
      ifelse(
        samples_match,
        "PASS",
        "FAIL"
      ),
    
    Sample_Order =
      ifelse(
        sample_order_identical,
        "PASS",
        "FAIL"
      ),
    
    Duplicated_GeneIDs =
      ifelse(
        duplicated_geneids_pass,
        "PASS",
        "FAIL"
      ),
    
    Numeric_Counts =
      ifelse(
        numeric_counts,
        "PASS",
        "FAIL"
      ),
    
    NA_Counts =
      ifelse(
        !na_count_values,
        "PASS",
        "FAIL"
      ),
    
    Negative_Counts =
      ifelse(
        is.na(negative_count_values),
        "NOT_EVALUATED",
        ifelse(
          !negative_count_values,
          "PASS",
          "FAIL"
        )
      ),
    
    Integer_Counts =
      ifelse(
        is.na(non_integer_values),
        "NOT_EVALUATED",
        ifelse(
          !non_integer_values,
          "PASS",
          "FAIL"
        )
      ),
    
    Zero_Count_Genes =
      zero_count_genes,
    
    Zero_Count_Samples =
      zero_count_samples,
    
    Condition_Information =
      ifelse(
        condition_pass,
        "PASS",
        "FAIL"
      )
    
  )
  
  
  all_results[[gse]] <- result
  
  
  # ==========================================================
  # 5.22 PRINT SUMMARY
  # ==========================================================
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("QC SUMMARY: ", gse, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  
  print(result)
  
}


# ============================================================
# 6. COMBINE DATASETS
# ============================================================

final_summary <- bind_rows(
  all_results
)


# ============================================================
# 7. SAVE MASTER QC SUMMARY
# ============================================================

summary_file <- file.path(
  
  qc_dir,
  
  "ALL_DATASETS_PRE_DESeq2_QC_SUMMARY.csv"
  
)


write_csv(
  
  final_summary,
  
  summary_file
  
)


# ============================================================
# 8. FINAL REPORT
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("ALL PRE-DESeq2 QC COMPLETED\n")
cat("============================================================\n\n")

print(
  final_summary
)

cat(
  "\nMaster QC summary saved to:\n",
  summary_file,
  "\n"
)

cat(
  "\nIndividual sample-column audits saved in:\n",
  qc_dir,
  "\n"
)