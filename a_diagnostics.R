# ============================================================
# Parkinson's Disease RNA-seq
# Raw Count Matrix + Metadata Diagnostic Inspection
# ============================================================
#
# Purpose:
#   1. Load raw count matrices for all three datasets
#   2. Load corresponding metadata
#   3. Inspect Ensembl gene IDs
#   4. Check duplicated Ensembl IDs
#   5. Inspect gene-symbol annotation
#   6. Check duplicated gene symbols
#   7. Check sample IDs between counts and metadata
#   8. Inspect tissue and condition information
#   9. Save one diagnostic summary for each dataset
#
# IMPORTANT:
#   We are NOT collapsing gene symbols yet.
#   This script is diagnostic only.
#
# Datasets:
#   GSE133101
#   GSE136666
#   GSE168496
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(tidyverse)
library(rio)


# ------------------------------------------------------------
# 2. Define directories
# ------------------------------------------------------------

raw_dir <- "/media/nabil/ssd/Parkinson-s/unedited_data/raw_counts"

metadata_dir <- "/media/nabil/ssd/Parkinson-s/unedited_data/metadata"

output_dir <- "/media/nabil/ssd/Parkinson-s/results/tables/preprocessing"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# 3. Define datasets
# ------------------------------------------------------------

geo_ids <- c(
  "GSE133101",
  "GSE136666",
  "GSE168496"
)


# ============================================================
# 4. Process each dataset
# ============================================================

for (geo_id in geo_ids) {
  
  cat("\n\n")
  cat("============================================================\n")
  cat(geo_id, "\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # 4.1 Define file paths
  # ----------------------------------------------------------
  
  count_file <- file.path(
    raw_dir,
    paste0(geo_id, "_raw_counts.csv")
  )
  
  metadata_file <- file.path(
    metadata_dir,
    paste0(geo_id, "_metadata.csv")
  )
  
  
  # ----------------------------------------------------------
  # 4.2 Check files exist
  # ----------------------------------------------------------
  
  if (!file.exists(count_file)) {
    
    cat("WARNING: Count file not found:\n")
    cat(count_file, "\n")
    
    next
  }
  
  if (!file.exists(metadata_file)) {
    
    cat("WARNING: Metadata file not found:\n")
    cat(metadata_file, "\n")
    
    next
  }
  
  
  # ----------------------------------------------------------
  # 4.3 Load count matrix and metadata
  # ----------------------------------------------------------
  
  counts <- import(count_file)
  
  metadata <- import(metadata_file)
  
  
  # ==========================================================
  # 5. BASIC DATASET STRUCTURE
  # ==========================================================
  
  cat("\n--- Basic dimensions ---\n")
  
  cat(
    "Count matrix:",
    nrow(counts),
    "rows x",
    ncol(counts),
    "columns\n"
  )
  
  cat(
    "Metadata:",
    nrow(metadata),
    "rows x",
    ncol(metadata),
    "columns\n"
  )
  
  
  # ----------------------------------------------------------
  # 5.1 Count-matrix columns
  # ----------------------------------------------------------
  
  cat("\nCount-matrix columns:\n")
  
  print(head(colnames(counts), 15))
  
  
  # ----------------------------------------------------------
  # 5.2 Metadata columns
  # ----------------------------------------------------------
  
  cat("\nMetadata columns:\n")
  
  print(colnames(metadata))
  
  
  # ==========================================================
  # 6. GENE ID INSPECTION
  # ==========================================================
  
  gene_id_column <- "GeneID"
  
  if (!gene_id_column %in% colnames(counts)) {
    
    cat(
      "\nWARNING: GeneID column not found.\n"
    )
    
    next
  }
  
  gene_ids <- as.character(
    counts[[gene_id_column]]
  )
  
  
  # ----------------------------------------------------------
  # 6.1 First gene IDs
  # ----------------------------------------------------------
  
  cat("\nFirst 20 gene IDs:\n")
  
  print(head(gene_ids, 20))
  
  
  # ----------------------------------------------------------
  # 6.2 Check Ensembl format
  # ----------------------------------------------------------
  
  ensembl_no_version <- grepl(
    "^ENSG[0-9]+$",
    gene_ids
  )
  
  ensembl_with_version <- grepl(
    "^ENSG[0-9]+\\.[0-9]+$",
    gene_ids
  )
  
  
  cat(
    "\nEnsembl IDs without version:",
    sum(ensembl_no_version),
    "\n"
  )
  
  cat(
    "Ensembl IDs with version:",
    sum(ensembl_with_version),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 6.3 Duplicate Ensembl IDs
  # ----------------------------------------------------------
  
  duplicated_gene_ids <- unique(
    gene_ids[
      duplicated(gene_ids)
    ]
  )
  
  
  cat(
    "Duplicated Ensembl IDs:",
    length(duplicated_gene_ids),
    "\n"
  )
  
  
  if (length(duplicated_gene_ids) > 0) {
    
    cat("\nFirst duplicated Ensembl IDs:\n")
    
    print(
      head(duplicated_gene_ids, 20)
    )
  }
  
  
  # ----------------------------------------------------------
  # 6.4 Unique Ensembl IDs
  # ----------------------------------------------------------
  
  cat(
    "Unique Ensembl IDs:",
    length(unique(gene_ids)),
    "\n"
  )
  
  
  # ==========================================================
  # 7. GENE SYMBOL INSPECTION
  # ==========================================================
  
  symbol_column <- "gene_symbol"
  
  if (!symbol_column %in% colnames(counts)) {
    
    cat(
      "\nWARNING: gene_symbol column not found.\n"
    )
    
    symbol_available <- FALSE
    
  } else {
    
    symbol_available <- TRUE
    
    symbols <- as.character(
      counts[[symbol_column]]
    )
    
    
    # --------------------------------------------------------
    # 7.1 Missing symbols
    # --------------------------------------------------------
    
    missing_symbols <- (
      is.na(symbols) |
        symbols == ""
    )
    
    
    cat(
      "\nRows with gene symbols:",
      sum(!missing_symbols),
      "\n"
    )
    
    cat(
      "Blank/NA gene symbols:",
      sum(missing_symbols),
      "\n"
    )
    
    
    # --------------------------------------------------------
    # 7.2 Unique gene symbols
    # --------------------------------------------------------
    
    symbols_clean <- symbols[
      !missing_symbols
    ]
    
    cat(
      "Unique gene symbols:",
      length(unique(symbols_clean)),
      "\n"
    )
    
    
    # --------------------------------------------------------
    # 7.3 Duplicate gene symbols
    # --------------------------------------------------------
    
    duplicated_symbols <- unique(
      symbols_clean[
        duplicated(symbols_clean)
      ]
    )
    
    
    cat(
      "Duplicated gene symbols:",
      length(duplicated_symbols),
      "\n"
    )
    
    
    # --------------------------------------------------------
    # 7.4 Show first duplicated symbols
    # --------------------------------------------------------
    
    if (length(duplicated_symbols) > 0) {
      
      cat(
        "\nFirst 20 duplicated gene symbols:\n"
      )
      
      print(
        head(duplicated_symbols, 20)
      )
    }
  }
  
  
  # ==========================================================
  # 8. SAMPLE ID INSPECTION
  # ==========================================================
  
  # Count matrix:
  #   column 1 = GeneID
  #   column 2 = gene_symbol
  #   remaining columns = samples
  
  count_sample_ids <- colnames(counts)[
    !colnames(counts) %in%
      c("GeneID", "gene_symbol")
  ]
  
  
  cat(
    "\nNumber of count-matrix samples:",
    length(count_sample_ids),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 8.1 Identify metadata sample column
  # ----------------------------------------------------------
  
  possible_sample_columns <- colnames(metadata)[
    grepl(
      "sample|gsm|accession",
      colnames(metadata),
      ignore.case = TRUE
    )
  ]
  
  
  cat("\nPossible metadata sample-ID columns:\n")
  
  print(possible_sample_columns)
  
  
  # ----------------------------------------------------------
  # 8.2 Compare sample IDs
  # ----------------------------------------------------------
  
  if ("sample" %in% colnames(metadata)) {
    
    metadata_sample_ids <- as.character(
      metadata$sample
    )
    
    
    count_not_metadata <- setdiff(
      count_sample_ids,
      metadata_sample_ids
    )
    
    
    metadata_not_count <- setdiff(
      metadata_sample_ids,
      count_sample_ids
    )
    
    
    cat(
      "\nCount samples missing from metadata:",
      length(count_not_metadata),
      "\n"
    )
    
    if (length(count_not_metadata) > 0) {
      
      print(
        head(count_not_metadata, 20)
      )
    }
    
    
    cat(
      "Metadata samples missing from count matrix:",
      length(metadata_not_count),
      "\n"
    )
    
    if (length(metadata_not_count) > 0) {
      
      print(
        head(metadata_not_count, 20)
      )
    }
    
    
    cat(
      "Matching sample IDs:",
      length(
        intersect(
          count_sample_ids,
          metadata_sample_ids
        )
      ),
      "\n"
    )
  }
  
  
  # ==========================================================
  # 9. METADATA INSPECTION
  # ==========================================================
  
  cat("\nMetadata preview:\n")
  
  print(head(metadata))
  
  
  # ----------------------------------------------------------
  # 9.1 Condition distribution
  # ----------------------------------------------------------
  
  if ("condition" %in% colnames(metadata)) {
    
    cat("\nCondition distribution:\n")
    
    print(
      table(
        metadata$condition,
        useNA = "ifany"
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 9.2 Tissue distribution
  # ----------------------------------------------------------
  
  if ("tissue" %in% colnames(metadata)) {
    
    cat("\nTissue distribution:\n")
    
    print(
      table(
        metadata$tissue,
        useNA = "ifany"
      )
    )
  }
  
  
  # ==========================================================
  # 10. CREATE DIAGNOSTIC SUMMARY
  # ==========================================================
  
  diagnostic_summary <- data.frame(
    
    GEO_ID = geo_id,
    
    count_rows = nrow(counts),
    
    count_columns = ncol(counts),
    
    metadata_rows = nrow(metadata),
    
    metadata_columns = ncol(metadata),
    
    gene_id_column = gene_id_column,
    
    unique_gene_ids = length(
      unique(gene_ids)
    ),
    
    duplicated_gene_ids = length(
      duplicated_gene_ids
    ),
    
    ensembl_no_version = sum(
      ensembl_no_version
    ),
    
    ensembl_with_version = sum(
      ensembl_with_version
    ),
    
    gene_symbol_column_present =
      symbol_available,
    
    rows_with_gene_symbols =
      ifelse(
        symbol_available,
        sum(!missing_symbols),
        NA
      ),
    
    unique_gene_symbols =
      ifelse(
        symbol_available,
        length(unique(symbols_clean)),
        NA
      ),
    
    duplicated_gene_symbols =
      ifelse(
        symbol_available,
        length(duplicated_symbols),
        NA
      ),
    
    missing_gene_symbols =
      ifelse(
        symbol_available,
        sum(missing_symbols),
        NA
      ),
    
    count_samples =
      length(count_sample_ids),
    
    metadata_samples =
      ifelse(
        "sample" %in% colnames(metadata),
        length(unique(metadata$sample)),
        NA
      ),
    
    matching_sample_ids =
      ifelse(
        "sample" %in% colnames(metadata),
        length(
          intersect(
            count_sample_ids,
            metadata$sample
          )
        ),
        NA
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ==========================================================
  # 11. SAVE SUMMARY
  # ==========================================================
  
  summary_file <- file.path(
    output_dir,
    paste0(
      geo_id,
      "_diagnostic_summary.csv"
    )
  )
  
  
  export(
    diagnostic_summary,
    summary_file
  )
  
  
  cat(
    "\nDiagnostic summary saved to:\n",
    summary_file,
    "\n"
  )
  
  
  cat(
    "\nFinished:",
    geo_id,
    "\n"
  )
}


# ============================================================
# END OF DIAGNOSTIC SCRIPT
# ============================================================