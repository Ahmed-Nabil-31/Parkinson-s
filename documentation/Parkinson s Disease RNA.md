# ============================================================
# Parkinson's Disease RNA-seq
# Raw Count Matrix + Metadata Inspection
# ============================================================
#
# Purpose:
#   1. Load the raw count matrices for all three datasets
#   2. Load their corresponding metadata
#   3. Inspect gene ID structure
#   4. Check sample IDs between counts and metadata
#   5. Identify duplicated Ensembl IDs
#   6. Identify the number of unique gene IDs
#   7. Inspect available metadata columns
#
# IMPORTANT:
#   We are NOT collapsing gene symbols yet.
#   We first verify how the original count matrices are annotated.
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
# 2. Define input directories
# ------------------------------------------------------------

raw_dir <- "/media/nabil/ssd/Parkinson-s/unedited_data/raw_counts"

metadata_dir <- "/media/nabil/ssd/Parkinson-s/unedited_data/metadata"


# ------------------------------------------------------------
# 3. Define GEO dataset IDs
# ------------------------------------------------------------

geo_ids <- c(
  "GSE133101",
  "GSE136666",
  "GSE168496"
)


# ------------------------------------------------------------
# 4. Create diagnostic output directory
# ------------------------------------------------------------

dir.create(
  "results/tables/preprocessing",
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 5. Process each dataset
# ============================================================

for (geo_id in geo_ids) {
  
  message("\n==========================================")
  message("Processing: ", geo_id)
  message("==========================================\n")
  
  
  # ----------------------------------------------------------
  # 5.1 Define file paths
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
  # 5.2 Check that files exist
  # ----------------------------------------------------------
  
  if (!file.exists(count_file)) {
    
    message(
      "WARNING: Raw count file not found: ",
      count_file
    )
    
    next
  }
  
  if (!file.exists(metadata_file)) {
    
    message(
      "WARNING: Metadata file not found: ",
      metadata_file
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # 5.3 Load raw counts and metadata
  # ----------------------------------------------------------
  
  counts <- import(count_file)
  
  metadata <- import(metadata_file)
  
  
  # ----------------------------------------------------------
  # 5.4 Basic dimensions
  # ----------------------------------------------------------
  
  message("Raw count dimensions:")
  message(
    "Rows (features): ",
    nrow(counts),
    " | Columns: ",
    ncol(counts)
  )
  
  message("\nMetadata dimensions:")
  message(
    "Rows: ",
    nrow(metadata),
    " | Columns: ",
    ncol(metadata)
  )
  
  
  # ----------------------------------------------------------
  # 5.5 Inspect column names
  # ----------------------------------------------------------
  
  message("\nRaw count column names:")
  print(head(colnames(counts), 15))
  
  
  message("\nMetadata column names:")
  print(colnames(metadata))
  
  
  # ----------------------------------------------------------
  # 5.6 Identify the gene ID column
  #
  # Usually this will be the first column.
  # ----------------------------------------------------------
  
  gene_id_column <- colnames(counts)[1]
  
  message(
    "\nAssuming first count-matrix column is gene ID: ",
    gene_id_column
  )
  
  
  gene_ids <- as.character(counts[[gene_id_column]])
  
  
  # ----------------------------------------------------------
  # 5.7 Inspect first gene IDs
  # ----------------------------------------------------------
  
  message("\nFirst 20 gene IDs:")
  print(head(gene_ids, 20))
  
  
  # ----------------------------------------------------------
  # 5.8 Check whether IDs look like Ensembl gene IDs
  # ----------------------------------------------------------
  
  ensembl_no_version <- grepl(
    "^ENSG[0-9]+$",
    gene_ids
  )
  
  ensembl_with_version <- grepl(
    "^ENSG[0-9]+\\.[0-9]+$",
    gene_ids
  )
  
  
  message(
    "\nNumber of IDs matching ENSG format: ",
    sum(ensembl_no_version)
  )
  
  message(
    "Number of IDs matching ENSG.version format: ",
    sum(ensembl_with_version)
  )
  
  
  # ----------------------------------------------------------
  # 5.9 Check duplicated gene IDs
  # ----------------------------------------------------------
  
  duplicated_gene_ids <- unique(
    gene_ids[duplicated(gene_ids)]
  )
  
  
  message(
    "\nNumber of duplicated gene IDs: ",
    length(duplicated_gene_ids)
  )
  
  
  if (length(duplicated_gene_ids) > 0) {
    
    message("First duplicated gene IDs:")
    
    print(
      head(duplicated_gene_ids, 20)
    )
  }
  
  
  # ----------------------------------------------------------
  # 5.10 Number of unique gene IDs
  # ----------------------------------------------------------
  
  message(
    "\nUnique gene IDs: ",
    length(unique(gene_ids))
  )
  
  
  # ----------------------------------------------------------
  # 5.11 Check for gene-symbol columns
  # ----------------------------------------------------------
  
  possible_symbol_columns <- colnames(counts)[
    grepl(
      "symbol|gene.?name|gene.?symbol",
      colnames(counts),
      ignore.case = TRUE
    )
  ]
  
  
  message("\nPossible gene-symbol columns:")
  
  if (length(possible_symbol_columns) == 0) {
    
    message("None found.")
    
  } else {
    
    print(possible_symbol_columns)
  }
  
  
  # ----------------------------------------------------------
  # 5.12 Check whether sample IDs in metadata
  #      occur in the count matrix
  # ----------------------------------------------------------
  
  count_sample_ids <- colnames(counts)[-1]
  
  
  message("\nNumber of count-matrix samples:")
  message(length(count_sample_ids))
  
  
  message("\nFirst metadata rows:")
  print(head(metadata))
  
  
  # ----------------------------------------------------------
  # Try to identify a likely sample-ID column
  # ----------------------------------------------------------
  
  possible_sample_columns <- colnames(metadata)[
    grepl(
      "sample|gsm|accession",
      colnames(metadata),
      ignore.case = TRUE
    )
  ]
  
  
  message("\nPossible metadata sample-ID columns:")
  
  if (length(possible_sample_columns) == 0) {
    
    message("None automatically identified.")
    
  } else {
    
    print(possible_sample_columns)
  }
  
  
  # ----------------------------------------------------------
  # 5.13 Save diagnostic summary
  # ----------------------------------------------------------
  
  diagnostic_summary <- data.frame(
    
    GEO_ID = geo_id,
    
    count_rows = nrow(counts),
    
    count_columns = ncol(counts),
    
    metadata_rows = nrow(metadata),
    
    metadata_columns = ncol(metadata),
    
    gene_id_column = gene_id_column,
    
    unique_gene_ids = length(unique(gene_ids)),
    
    duplicated_gene_ids = length(
      duplicated_gene_ids
    ),
    
    ensembl_no_version = sum(
      ensembl_no_version
    ),
    
    ensembl_with_version = sum(
      ensembl_with_version
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  export(
    diagnostic_summary,
    file.path(
      "results/tables/preprocessing",
      paste0(
        geo_id,
        "_diagnostic_summary.csv"
      )
    )
  )
  
  
  message(
    "\nFinished diagnostic inspection for ",
    geo_id
  )
}

# GSE133101
setdiff(
  colnames(import("unedited_data/raw_counts/GSE133101_raw_counts.csv"))[-c(1,2)],
  import("unedited_data/metadata/GSE133101_metadata.csv")$sample
)

# GSE136666
setdiff(
  colnames(import("unedited_data/raw_counts/GSE136666_raw_counts.csv"))[-c(1,2)],
  import("unedited_data/metadata/GSE136666_metadata.csv")$sample
)

# GSE168496
setdiff(
  colnames(import("unedited_data/raw_counts/GSE168496_raw_counts.csv"))[-c(1,2)],
  import("unedited_data/metadata/GSE168496_metadata.csv")$sample
)

library(rio)

raw_dir <- "/media/nabil/ssd/Parkinson-s/unedited_data/raw_counts"

geo_ids <- c("GSE133101", "GSE136666", "GSE168496")

for (geo in geo_ids) {
  
  file <- file.path(raw_dir, paste0(geo, "_raw_counts.csv"))
  counts <- import(file)
  
  symbols <- counts$gene_symbol
  
  # Remove blank/NA symbols for the duplicate check
  symbols_clean <- symbols[!is.na(symbols) & symbols != ""]
  
  dup_symbols <- unique(symbols_clean[duplicated(symbols_clean)])
  
  cat("\n============================\n")
  cat(geo, "\n")
  cat("============================\n")
  
  cat("Total rows:", nrow(counts), "\n")
  cat("Rows with gene symbols:", length(symbols_clean), "\n")
  cat("Unique gene symbols:", length(unique(symbols_clean)), "\n")
  cat("Duplicated gene symbols:", length(dup_symbols), "\n")
  cat("Blank/NA gene symbols:",
      sum(is.na(symbols) | symbols == ""), "\n")
  
  if (length(dup_symbols) > 0) {
    cat("\nFirst duplicated symbols:\n")
    print(head(dup_symbols, 20))
  }
}

# ============================================================
# END OF DIAGNOSTIC SCRIPT
# ============================================================
