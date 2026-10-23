# ============================================================
# FINAL TISSUE FILTER + STRICT HAPLOTYPE COLLAPSE
# ============================================================
#
# Datasets:
#   GSE133101
#   GSE136666
#   GSE168496
#
# Tissue:
#   Substantia nigra
#
# STRICT HAPLOTYPE MERGE:
#
#   1 primary chromosome record
#   + >=1 alternative haplotype record
#   + same Gene_Symbol
#   + same Gene_Biotype
#   + same Chromosome_Family
#   + no OTHER chromosome records
#
# Primary GeneID is retained.
# Alternative haplotype counts are summed into the
# corresponding primary GeneID.
#
# IMPORTANT:
# Ambiguous groups are detected automatically.
# No gene symbols are hardcoded.
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(tidyverse)


# ============================================================
# 2. PATHS
# ============================================================

base_dir <- "/media/nabil/ssd/Parkinson-s"


raw_count_dir <- file.path(
  base_dir,
  "unedited_data",
  "raw_counts"
)


metadata_dir <- file.path(
  base_dir,
  "unedited_data",
  "metadata"
)


annotation_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "annotation_audit"
)


output_dir <- file.path(
  base_dir,
  "results",
  "tables",
  "final_tissue_haplotype_collapsed"
)


audit_dir <- file.path(
  output_dir,
  "audit"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  audit_dir,
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


# ============================================================
# 4. TARGET TISSUE
# ============================================================

TARGET_TISSUE <- "Substantia nigra"


# ============================================================
# 5. HELPER FUNCTION
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
      )
    )
  }
  
  found[1]
}


# ============================================================
# 6. PROCESS ONE DATASET
# ============================================================

process_dataset <- function(gse) {
  
  cat("\n\n")
  cat("============================================================\n")
  cat("PROCESSING ", gse, "\n", sep = "")
  cat("============================================================\n")
  
  
  # ==========================================================
  # 6.1 FILE PATHS
  # ==========================================================
  
  count_file <- file.path(
    raw_count_dir,
    paste0(
      gse,
      "_raw_counts.csv"
    )
  )
  
  
  metadata_file <- file.path(
    metadata_dir,
    paste0(
      gse,
      "_metadata.csv"
    )
  )
  
  
  annotation_file <- file.path(
    annotation_dir,
    paste0(
      gse,
      "_Ensembl_annotation.csv"
    )
  )
  
  
  # ==========================================================
  # 6.2 CHECK FILES
  # ==========================================================
  
  required_files <- c(
    count_file,
    metadata_file,
    annotation_file
  )
  
  
  missing_files <- required_files[
    !file.exists(required_files)
  ]
  
  
  if (length(missing_files) > 0) {
    
    stop(
      "\nMissing required files for ",
      gse,
      ":\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  }
  
  
  # ==========================================================
  # 6.3 READ FILES
  # ==========================================================
  
  counts <- readr::read_csv(
    count_file,
    show_col_types = FALSE
  )
  
  
  metadata <- readr::read_csv(
    metadata_file,
    show_col_types = FALSE
  )
  
  
  annotation <- readr::read_csv(
    annotation_file,
    show_col_types = FALSE
  )
  
  
  cat(
    "Original count dimensions: ",
    nrow(counts),
    " rows x ",
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
  # 6.4 CHECK REQUIRED ANNOTATION COLUMNS
  # ==========================================================
  
  required_annotation_columns <- c(
    "GeneID",
    "Gene_Symbol",
    "Gene_Biotype",
    "Chromosome"
  )
  
  
  missing_annotation_columns <- setdiff(
    required_annotation_columns,
    colnames(annotation)
  )
  
  
  if (length(missing_annotation_columns) > 0) {
    
    stop(
      "\nMissing annotation columns in ",
      gse,
      ":\n",
      paste(
        missing_annotation_columns,
        collapse = ", "
      ),
      "\n\nAvailable annotation columns:\n",
      paste(
        colnames(annotation),
        collapse = ", "
      )
    )
  }
  
  
  # ==========================================================
  # 6.5 STANDARDISE ANNOTATION
  # ==========================================================
  
  annotation <- annotation |>
    
    dplyr::mutate(
      
      GeneID =
        as.character(GeneID),
      
      Gene_Symbol =
        as.character(Gene_Symbol),
      
      Gene_Biotype =
        as.character(Gene_Biotype),
      
      Chromosome =
        as.character(Chromosome)
      
    )
  
  
  # ==========================================================
  # 6.6 CHECK ANNOTATION GENEID UNIQUENESS
  # ==========================================================
  
  duplicated_annotation <- annotation |>
    
    dplyr::count(
      GeneID,
      name = "n"
    ) |>
    
    dplyr::filter(
      n > 1
    )
  
  
  if (nrow(duplicated_annotation) > 0) {
    
    readr::write_csv(
      duplicated_annotation,
      file.path(
        audit_dir,
        paste0(
          gse,
          "_DUPLICATED_ANNOTATION_GENEIDS.csv"
        )
      )
    )
    
    
    stop(
      "\nAnnotation contains duplicated GeneIDs in ",
      gse
    )
  }
  
  
  # ==========================================================
  # 6.7 DERIVE CHROMOSOME TYPE
  # ==========================================================
  
  annotation <- annotation |>
    
    dplyr::mutate(
      
      Chromosome_Type = dplyr::case_when(
        
        Chromosome %in%
          c(
            as.character(1:22),
            "X",
            "Y",
            "MT"
          ) ~
          "PRIMARY_CHROMOSOME",
        
        stringr::str_detect(
          Chromosome,
          "^HSCHR"
        ) ~
          "ALTERNATIVE_HAPLOTYPE",
        
        TRUE ~
          "OTHER"
      )
    )
  
  
  # ==========================================================
  # 6.8 DERIVE CHROMOSOME FAMILY
  #
  # IMPORTANT:
  # Keep this as CHARACTER throughout.
  # This prevents coercion warnings.
  # ==========================================================
  
  annotation <- annotation |>
    
    dplyr::mutate(
      
      Chromosome_Family = dplyr::case_when(
        
        Chromosome %in%
          as.character(1:22) ~
          Chromosome,
        
        Chromosome == "X" ~
          "X",
        
        Chromosome == "Y" ~
          "Y",
        
        Chromosome == "MT" ~
          "MT",
        
        stringr::str_detect(
          Chromosome,
          "^HSCHRX"
        ) ~
          "X",
        
        stringr::str_detect(
          Chromosome,
          "^HSCHRY"
        ) ~
          "Y",
        
        stringr::str_detect(
          Chromosome,
          "^HSCHRMT"
        ) ~
          "MT",
        
        stringr::str_detect(
          Chromosome,
          "^HSCHR[0-9]+"
        ) ~
          stringr::str_extract(
            Chromosome,
            "(?<=HSCHR)[0-9]+"
          ),
        
        TRUE ~
          NA_character_
      )
    )
  
  
  # ==========================================================
  # 6.9 CHECK ALTERNATIVE HAPLOTYPES WITHOUT FAMILY
  # ==========================================================
  
  haplotype_family_problems <- annotation |>
    
    dplyr::filter(
      
      Chromosome_Type ==
        "ALTERNATIVE_HAPLOTYPE",
      
      is.na(Chromosome_Family)
    ) |>
    
    dplyr::select(
      GeneID,
      Gene_Symbol,
      Chromosome
    )
  
  
  if (nrow(haplotype_family_problems) > 0) {
    
    readr::write_csv(
      haplotype_family_problems,
      file.path(
        audit_dir,
        paste0(
          gse,
          "_HAPLOTYPE_FAMILY_PROBLEMS.csv"
        )
      )
    )
    
    
    stop(
      "\nAlternative haplotypes without chromosome family found in ",
      gse,
      ".\nCheck:\n",
      file.path(
        audit_dir,
        paste0(
          gse,
          "_HAPLOTYPE_FAMILY_PROBLEMS.csv"
        )
      )
    )
  }
  
  
  # ==========================================================
  # 6.10 STANDARDISE COUNTS
  # ==========================================================
  
  if (
    !"GeneID" %in%
    colnames(counts)
  ) {
    
    stop(
      "\nGeneID column missing from count matrix for ",
      gse
    )
  }
  
  
  counts <- counts |>
    
    dplyr::mutate(
      GeneID =
        as.character(GeneID)
    )
  
  
  # ==========================================================
  # 6.11 IDENTIFY GENE SYMBOL COLUMN
  # ==========================================================
  
  if (
    "gene_symbol" %in%
    colnames(counts)
  ) {
    
    count_symbol_column <- "gene_symbol"
    
  } else if (
    "Gene_Symbol" %in%
    colnames(counts)
  ) {
    
    counts <- counts |>
      
      dplyr::rename(
        gene_symbol = Gene_Symbol
      )
    
    count_symbol_column <- "gene_symbol"
    
  } else {
    
    stop(
      "\nNo gene symbol column found in counts for ",
      gse
    )
  }
  
  
  counts <- counts |>
    
    dplyr::mutate(
      gene_symbol =
        as.character(gene_symbol)
    )
  
  
  # ==========================================================
  # 6.12 RAW COUNT GENEID DUPLICATE CHECK
  # ==========================================================
  
  duplicated_raw <- counts |>
    
    dplyr::count(
      GeneID,
      name = "n"
    ) |>
    
    dplyr::filter(
      n > 1
    )
  
  
  if (nrow(duplicated_raw) > 0) {
    
    readr::write_csv(
      duplicated_raw,
      file.path(
        audit_dir,
        paste0(
          gse,
          "_DUPLICATED_RAW_COUNT_GENEIDS.csv"
        )
      )
    )
    
    
    stop(
      "\nRaw count matrix already contains duplicated GeneIDs in ",
      gse
    )
  }
  
  
  # ==========================================================
  # 6.13 FIND METADATA COLUMNS
  # ==========================================================
  
  sample_col <- find_column(
    metadata,
    c(
      "sample",
      "Sample",
      "SampleID",
      "sample_id",
      "GSM"
    ),
    gse,
    "sample"
  )
  
  
  tissue_col <- find_column(
    metadata,
    c(
      "tissue",
      "Tissue",
      "Tissue_Type",
      "tissue_type"
    ),
    gse,
    "tissue"
  )
  
  
  # ==========================================================
  # 6.14 STANDARDISE METADATA
  # ==========================================================
  
  metadata <- metadata |>
    
    dplyr::mutate(
      
      sample_id =
        as.character(
          .data[[sample_col]]
        ),
      
      tissue_value =
        as.character(
          .data[[tissue_col]]
        )
    )
  
  
  # ==========================================================
  # 6.15 DISPLAY AVAILABLE TISSUES
  # ==========================================================
  
  cat("\nAvailable tissues:\n")
  
  
  print(
    table(
      metadata$tissue_value,
      useNA = "ifany"
    )
  )
  
  
  # ==========================================================
  # 6.16 FILTER TO SUBSTANTIA NIGRA
  # ==========================================================
  
  metadata_filtered <- metadata |>
    
    dplyr::filter(
      
      stringr::str_to_lower(
        stringr::str_trim(
          tissue_value
        )
      ) ==
        stringr::str_to_lower(
          TARGET_TISSUE
        )
    )
  
  
  if (
    nrow(metadata_filtered) == 0
  ) {
    
    stop(
      "\nNo ",
      TARGET_TISSUE,
      " samples found in ",
      gse
    )
  }
  
  
  cat(
    "\nTissue selected for analysis: ",
    TARGET_TISSUE,
    "\n",
    sep = ""
  )
  
  
  cat(
    "Samples retained after tissue filtering: ",
    nrow(metadata_filtered),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.17 MATCH SAMPLES
  # ==========================================================
  
  count_sample_columns <- setdiff(
    colnames(counts),
    c(
      "GeneID",
      "gene_symbol"
    )
  )
  
  
  selected_samples <- intersect(
    metadata_filtered$sample_id,
    count_sample_columns
  )
  
  
  if (
    length(selected_samples) == 0
  ) {
    
    stop(
      "\nNo tissue-filtered samples matched counts for ",
      gse
    )
  }
  
  
  cat(
    "Matched count samples: ",
    length(selected_samples),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.18 TISSUE FILTER COUNT MATRIX
  # ==========================================================
  
  counts_filtered <- counts |>
    
    dplyr::select(
      GeneID,
      gene_symbol,
      dplyr::all_of(selected_samples)
    )
  
  
  cat(
    "Count matrix after tissue filtering: ",
    nrow(counts_filtered),
    " genes x ",
    length(selected_samples),
    " samples\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.19 CREATE ANNOTATION GROUPS
  # ==========================================================
  
  annotation_groups <- annotation |>
    
    dplyr::filter(
      
      !is.na(Gene_Symbol),
      
      Gene_Symbol != "",
      
      !is.na(Chromosome_Family)
    ) |>
    
    dplyr::group_by(
      
      Gene_Symbol,
      
      Gene_Biotype,
      
      Chromosome_Family
    ) |>
    
    dplyr::summarise(
      
      n_records = dplyr::n(),
      
      n_primary =
        sum(
          Chromosome_Type ==
            "PRIMARY_CHROMOSOME"
        ),
      
      n_haplotype =
        sum(
          Chromosome_Type ==
            "ALTERNATIVE_HAPLOTYPE"
        ),
      
      n_other =
        sum(
          Chromosome_Type ==
            "OTHER"
        ),
      
      chromosomes =
        paste(
          unique(Chromosome),
          collapse = "; "
        ),
      
      gene_ids =
        paste(
          unique(GeneID),
          collapse = "; "
        ),
      
      .groups = "drop"
    )
  
  
  # ==========================================================
  # 6.20 IDENTIFY AMBIGUOUS GROUPS
  #
  # Ambiguous:
  # More than one primary GeneID exists within
  # the same symbol/biotype/chromosome-family group.
  #
  # These groups are NOT collapsed.
  # ==========================================================
  
  ambiguous_groups <- annotation_groups |>
    
    dplyr::filter(
      
      n_primary > 1
    ) |>
    
    dplyr::mutate(
      
      Merge_Status =
        "MULTIPLE_PRIMARY_IDS"
    )
  
  
  # Save ambiguous groups even if empty.
  
  readr::write_csv(
    
    ambiguous_groups,
    
    file.path(
      audit_dir,
      paste0(
        gse,
        "_AMBIGUOUS_GROUPS.csv"
      )
    )
  )
  
  
  cat(
    "Ambiguous groups excluded: ",
    nrow(ambiguous_groups),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.21 CREATE VALID STRICT MERGE GROUPS
  # ==========================================================
  
  valid_group_keys <- annotation_groups |>
    
    dplyr::filter(
      
      n_primary == 1,
      
      n_haplotype >= 1,
      
      n_other == 0,
      
      !(
        Gene_Symbol %in%
          ambiguous_groups$Gene_Symbol
      )
    ) |>
    
    dplyr::select(
      
      Gene_Symbol,
      
      Gene_Biotype,
      
      Chromosome_Family
    )
  
  
  cat(
    "Number of valid haplotype merge groups: ",
    nrow(valid_group_keys),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.22 BUILD MERGE MAP
  # ==========================================================
  
  merge_map <- annotation |>
    
    dplyr::inner_join(
      
      valid_group_keys,
      
      by = c(
        "Gene_Symbol",
        "Gene_Biotype",
        "Chromosome_Family"
      )
    ) |>
    
    dplyr::group_by(
      
      Gene_Symbol,
      
      Gene_Biotype,
      
      Chromosome_Family
    ) |>
    
    dplyr::mutate(
      
      Representative_GeneID =
        GeneID[
          Chromosome_Type ==
            "PRIMARY_CHROMOSOME"
        ][1],
      
      Action = dplyr::case_when(
        
        Chromosome_Type ==
          "PRIMARY_CHROMOSOME" ~
          "REPRESENTATIVE_PRIMARY",
        
        Chromosome_Type ==
          "ALTERNATIVE_HAPLOTYPE" ~
          "MERGE_INTO_PRIMARY",
        
        TRUE ~
          "DO_NOT_MERGE"
      )
    ) |>
    
    dplyr::ungroup() |>
    
    dplyr::select(
      
      GeneID,
      
      Gene_Symbol,
      
      Gene_Biotype,
      
      Chromosome,
      
      Chromosome_Type,
      
      Chromosome_Family,
      
      Representative_GeneID,
      
      Action
    )
  
  
  # ==========================================================
  # 6.23 VERIFY MERGE MAP
  # ==========================================================
  
  duplicate_mapping <- merge_map |>
    
    dplyr::count(
      GeneID,
      name = "n"
    ) |>
    
    dplyr::filter(
      n > 1
    )
  
  
  if (
    nrow(duplicate_mapping) > 0
  ) {
    
    readr::write_csv(
      
      duplicate_mapping,
      
      file.path(
        audit_dir,
        paste0(
          gse,
          "_DUPLICATED_MERGE_MAP_GENEIDS.csv"
        )
      )
    )
    
    
    stop(
      "\nERROR: one GeneID maps to multiple merge groups in ",
      gse
    )
  }
  
  
  # ==========================================================
  # 6.24 SAVE MERGE MAP
  # ==========================================================
  
  readr::write_csv(
    
    merge_map,
    
    file.path(
      audit_dir,
      paste0(
        gse,
        "_FINAL_MERGE_MAP.csv"
      )
    )
  )
  
  
  # ==========================================================
  # 6.25 CREATE HAPLOTYPE LOOKUP
  # ==========================================================
  
  collapse_ids <- merge_map |>
    
    dplyr::filter(
      
      Action ==
        "MERGE_INTO_PRIMARY"
    ) |>
    
    dplyr::select(
      
      GeneID,
      
      Representative_GeneID
    )
  
  
  cat(
    "Alternative haplotype rows to collapse: ",
    nrow(collapse_ids),
    "\n",
    sep = ""
  )
  
  
  # ==========================================================
  # 6.26 SAFE LOOKUP VECTOR
  # ==========================================================
  
  collapse_lookup <- stats::setNames(
    
    collapse_ids$Representative_GeneID,
    
    collapse_ids$GeneID
  )
  
  
  # ==========================================================
  # 6.27 ASSIGN FINAL GENE IDs
  # ==========================================================
  
  counts_for_collapse <- counts_filtered |>
    
    dplyr::mutate(
      
      Final_GeneID =
        dplyr::if_else(
          
          GeneID %in%
            names(collapse_lookup),
          
          unname(
            collapse_lookup[GeneID]
          ),
          
          GeneID
        )
    )
  
  
  # ==========================================================
  # 6.28 COLLAPSE COUNTS
  # ==========================================================
  
  final_counts <- counts_for_collapse |>
    
    dplyr::group_by(
      Final_GeneID
    ) |>
    
    dplyr::summarise(
      
      dplyr::across(
        
        dplyr::all_of(
          selected_samples
        ),
        
        ~ sum(
          as.numeric(.x),
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
    ) |>
    
    dplyr::rename(
      GeneID = Final_GeneID
    )
  
  
  # ==========================================================
  # 6.29 CREATE SYMBOL LOOKUP
  # ==========================================================
  
  symbol_lookup <- counts_filtered |>
    
    dplyr::select(
      GeneID,
      gene_symbol
    ) |>
    
    dplyr::distinct(
      GeneID,
      .keep_all = TRUE
    )
  
  
  # ==========================================================
  # 6.30 ADD SYMBOLS
  # ==========================================================
  
  final_counts <- final_counts |>
    
    dplyr::left_join(
      symbol_lookup,
      by = "GeneID"
    ) |>
    
    dplyr::select(
      GeneID,
      gene_symbol,
      dplyr::everything()
    )
  
  
  # ==========================================================
  # 6.31 FINAL DUPLICATE GeneID CHECK
  # ==========================================================
  
  duplicated_final <- final_counts |>
    
    dplyr::count(
      GeneID,
      name = "n"
    ) |>
    
    dplyr::filter(
      n > 1
    )
  
  
  if (
    nrow(duplicated_final) > 0
  ) {
    
    readr::write_csv(
      
      duplicated_final,
      
      file.path(
        audit_dir,
        paste0(
          gse,
          "_FINAL_DUPLICATED_GENEIDS.csv"
        )
      )
    )
    
    
    stop(
      "\nERROR: duplicated GeneIDs remain after final collapse in ",
      gse,
      ".\nCheck:\n",
      file.path(
        audit_dir,
        paste0(
          gse,
          "_FINAL_DUPLICATED_GENEIDS.csv"
        )
      )
    )
  }
  
  
  # ==========================================================
  # 6.32 CHECK DUPLICATE SYMBOLS
  #
  # Duplicate symbols are NOT automatically an error.
  # ==========================================================
  
  duplicated_symbols <- final_counts |>
    
    dplyr::count(
      gene_symbol,
      name = "n"
    ) |>
    
    dplyr::filter(
      
      !is.na(gene_symbol),
      
      gene_symbol != "",
      
      n > 1
    )
  
  
  readr::write_csv(
    
    duplicated_symbols,
    
    file.path(
      audit_dir,
      paste0(
        gse,
        "_FINAL_DUPLICATED_SYMBOLS.csv"
      )
    )
  )
  
  
  # ==========================================================
  # 6.33 SAVE FINAL COUNTS
  # ==========================================================
  
  final_count_file <- file.path(
    
    output_dir,
    
    paste0(
      gse,
      "_Substantia_nigra_FINAL_counts.csv"
    )
  )
  
  
  readr::write_csv(
    final_counts,
    final_count_file
  )
  
  
  # ==========================================================
  # 6.34 SAVE FILTERED METADATA
  # ==========================================================
  
  metadata_output <- metadata_filtered |>
    
    dplyr::select(
      -tissue_value
    )
  
  
  metadata_file_out <- file.path(
    
    output_dir,
    
    paste0(
      gse,
      "_Substantia_nigra_metadata.csv"
    )
  )
  
  
  readr::write_csv(
    metadata_output,
    metadata_file_out
  )
  
  
  # ==========================================================
  # 6.35 SAVE PROCESSING SUMMARY
  # ==========================================================
  
  summary <- tibble::tibble(
    
    Dataset =
      gse,
    
    Tissue =
      TARGET_TISSUE,
    
    Original_Genes =
      nrow(counts),
    
    Tissue_Filtered_Genes =
      nrow(counts_filtered),
    
    Final_Genes =
      nrow(final_counts),
    
    Rows_Removed =
      nrow(counts_filtered) -
      nrow(final_counts),
    
    Valid_Merge_Groups =
      nrow(valid_group_keys),
    
    Haplotype_Rows_Collapsed =
      nrow(collapse_ids),
    
    Ambiguous_Groups_Excluded =
      nrow(ambiguous_groups),
    
    Final_Duplicated_GeneIDs =
      nrow(duplicated_final),
    
    Final_Duplicated_Symbols =
      nrow(duplicated_symbols)
  )
  
  
  readr::write_csv(
    
    summary,
    
    file.path(
      audit_dir,
      paste0(
        gse,
        "_FINAL_PROCESSING_SUMMARY.csv"
      )
    )
  )
  
  
  # ==========================================================
  # 6.36 REPORT
  # ==========================================================
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("FINAL RESULT: ", gse, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  
  
  cat(
    "Original genes: ",
    nrow(counts),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Substantia nigra genes: ",
    nrow(counts_filtered),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Final genes: ",
    nrow(final_counts),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Rows removed by haplotype collapse: ",
    nrow(counts_filtered) -
      nrow(final_counts),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Valid merge groups: ",
    nrow(valid_group_keys),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Haplotype rows collapsed: ",
    nrow(collapse_ids),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Ambiguous groups excluded: ",
    nrow(ambiguous_groups),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Final duplicated GeneIDs: ",
    nrow(duplicated_final),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Final duplicated symbols: ",
    nrow(duplicated_symbols),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Finished: ",
    gse,
    "\n",
    sep = ""
  )
  
  
  return(summary)
}


# ============================================================
# 7. RUN ALL DATASETS
# ============================================================

results <- purrr::map(
  datasets,
  process_dataset
)


# ============================================================
# 8. COMBINE SUMMARIES
# ============================================================

summary_files <- file.path(
  
  audit_dir,
  
  paste0(
    datasets,
    "_FINAL_PROCESSING_SUMMARY.csv"
  )
)


all_summaries <- purrr::map_dfr(
  
  summary_files,
  
  ~ readr::read_csv(
    .x,
    show_col_types = FALSE
  )
)


# ============================================================
# 9. SAVE MASTER SUMMARY
# ============================================================

readr::write_csv(
  
  all_summaries,
  
  file.path(
    audit_dir,
    "ALL_DATASETS_FINAL_PROCESSING_SUMMARY.csv"
  )
)


# ============================================================
# 10. DONE
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("ALL DATASETS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")


print(
  all_summaries
)