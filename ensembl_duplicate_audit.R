# ============================================================
# Ensembl ID / Gene Symbol Duplicate Audit
# Parkinson's disease RNA-seq datasets
#
# Purpose:
# Determine what multiple Ensembl IDs sharing the same
# gene symbol actually represent.
#
# For each Ensembl ID, retrieve:
#   - Ensembl gene ID
#   - Gene symbol
#   - Gene biotype
#   - Chromosome
#   - Start position
#   - End position
#   - Strand
#
# Then identify gene symbols associated with multiple
# Ensembl records.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(rio)
library(biomaRt)


# ------------------------------------------------------------
# 2. Define input/output directories
# ------------------------------------------------------------

raw_dir <- "/media/nabil/ssd/Parkinson-s/unedited_data/raw_counts"

output_dir <- "/media/nabil/ssd/Parkinson-s/results/tables/annotation_audit"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Define datasets
# ------------------------------------------------------------

geo_ids <- c(
  "GSE133101",
  "GSE136666",
  "GSE168496"
)


# ------------------------------------------------------------
# 4. Connect to Ensembl
# ------------------------------------------------------------

message("Connecting to Ensembl...")

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

message("Connected to Ensembl.")


# ------------------------------------------------------------
# 5. Retrieve annotation for one dataset
# ------------------------------------------------------------

for (geo_id in geo_ids) {
  
  message("\n============================")
  message("Processing ", geo_id)
  message("============================")
  
  
  # ----------------------------------------------------------
  # 5.1 Read raw count matrix
  # ----------------------------------------------------------
  
  input_file <- file.path(
    raw_dir,
    paste0(geo_id, "_raw_counts.csv")
  )
  
  if (!file.exists(input_file)) {
    
    message(
      "File not found: ",
      input_file
    )
    
    next
  }
  
  counts <- import(input_file)
  
  
  # ----------------------------------------------------------
  # 5.2 Identify Ensembl ID column
  # ----------------------------------------------------------
  
  if (!"GeneID" %in% colnames(counts)) {
    
    stop(
      "GeneID column not found in ",
      geo_id
    )
  }
  
  gene_ids <- counts$GeneID
  
  
  # ----------------------------------------------------------
  # 5.3 Remove Ensembl version numbers if present
  #
  # Example:
  # ENSG00000123456.7
  # becomes
  # ENSG00000123456
  # ----------------------------------------------------------
  
  gene_ids <- sub(
    "\\..*$",
    "",
    gene_ids
  )
  
  gene_ids <- unique(gene_ids)
  
  message(
    "Unique Ensembl IDs: ",
    length(gene_ids)
  )
  
  
  # ----------------------------------------------------------
  # 5.4 Query Ensembl
  # ----------------------------------------------------------
  
  annotation <- getBM(
    
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "gene_biotype",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    ),
    
    filters = "ensembl_gene_id",
    
    values = gene_ids,
    
    mart = ensembl
  )
  
  
  # ----------------------------------------------------------
  # 5.5 Rename columns
  # ----------------------------------------------------------
  
  annotation <- annotation %>%
    
    rename(
      
      GeneID = ensembl_gene_id,
      
      Gene_Symbol = external_gene_name,
      
      Gene_Biotype = gene_biotype,
      
      Chromosome = chromosome_name,
      
      Start = start_position,
      
      End = end_position,
      
      Strand = strand
      
    )
  
  
  # ----------------------------------------------------------
  # 5.6 Add annotation status
  # ----------------------------------------------------------
  
  annotation <- annotation %>%
    
    mutate(
      
      Annotation_Found =
        GeneID %in% annotation$GeneID
      
    )
  
  
  # ----------------------------------------------------------
  # 5.7 Determine how many Ensembl records share
  #     each gene symbol
  # ----------------------------------------------------------
  
  annotation <- annotation %>%
    mutate(
      Gene_Symbol = na_if(Gene_Symbol, ""),
      Number_of_Ensembl_IDs_for_Symbol = NA_integer_,
      Symbol_Duplicated = FALSE
    ) %>%
    group_by(Gene_Symbol) %>%
    mutate(
      Number_of_Ensembl_IDs_for_Symbol =
        if_else(
          !is.na(Gene_Symbol),
          n(),
          NA_integer_
        ),
      Symbol_Duplicated =
        !is.na(Gene_Symbol) & n() > 1
    ) %>%
    ungroup()
  
  
  # ----------------------------------------------------------
  # 5.8 Order results
  # ----------------------------------------------------------
  
  annotation <- annotation %>%
    
    arrange(
      
      desc(Symbol_Duplicated),
      
      Gene_Symbol,
      
      GeneID
      
    )
  
  
  # ----------------------------------------------------------
  # 5.9 Save complete annotation
  # ----------------------------------------------------------
  
  complete_file <- file.path(
    
    output_dir,
    
    paste0(
      geo_id,
      "_Ensembl_annotation.csv"
    )
    
  )
  
  export(
    annotation,
    complete_file
  )
  
  
  # ----------------------------------------------------------
  # 5.10 Extract duplicated gene symbols only
  # ----------------------------------------------------------
  
  duplicated_symbols <- annotation %>%
    
    filter(
      Symbol_Duplicated == TRUE,
      !is.na(Gene_Symbol),
      Gene_Symbol != ""
    )
  
  
  # ----------------------------------------------------------
  # 5.11 Save duplicate report
  # ----------------------------------------------------------
  
  duplicate_file <- file.path(
    
    output_dir,
    
    paste0(
      geo_id,
      "_duplicated_symbols.csv"
    )
    
  )
  
  export(
    duplicated_symbols,
    duplicate_file
  )
  
  
  # ----------------------------------------------------------
  # 5.12 Print summary
  # ----------------------------------------------------------
  
  message(
    "Total Ensembl IDs queried: ",
    length(gene_ids)
  )
  
  message(
    "Ensembl IDs returned by Ensembl: ",
    nrow(annotation)
  )
  
  message(
    "Unique gene symbols: ",
    n_distinct(
      annotation$Gene_Symbol,
      na.rm = TRUE
    )
  )
  
  message(
    "Duplicated gene symbols: ",
    n_distinct(
      duplicated_symbols$Gene_Symbol
    )
  )
  
  message(
    "Saved: ",
    complete_file
  )
  
  message(
    "Saved: ",
    duplicate_file
  )
  
}


# ============================================================
# END
# ============================================================