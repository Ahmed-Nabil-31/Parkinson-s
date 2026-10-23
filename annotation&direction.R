# ============================================================
# ANNOTATE 9 COMMON GENES + 4 GSE133101-SPECIFIC GENES
# Parkinson's disease RNA-seq analysis
#
# GSE136666 ∩ GSE168496 = common genes
# GSE133101 = genes unique to that dataset
#
# Uses EXISTING annotation files.
# No new Ensembl query.
# ============================================================


library(tidyverse)
library(rio)


# ------------------------------------------------------------
# 1. Directories
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

output_dir <- file.path(
  base_dir,
  "results/tables/final_recurrent_DEGs"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. DEG thresholds
# ------------------------------------------------------------

padj_cutoff <- 0.05
lfc_cutoff <- 1


# ------------------------------------------------------------
# 3. Read DESeq2 results
# ------------------------------------------------------------

gse133101 <- import(
  file.path(
    deseq_dir,
    "GSE133101.csv"
  )
)

gse136666 <- import(
  file.path(
    deseq_dir,
    "GSE136666.csv"
  )
)

gse168496 <- import(
  file.path(
    deseq_dir,
    "GSE168496.csv"
  )
)


# ------------------------------------------------------------
# 4. Standardize Gene_ID
# ------------------------------------------------------------

gse133101 <- gse133101 |>
  mutate(
    Gene_ID = sub(
      "\\..*$",
      "",
      as.character(Gene_ID)
    )
  )

gse136666 <- gse136666 |>
  mutate(
    Gene_ID = sub(
      "\\..*$",
      "",
      as.character(Gene_ID)
    )
  )

gse168496 <- gse168496 |>
  mutate(
    Gene_ID = sub(
      "\\..*$",
      "",
      as.character(Gene_ID)
    )
  )


# ------------------------------------------------------------
# 5. Apply DEG cutoff
# ------------------------------------------------------------

deg133101 <- gse133101 |>
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > lfc_cutoff
  )

deg136666 <- gse136666 |>
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > lfc_cutoff
  )

deg168496 <- gse168496 |>
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > lfc_cutoff
  )


cat("\nDEG counts after filtering:\n")
cat("GSE133101:", nrow(deg133101), "\n")
cat("GSE136666:", nrow(deg136666), "\n")
cat("GSE168496:", nrow(deg168496), "\n")


# ------------------------------------------------------------
# 6. Find pairwise common genes
#
# GSE136666 ∩ GSE168496
# ------------------------------------------------------------

common_ids <- intersect(
  deg136666$Gene_ID,
  deg168496$Gene_ID
)


cat(
  "\nGSE136666 ∩ GSE168496:",
  length(common_ids),
  "genes\n"
)


# ------------------------------------------------------------
# 7. Find GSE133101-specific genes
#
# Genes significant in GSE133101 but NOT in either
# of the other two datasets.
# ------------------------------------------------------------

specific_133101_ids <- setdiff(
  deg133101$Gene_ID,
  union(
    deg136666$Gene_ID,
    deg168496$Gene_ID
  )
)


cat(
  "GSE133101-specific genes:",
  length(specific_133101_ids),
  "\n"
)


# ------------------------------------------------------------
# 8. Read existing annotation
#
# Use GSE136666 annotation for the common genes
# and GSE133101 annotation for its specific genes.
# ------------------------------------------------------------

annotation_136666 <- import(
  file.path(
    annotation_dir,
    "GSE136666_Ensembl_annotation.csv"
  )
)

annotation_133101 <- import(
  file.path(
    annotation_dir,
    "GSE133101_Ensembl_annotation.csv"
  )
)


# Standardize GeneID

annotation_136666 <- annotation_136666 |>
  mutate(
    GeneID = sub(
      "\\..*$",
      "",
      as.character(GeneID)
    )
  )

annotation_133101 <- annotation_133101 |>
  mutate(
    GeneID = sub(
      "\\..*$",
      "",
      as.character(GeneID)
    )
  )


# ------------------------------------------------------------
# 9. Annotate the 9 common genes
# ------------------------------------------------------------

common_annotation <- annotation_136666 |>
  filter(
    GeneID %in% common_ids
  ) |>
  select(
    GeneID,
    Gene_Symbol,
    Gene_Biotype,
    Chromosome,
    Start,
    End,
    Strand
  ) |>
  rename(
    Gene_ID = GeneID
  )


# ------------------------------------------------------------
# 10. Add DESeq2 results for both datasets
# ------------------------------------------------------------

common_final <- common_annotation |>
  
  left_join(
    deg136666 |>
      filter(
        Gene_ID %in% common_ids
      ) |>
      select(
        Gene_ID,
        log2FoldChange,
        padj
      ) |>
      rename(
        log2FC_GSE136666 = log2FoldChange,
        padj_GSE136666 = padj
      ),
    by = "Gene_ID"
  ) |>
  
  left_join(
    deg168496 |>
      filter(
        Gene_ID %in% common_ids
      ) |>
      select(
        Gene_ID,
        log2FoldChange,
        padj
      ) |>
      rename(
        log2FC_GSE168496 = log2FoldChange,
        padj_GSE168496 = padj
      ),
    by = "Gene_ID"
  ) |>
  
  mutate(
    
    Direction_GSE136666 = case_when(
      log2FC_GSE136666 > 0 ~ "Up",
      log2FC_GSE136666 < 0 ~ "Down",
      TRUE ~ NA_character_
    ),
    
    Direction_GSE168496 = case_when(
      log2FC_GSE168496 > 0 ~ "Up",
      log2FC_GSE168496 < 0 ~ "Down",
      TRUE ~ NA_character_
    ),
    
    Direction = case_when(
      
      Direction_GSE136666 == "Up" &
        Direction_GSE168496 == "Up" ~
        "Consistently Up",
      
      Direction_GSE136666 == "Down" &
        Direction_GSE168496 == "Down" ~
        "Consistently Down",
      
      Direction_GSE136666 != Direction_GSE168496 ~
        "Discordant",
      
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 11. Save common genes
# ------------------------------------------------------------

export(
  common_final,
  file.path(
    output_dir,
    "GSE136666_GSE168496_common_genes_annotated.csv"
  )
)


# ------------------------------------------------------------
# 12. Annotate the 4 GSE133101-specific genes
# ------------------------------------------------------------

specific_133101 <- annotation_133101 |>
  filter(
    GeneID %in% specific_133101_ids
  ) |>
  select(
    GeneID,
    Gene_Symbol,
    Gene_Biotype,
    Chromosome,
    Start,
    End,
    Strand
  ) |>
  rename(
    Gene_ID = GeneID
  ) |>
  
  left_join(
    
    deg133101 |>
      filter(
        Gene_ID %in% specific_133101_ids
      ) |>
      select(
        Gene_ID,
        log2FoldChange,
        padj
      ) |>
      rename(
        log2FC_GSE133101 = log2FoldChange,
        padj_GSE133101 = padj
      ),
    
    by = "Gene_ID"
  ) |>
  
  mutate(
    
    Direction_GSE133101 = case_when(
      
      log2FC_GSE133101 > 0 ~
        "Up",
      
      log2FC_GSE133101 < 0 ~
        "Down",
      
      TRUE ~
        NA_character_
    )
  )


# ------------------------------------------------------------
# 13. Save GSE133101-specific genes
# ------------------------------------------------------------

export(
  specific_133101,
  file.path(
    output_dir,
    "GSE133101_specific_DEGs_annotated.csv"
  )
)


# ------------------------------------------------------------
# 14. Print results
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FINAL RECURRENT DEG RESULTS\n")
cat("============================================================\n")

cat(
  "\nGSE136666 ∩ GSE168496:",
  nrow(common_final),
  "genes\n"
)

cat(
  "GSE133101-specific:",
  nrow(specific_133101),
  "genes\n"
)


cat("\nCommon genes:\n")
print(
  common_final |>
    select(
      Gene_ID,
      Gene_Symbol,
      log2FC_GSE136666,
      log2FC_GSE168496,
      Direction
    )
)


cat("\nGSE133101-specific genes:\n")
print(
  specific_133101 |>
    select(
      Gene_ID,
      Gene_Symbol,
      log2FC_GSE133101,
      Direction_GSE133101
    )
)


cat("\n============================================================\n")
cat("OUTPUT FILES\n")
cat("============================================================\n")

cat(
  "\nCommon genes:\n",
  file.path(
    output_dir,
    "GSE136666_GSE168496_common_genes_annotated.csv"
  ),
  "\n"
)

cat(
  "\nGSE133101-specific genes:\n",
  file.path(
    output_dir,
    "GSE133101_specific_DEGs_annotated.csv"
  ),
  "\n"
)

cat("\n============================================================\n")
cat("COMPLETED\n")
cat("============================================================\n")