# ============================================================
# Parkinson's Disease Meta-analysis Functional Enrichment
#
# Input:
#   Random-effects meta-analysis results
#
# Analysis:
#   1. ORA
#      - GO Biological Process
#      - GO Molecular Function
#      - GO Cellular Component
#      - KEGG
#      - Reactome
#
#   2. GSEA
#      - GO Biological Process
#      - KEGG
#      - MSigDB Hallmark
#
# Gene identifier in meta-analysis:
#   Ensembl Gene ID
#
# Conversion:
#   Ensembl -> Entrez
#
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
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
distinct  <- dplyr::distinct
left_join <- dplyr::left_join


# ------------------------------------------------------------
# 3. Analysis thresholds
# ------------------------------------------------------------

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1


# ------------------------------------------------------------
# 4. Input / output directories
# ------------------------------------------------------------

BASE_DIR <- "/media/nabil/ssd/Parkinson-s"

META_FILE <- file.path(
  BASE_DIR,
  "results/tables/meta-analysis",
  "Parkinson_random_effects_meta_analysis.csv"
)

OUT_CSV <- file.path(
  BASE_DIR,
  "results/tables/enrichment"
)

OUT_FIG <- file.path(
  BASE_DIR,
  "results/figures/enrichment"
)

PREFIX <- "Parkinson"


dir.create(
  OUT_CSV,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  OUT_FIG,
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# 5. Load meta-analysis results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PARKINSON'S DISEASE FUNCTIONAL ENRICHMENT\n")
cat("============================================================\n\n")


meta <- import(META_FILE)


cat(
  "Genes in meta-analysis result: ",
  nrow(meta),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "Gene_ID",
  "randomSummary",
  "randomP"
)


missing_columns <- setdiff(
  required_columns,
  colnames(meta)
)


if (length(missing_columns) > 0) {
  
  stop(
    "Missing required column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
  
}


# ------------------------------------------------------------
# 7. Standardize Ensembl Gene IDs
# ------------------------------------------------------------

meta <- meta |>
  
  mutate(
    
    Gene_ID = as.character(Gene_ID),
    
    # Remove Ensembl version suffix
    # Example:
    # ENSG00000123456.7 -> ENSG00000123456
    
    Gene_ID = sub(
      "\\..*$",
      "",
      Gene_ID
    ),
    
    randomSummary =
      as.numeric(randomSummary),
    
    randomP =
      as.numeric(randomP)
    
  )


# ------------------------------------------------------------
# 8. Calculate BH-adjusted meta-analysis p-values
#
# MetaVolcanoR gives randomP.
# We calculate BH-adjusted values across the
# complete meta-analysis result.
# ------------------------------------------------------------

meta <- meta |>
  
  mutate(
    
    randomPadj =
      p.adjust(
        randomP,
        method = "BH"
      )
    
  )


cat(
  "Genes with valid meta-analysis p-values: ",
  sum(!is.na(meta$randomP)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Convert Ensembl -> Entrez
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONVERTING ENSEMBL IDs TO ENTREZ IDs\n")
cat("============================================================\n\n")


id_conversion <- tryCatch(
  
  bitr(
    unique(meta$Gene_ID),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ),
  
  error = function(e) {
    
    stop(
      "Ensembl -> Entrez conversion failed:\n",
      e$message
    )
    
  }
  
)


# Remove duplicate Ensembl IDs
id_conversion <- id_conversion |>
  
  distinct(
    ENSEMBL,
    .keep_all = TRUE
  )


cat(
  "Ensembl IDs available: ",
  length(unique(meta$Gene_ID)),
  "\n",
  sep = ""
)

cat(
  "Ensembl IDs successfully converted: ",
  nrow(id_conversion),
  "\n",
  sep = ""
)

cat(
  "Conversion rate: ",
  round(
    nrow(id_conversion) /
      length(unique(meta$Gene_ID)) *
      100,
    2
  ),
  "%\n",
  sep = ""
)


# Save ID conversion table
export(
  id_conversion,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_Ensembl_to_Entrez_conversion.csv"
    )
  )
)


# ------------------------------------------------------------
# 10. Add Entrez IDs to meta-analysis table
# ------------------------------------------------------------

meta <- meta |>
  
  left_join(
    
    id_conversion,
    
    by = c(
      "Gene_ID" = "ENSEMBL"
    )
    
  )


# ------------------------------------------------------------
# 11. Define background universe
#
# All genes that were actually tested in the meta-analysis
# and successfully converted to Entrez IDs.
# ------------------------------------------------------------

universe_entrez <- meta |>
  
  filter(
    !is.na(ENTREZID)
  ) |>
  
  pull(
    ENTREZID
  ) |>
  
  unique()


cat(
  "\nBackground Entrez gene universe: ",
  length(universe_entrez),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Define significant meta-analysis DEGs
#
# Criteria:
#
#   BH-adjusted meta-analysis P < 0.05
#   |randomSummary| >= 1
#
# randomSummary = pooled log2 fold change
# ------------------------------------------------------------

sig <- meta |>
  
  filter(
    
    !is.na(randomPadj),
    
    randomPadj < PADJ_CUTOFF,
    
    !is.na(randomSummary),
    
    abs(randomSummary) >= LFC_CUTOFF,
    
    !is.na(ENTREZID)
    
  )


cat(
  "\nSignificant meta-DEGs: ",
  nrow(sig),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Separate UP and DOWN genes
# ------------------------------------------------------------

entrez_up <- sig |>
  
  filter(
    randomSummary > 0
  ) |>
  
  pull(
    ENTREZID
  ) |>
  
  unique()


entrez_down <- sig |>
  
  filter(
    randomSummary < 0
  ) |>
  
  pull(
    ENTREZID
  ) |>
  
  unique()


entrez_all <- sig |>
  
  pull(
    ENTREZID
  ) |>
  
  unique()


cat(
  "Upregulated genes: ",
  length(entrez_up),
  "\n",
  sep = ""
)

cat(
  "Downregulated genes: ",
  length(entrez_down),
  "\n",
  sep = ""
)


# ============================================================
# 14. ORA FUNCTION
# ============================================================

run_ora <- function(
    entrez_ids,
    universe_entrez
) {
  
  if (length(entrez_ids) < 5) {
    
    message(
      "Too few genes for ORA: n = ",
      length(entrez_ids)
    )
    
    return(NULL)
    
  }
  
  
  results <- list()
  
  
  # ----------------------------------------------------------
  # GO Biological Process
  # ----------------------------------------------------------
  
  results[["GO_BP"]] <- tryCatch(
    
    enrichGO(
      
      gene = entrez_ids,
      
      universe = universe_entrez,
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "BP",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.20,
      
      readable = TRUE
      
    ),
    
    error = function(e) {
      
      message(
        "GO-BP error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # GO Molecular Function
  # ----------------------------------------------------------
  
  results[["GO_MF"]] <- tryCatch(
    
    enrichGO(
      
      gene = entrez_ids,
      
      universe = universe_entrez,
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "MF",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.20,
      
      readable = TRUE
      
    ),
    
    error = function(e) {
      
      message(
        "GO-MF error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # GO Cellular Component
  # ----------------------------------------------------------
  
  results[["GO_CC"]] <- tryCatch(
    
    enrichGO(
      
      gene = entrez_ids,
      
      universe = universe_entrez,
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "CC",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.20,
      
      readable = TRUE
      
    ),
    
    error = function(e) {
      
      message(
        "GO-CC error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # KEGG
  # ----------------------------------------------------------
  
  results[["KEGG"]] <- tryCatch(
    
    enrichKEGG(
      
      gene = entrez_ids,
      
      universe = universe_entrez,
      
      organism = "hsa",
      
      keyType = "kegg",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05
      
    ),
    
    error = function(e) {
      
      message(
        "KEGG error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # Reactome
  # ----------------------------------------------------------
  
  results[["Reactome"]] <- tryCatch(
    
    enrichPathway(
      
      gene = entrez_ids,
      
      universe = universe_entrez,
      
      organism = "human",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      readable = TRUE
      
    ),
    
    error = function(e) {
      
      message(
        "Reactome error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  return(results)
  
}


# ============================================================
# 15. RUN ORA
# ============================================================

cat("\n============================================================\n")
cat("RUNNING ORA\n")
cat("============================================================\n\n")


ora_up <- run_ora(
  entrez_up,
  universe_entrez
)


ora_down <- run_ora(
  entrez_down,
  universe_entrez
)


ora_all <- run_ora(
  entrez_all,
  universe_entrez
)


# ============================================================
# 16. EXPORT ORA RESULTS
# ============================================================

export_enrichment <- function(
    result_list,
    out_dir,
    prefix
) {
  
  for (name in names(result_list)) {
    
    res <- result_list[[name]]
    
    if (is.null(res)) {
      next
    }
    
    df <- tryCatch(
      as.data.frame(res),
      error = function(e) NULL
    )
    
    if (
      is.null(df) ||
      nrow(df) == 0
    ) {
      next
    }
    
    export(
      df,
      file.path(
        out_dir,
        paste0(
          prefix,
          "_",
          name,
          ".csv"
        )
      )
    )
    
  }
  
}


export_enrichment(
  ora_up,
  OUT_CSV,
  paste0(PREFIX, "_ORA_UP")
)


export_enrichment(
  ora_down,
  OUT_CSV,
  paste0(PREFIX, "_ORA_DOWN")
)


export_enrichment(
  ora_all,
  OUT_CSV,
  paste0(PREFIX, "_ORA_ALL")
)


# ============================================================
# 17. BUILD GSEA RANKED LIST
#
# Ranking:
#
# sign(randomSummary) * -log10(randomP)
#
# Positive values = Parkinson's-up
# Negative values = Parkinson's-down
#
# IMPORTANT:
# GSEA uses the full meta-analysis gene list,
# NOT only significant DEGs.
# ============================================================

build_ranked_list <- function(meta_df) {
  
  df <- meta_df |>
    
    filter(
      
      !is.na(randomP),
      
      randomP > 0,
      
      !is.na(randomSummary),
      
      !is.na(ENTREZID)
      
    ) |>
    
    mutate(
      
      rank_score =
        sign(randomSummary) *
        -log10(randomP),
      
      ENTREZID =
        as.character(ENTREZID)
      
    ) |>
    
    filter(
      !duplicated(ENTREZID)
    ) |>
    
    arrange(
      desc(rank_score)
    )
  
  
  ranked <- setNames(
    df$rank_score,
    df$ENTREZID
  )
  
  
  ranked
  
}


# ============================================================
# 18. RUN GSEA
# ============================================================

run_gsea <- function(
    ranked_list
) {
  
  if (
    length(ranked_list) < 100
  ) {
    
    message(
      "Too few genes for GSEA: n = ",
      length(ranked_list)
    )
    
    return(NULL)
    
  }
  
  
  results <- list()
  
  
  # ----------------------------------------------------------
  # GO Biological Process
  # ----------------------------------------------------------
  
  results[["GSEA_GO_BP"]] <- tryCatch(
    
    gseGO(
      
      geneList = ranked_list,
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "BP",
      
      minGSSize = 10,
      
      maxGSSize = 500,
      
      pvalueCutoff = 0.05,
      
      verbose = FALSE
      
    ),
    
    error = function(e) {
      
      message(
        "GSEA GO error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # KEGG
  # ----------------------------------------------------------
  
  results[["GSEA_KEGG"]] <- tryCatch(
    
    gseKEGG(
      
      geneList = ranked_list,
      
      organism = "hsa",
      
      minGSSize = 10,
      
      maxGSSize = 500,
      
      pvalueCutoff = 0.05,
      
      verbose = FALSE
      
    ),
    
    error = function(e) {
      
      message(
        "GSEA KEGG error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # MSigDB Hallmark
  # ----------------------------------------------------------
  
  hallmark <- msigdbr(
    
    species = "Homo sapiens",
    
    collection = "H"
    
  ) |>
    
    select(
      gs_name,
      ncbi_gene
    ) |>
    
    mutate(
      ncbi_gene =
        as.character(ncbi_gene)
    )
  
  
  results[["GSEA_Hallmark"]] <- tryCatch(
    
    GSEA(
      
      geneList = ranked_list,
      
      TERM2GENE = hallmark,
      
      minGSSize = 10,
      
      maxGSSize = 500,
      
      pvalueCutoff = 0.05,
      
      verbose = FALSE
      
    ),
    
    error = function(e) {
      
      message(
        "GSEA Hallmark error: ",
        e$message
      )
      
      NULL
      
    }
    
  )
  
  
  return(results)
  
}


# ============================================================
# 19. RUN GSEA
# ============================================================

cat("\n============================================================\n")
cat("RUNNING GSEA\n")
cat("============================================================\n\n")


ranked <- build_ranked_list(
  meta
)


cat(
  "Genes ranked for GSEA: ",
  length(ranked),
  "\n",
  sep = ""
)


gsea_res <- run_gsea(
  ranked
)


# ============================================================
# 20. EXPORT GSEA RESULTS
# ============================================================

export_enrichment(
  gsea_res,
  OUT_CSV,
  paste0(PREFIX, "_GSEA")
)


# ============================================================
# 21. Dotplot function
# ============================================================

save_dotplot <- function(
    res,
    title,
    path,
    width = 10,
    height = 8,
    show = 20
) {
  
  if (is.null(res)) {
    return(invisible(NULL))
  }
  
  
  df <- tryCatch(
    as.data.frame(res),
    error = function(e) NULL
  )
  
  
  if (
    is.null(df) ||
    nrow(df) == 0
  ) {
    return(invisible(NULL))
  }
  
  
  # GSEA result
  if (
    all(
      c(
        "NES",
        "setSize"
      ) %in% colnames(df)
    )
  ) {
    
    df <- df |>
      arrange(
        p.adjust
      ) |>
      slice_head(
        n = show
      )
    
    
    df$Label <- stringr::str_wrap(
      df$Description,
      35
    )
    
    
    p <- ggplot(
      df,
      aes(
        x = NES,
        y = reorder(
          Label,
          NES
        )
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey60"
      ) +
      
      geom_point(
        aes(
          size = setSize,
          fill = p.adjust
        ),
        shape = 21,
        color = "grey30"
      ) +
      
      scale_fill_gradient(
        low = "#b2182b",
        high = "#2166ac",
        trans = "log10",
        name = "Adjusted P"
      ) +
      
      scale_size_continuous(
        range = c(
          3,
          10
        ),
        name = "Gene set size"
      ) +
      
      labs(
        title = title,
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      ) +
      
      theme_bw(
        base_size = 13
      ) +
      
      theme(
        plot.title =
          element_text(
            face = "bold",
            hjust = 0.5
          ),
        panel.grid.minor =
          element_blank()
      )
    
    
  } else {
    
    p <- dotplot(
      res,
      showCategory = show
    ) +
      
      ggtitle(
        title
      )
    
  }
  
  
  ggsave(
    path,
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  
  invisible(p)
  
}


# ============================================================
# 22. Save ORA dotplots
# ============================================================

for (
  name in names(ora_all)
) {
  
  save_dotplot(
    
    ora_all[[name]],
    
    paste(
      PREFIX,
      "ORA",
      name
    ),
    
    file.path(
      OUT_FIG,
      paste0(
        PREFIX,
        "_ORA_ALL_",
        name,
        "_dotplot.png"
      )
    )
    
  )
  
}


# ============================================================
# 23. Save GSEA dotplots
# ============================================================

for (
  name in names(gsea_res)
) {
  
  save_dotplot(
    
    gsea_res[[name]],
    
    paste(
      PREFIX,
      name
    ),
    
    file.path(
      OUT_FIG,
      paste0(
        PREFIX,
        "_",
        name,
        "_dotplot.png"
      )
    )
    
  )
  
}


# ============================================================
# 24. GSEA enrichment plots
# ============================================================

save_gsea_plot <- function(
    res,
    title,
    path,
    n = 3
) {
  
  if (is.null(res)) {
    return(invisible(NULL))
  }
  
  
  df <- tryCatch(
    as.data.frame(res),
    error = function(e) NULL
  )
  
  
  if (
    is.null(df) ||
    nrow(df) == 0
  ) {
    return(invisible(NULL))
  }
  
  
  ids <- head(
    df$ID,
    n
  )
  
  
  p <- gseaplot2(
    res,
    geneSetID = ids,
    title = title
  )
  
  
  ggsave(
    path,
    plot = p,
    width = 12,
    height = 6 * ceiling(n / 2),
    dpi = 300,
    bg = "white"
  )
  
  
  invisible(p)
  
}


for (
  name in names(gsea_res)
) {
  
  save_gsea_plot(
    
    gsea_res[[name]],
    
    paste(
      PREFIX,
      name,
      "- Top Pathways"
    ),
    
    file.path(
      OUT_FIG,
      paste0(
        PREFIX,
        "_",
        name,
        "_enrichment_plot.png"
      )
    )
    
  )
  
}


# ============================================================
# 25. Save annotated significant meta-DEGs
# ============================================================

sig_output <- sig |>
  
  select(
    Gene_ID,
    ENTREZID,
    any_of(
      c(
        "Gene_Symbol",
        "Gene_Biotype",
        "Chromosome",
        "Start",
        "End",
        "Strand"
      )
    ),
    randomSummary,
    randomP,
    randomPadj
  ) |>
  
  arrange(
    randomPadj
  )


export(
  sig_output,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_significant_meta_DEGs_for_enrichment.csv"
    )
  )
)


# ============================================================
# 26. Final summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PARKINSON'S ENRICHMENT ANALYSIS COMPLETED\n")
cat("============================================================\n\n")


cat(
  "Meta-analysis genes: ",
  nrow(meta),
  "\n",
  sep = ""
)


cat(
  "Genes converted to Entrez: ",
  length(
    unique(
      meta$ENTREZID[
        !is.na(meta$ENTREZID)
      ]
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Significant meta-DEGs: ",
  nrow(sig),
  "\n",
  sep = ""
)


cat(
  "Upregulated: ",
  length(entrez_up),
  "\n",
  sep = ""
)


cat(
  "Downregulated: ",
  length(entrez_down),
  "\n",
  sep = ""
)


cat(
  "GSEA ranked genes: ",
  length(ranked),
  "\n",
  sep = ""
)


cat("\nOutput directory:\n")
cat(
  OUT_CSV,
  "\n"
)

cat(
  OUT_FIG,
  "\n"
)


cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")