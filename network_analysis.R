# ==============================================================================
# PPI network analysis of colorectal cancer meta-analysis
# STRING REST API + igraph
#
# ANALYSIS LOGIC
# ------------------------------------------------------------------------------
# 1. Try strict meta-analysis DEGs:
#       BH-adjusted randomP < 0.05
#       |randomSummary| >= 1
#
# 2. If >= MIN_STRICT_GENES are available:
#       Build a STRICT DEG-based PPI network.
#
# 3. If fewer than MIN_STRICT_GENES are available:
#       Automatically switch to an EXPLORATORY network using the top
#       TOP_EXPLORATORY_GENES genes ranked by meta-analysis P value,
#       while retaining |log2FC| >= 1.
#
# IMPORTANT:
# Exploratory networks are NOT equivalent to statistically significant
# DEG networks. They should be described as exploratory PPI analysis.
# ==============================================================================


# ==============================================================================
# 01. PACKAGES
# ==============================================================================

library(httr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidyverse)
library(rio)
library(scales)


# Pin dplyr verbs in case another package masks them
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize
arrange   <- dplyr::arrange
rename    <- dplyr::rename
count     <- dplyr::count
desc      <- dplyr::desc
select    <- dplyr::select
filter    <- dplyr::filter


# ==============================================================================
# 02. PARAMETERS
# ==============================================================================

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

# Minimum number of genes required for strict PPI analysis
MIN_STRICT_GENES <- 5

# Number of genes to use if exploratory mode is triggered
TOP_EXPLORATORY_GENES <- 100

# STRING combined score
# 700 = high confidence
STRING_SCORE <- 700

STRING_SPECIES <- 9606

# STRING API practical node limit
STRING_MAX_NODES <- 2000

# Top 10% degree = hub genes
HUB_QUANTILE <- 0.90

# Minimum module size for plotting
MIN_MODULE_SIZE <- 5


# ==============================================================================
# 03. INPUT / OUTPUT
# ==============================================================================

META_FILE <- "results/tables/meta-analysis/Parkinson_random_effects_meta_analysis.csv"

OUT_CSV <- "results/tables/network"
OUT_FIG <- "results/figures/network"

PREFIX <- "Parkinson"

dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 04. CHECK INPUT FILE
# ==============================================================================

if (!file.exists(META_FILE)) {
  stop(
    "\nMeta-analysis file not found:\n",
    META_FILE,
    "\n\nMake sure the random-effect meta-analysis has been successfully run first."
  )
}


# ==============================================================================
# 05. STRING REST API
# ==============================================================================

string_api_network <- function(
    symbols,
    required_score = STRING_SCORE,
    species = STRING_SPECIES
) {
  
  symbols <- unique(
    symbols[
      !is.na(symbols) &
        symbols != ""
    ]
  )
  
  if (length(symbols) < 2) {
    return(NULL)
  }
  
  message(
    "  Sending ",
    length(symbols),
    " genes to STRING..."
  )
  
  resp <- tryCatch(
    
    httr::POST(
      "https://string-db.org/api/tsv/network",
      
      body = list(
        identifiers = paste(symbols, collapse = "\n"),
        species = species,
        required_score = required_score,
        caller_identity = "crc_meta_pipeline"
      ),
      
      encode = "form",
      
      httr::timeout(120)
    ),
    
    error = function(e) {
      message(
        "  STRING API error: ",
        e$message
      )
      NULL
    }
  )
  
  if (is.null(resp) || httr::http_error(resp)) {
    
    message(
      "  STRING API request failed."
    )
    
    return(NULL)
  }
  
  txt <- httr::content(
    resp,
    as = "text",
    encoding = "UTF-8"
  )
  
  df <- tryCatch(
    
    utils::read.delim(
      text = txt,
      stringsAsFactors = FALSE
    ),
    
    error = function(e) NULL
  )
  
  if (
    is.null(df) ||
    !all(
      c(
        "preferredName_A",
        "preferredName_B",
        "score"
      ) %in% names(df)
    )
  ) {
    
    message(
      "  STRING returned no usable network."
    )
    
    return(NULL)
  }
  
  df
}


# ==============================================================================
# 06. BUILD IGRAPH NETWORK
# ==============================================================================

build_ppi_graph <- function(
    net_df,
    lfc_map
) {
  
  if (
    is.null(net_df) ||
    nrow(net_df) == 0
  ) {
    return(NULL)
  }
  
  edges <- net_df |>
    
    transmute(
      gene_from = preferredName_A,
      gene_to = preferredName_B,
      combined_score = score
    ) |>
    
    filter(
      gene_from %in% names(lfc_map),
      gene_to %in% names(lfc_map)
    )
  
  if (nrow(edges) == 0) {
    return(NULL)
  }
  
  nodes <- unique(
    c(
      edges$gene_from,
      edges$gene_to
    )
  )
  
  vertex_df <- tibble(
    Symbol = nodes,
    pooled_log2FC = unname(
      lfc_map[nodes]
    )
  )
  
  g <- graph_from_data_frame(
    d = edges,
    vertices = vertex_df,
    directed = FALSE
  )
  
  g <- simplify(g)
  
  # Remove isolated nodes
  iso <- which(
    degree(g) == 0
  )
  
  if (length(iso) > 0) {
    g <- delete_vertices(
      g,
      iso
    )
  }
  
  g
}


# ==============================================================================
# 07. CENTRALITY METRICS
# ==============================================================================

compute_metrics <- function(g) {
  
  tibble(
    
    Symbol = V(g)$name,
    
    Degree = degree(g),
    
    Betweenness =
      betweenness(
        g,
        normalized = TRUE
      ),
    
    Closeness =
      closeness(
        g,
        normalized = TRUE
      ),
    
    Eigenvector =
      eigen_centrality(g)$vector,
    
    Log2FC =
      V(g)$pooled_log2FC
    
  ) |>
    
    arrange(
      desc(Degree)
    )
}


# ==============================================================================
# 08. HUB GENES
# ==============================================================================

get_hub_genes <- function(
    metrics_df,
    quantile_cut = HUB_QUANTILE
) {
  
  if (
    is.null(metrics_df) ||
    nrow(metrics_df) == 0
  ) {
    return(metrics_df)
  }
  
  thresh <- quantile(
    metrics_df$Degree,
    quantile_cut,
    na.rm = TRUE
  )
  
  metrics_df |>
    filter(
      Degree >= thresh
    )
}


# ==============================================================================
# 09. COMMUNITY DETECTION
# ==============================================================================

detect_communities <- function(g) {
  
  comm <- cluster_louvain(g)
  
  V(g)$community <-
    membership(comm)
  
  list(
    graph = g,
    communities = comm
  )
}


# ==============================================================================
# 10. VISUAL ATTRIBUTES
# ==============================================================================

set_visual_attrs <- function(
    g,
    metrics_df,
    label_top = 15
) {
  
  deg <- degree(g)
  
  top_hubs <-
    names(
      sort(
        deg,
        decreasing = TRUE
      )
    )[
      seq_len(
        min(
          label_top,
          length(deg)
        )
      )
    ]
  
  V(g)$size <-
    rescale(
      deg,
      to = c(3, 16)
    )
  
  V(g)$color <-
    ifelse(
      V(g)$pooled_log2FC > 0,
      "#de2d26",
      "#2171b5"
    )
  
  V(g)$label <-
    ifelse(
      V(g)$name %in% top_hubs,
      V(g)$name,
      NA
    )
  
  g
}


# ==============================================================================
# 11. BASE IGRAPH PLOT
# ==============================================================================

save_network_plot <- function(
    g,
    title,
    path,
    seed = 42
) {
  
  if (
    is.null(g) ||
    vcount(g) < 2
  ) {
    return(invisible(NULL))
  }
  
  dir.create(
    dirname(path),
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  png(
    path,
    width = 12,
    height = 12,
    units = "in",
    res = 300
  )
  
  set.seed(seed)
  
  plot(
    g,
    
    vertex.size =
      V(g)$size,
    
    vertex.color =
      V(g)$color,
    
    vertex.label =
      V(g)$label,
    
    vertex.label.cex = 0.65,
    
    vertex.label.color = "black",
    
    edge.width = 0.4,
    
    edge.color = "gray75",
    
    layout =
      layout_with_fr(g),
    
    main = title
  )
  
  legend(
    "topright",
    
    legend =
      c(
        "Upregulated",
        "Downregulated"
      ),
    
    col =
      c(
        "#de2d26",
        "#2171b5"
      ),
    
    pch = 19,
    
    bty = "n",
    
    pt.cex = 1.5
  )
  
  dev.off()
  
  invisible(NULL)
}


# ==============================================================================
# 12. GGRAPH PLOT
# ==============================================================================

save_ggraph_plot <- function(
    g,
    metrics_df,
    title,
    path,
    seed = 42
) {
  
  if (
    is.null(g) ||
    vcount(g) < 2
  ) {
    return(invisible(NULL))
  }
  
  dir.create(
    dirname(path),
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  deg <- degree(g)
  
  top_hubs <-
    names(
      sort(
        deg,
        decreasing = TRUE
      )
    )[
      seq_len(
        min(
          15,
          length(deg)
        )
      )
    ]
  
  V(g)$deg <- deg
  
  V(g)$direction <-
    ifelse(
      V(g)$pooled_log2FC > 0,
      "Up",
      "Down"
    )
  
  V(g)$lab <-
    ifelse(
      V(g)$name %in% top_hubs,
      V(g)$name,
      NA_character_
    )
  
  set.seed(seed)
  
  p <-
    ggraph(
      g,
      layout = "fr"
    ) +
    
    geom_edge_link(
      alpha = 0.25,
      color = "gray60",
      width = 0.3
    ) +
    
    geom_node_point(
      aes(
        size = deg,
        color = direction
      ),
      alpha = 0.85
    ) +
    
    geom_node_text(
      aes(
        label = lab
      ),
      repel = TRUE,
      size = 3,
      max.overlaps = Inf,
      na.rm = TRUE,
      fontface = "italic"
    ) +
    
    scale_color_manual(
      values =
        c(
          Up = "#de2d26",
          Down = "#2171b5"
        ),
      name = "Direction"
    ) +
    
    scale_size_continuous(
      range = c(1.5, 9),
      name = "Degree"
    ) +
    
    labs(
      title = title
    ) +
    
    theme_graph(
      base_family = "sans"
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5,
          size = 14
        ),
      
      legend.position = "right"
    )
  
  ggsave(
    path,
    plot = p,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  invisible(p)
}


# ==============================================================================
# 13. HUB GENE BARPLOT
# ==============================================================================

save_hub_barplot <- function(
    hub_df,
    title,
    path,
    n = 25
) {
  
  if (
    is.null(hub_df) ||
    nrow(hub_df) == 0
  ) {
    return(invisible(NULL))
  }
  
  dir.create(
    dirname(path),
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  df <- head(
    hub_df,
    n
  )
  
  p <-
    ggplot(
      df,
      aes(
        x = reorder(Symbol, Degree),
        y = Degree,
        fill = Log2FC
      )
    ) +
    
    geom_bar(
      stat = "identity"
    ) +
    
    scale_fill_gradient2(
      low = "#2171b5",
      mid = "white",
      high = "#de2d26",
      midpoint = 0,
      name = "Log2FC"
    ) +
    
    coord_flip() +
    
    labs(
      title = title,
      x = NULL,
      y = "Degree Centrality"
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        )
    )
  
  ggsave(
    path,
    plot = p,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  invisible(p)
}


# ==============================================================================
# 14. MODULE PLOTS
# ==============================================================================

save_module_plots <- function(
    g,
    fig_dir,
    prefix,
    n_modules = 5,
    seed = 42
) {
  
  if (
    is.null(g) ||
    is.null(V(g)$community)
  ) {
    return(invisible(NULL))
  }
  
  module_sizes <-
    sort(
      table(V(g)$community),
      decreasing = TRUE
    )
  
  top_ids <-
    names(
      head(
        module_sizes,
        n_modules
      )
    )
  
  for (mid in top_ids) {
    
    vids <-
      which(
        V(g)$community == mid
      )
    
    subg <-
      induced_subgraph(
        g,
        vids = vids
      )
    
    if (
      vcount(subg) < MIN_MODULE_SIZE
    ) {
      next
    }
    
    path <-
      file.path(
        fig_dir,
        paste0(
          prefix,
          "_module_",
          mid,
          "_network.png"
        )
      )
    
    png(
      path,
      width = 8,
      height = 8,
      units = "in",
      res = 300
    )
    
    set.seed(seed)
    
    sub_deg <-
      degree(subg)
    
    plot(
      subg,
      
      vertex.size =
        rescale(
          sub_deg,
          to = c(5, 18)
        ),
      
      vertex.color =
        V(subg)$color,
      
      vertex.label =
        V(subg)$name,
      
      vertex.label.cex = 0.6,
      
      vertex.label.color = "black",
      
      edge.width = 1,
      
      edge.color = "gray70",
      
      layout =
        layout_with_fr(subg),
      
      main =
        paste0(
          prefix,
          " - Module ",
          mid,
          " (",
          vcount(subg),
          " genes)"
        )
    )
    
    dev.off()
  }
  
  invisible(NULL)
}


# ==============================================================================
# 15. EXPORT NETWORK RESULTS
# ==============================================================================

export_network <- function(
    metrics_df,
    hub_df,
    modules_df,
    out_dir,
    prefix
) {
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  export(
    metrics_df,
    file.path(
      out_dir,
      paste0(
        prefix,
        "_network_metrics.csv"
      )
    )
  )
  
  export(
    hub_df,
    file.path(
      out_dir,
      paste0(
        prefix,
        "_hub_genes.csv"
      )
    )
  )
  
  export(
    modules_df,
    file.path(
      out_dir,
      paste0(
        prefix,
        "_network_modules.csv"
      )
    )
  )
}


# ==============================================================================
# 16. LOAD META-ANALYSIS RESULTS
# ==============================================================================

message(
  "\n============================================================"
)

message(
  "COLORECTAL CANCER META-ANALYSIS PPI"
)

message(
  "============================================================"
)

meta <-
  import(META_FILE)


# Convert Gene_ID to character
meta <-
  meta |>
  mutate(
    Gene_ID =
      as.character(Gene_ID)
  )


message(
  "Genes in meta-analysis: ",
  nrow(meta)
)


# ==============================================================================
# 17. DERIVE BH-ADJUSTED P VALUES
# ==============================================================================

meta <-
  meta |>
  mutate(
    rem_padj =
      p.adjust(
        randomP,
        method = "BH"
      )
  )


# ==============================================================================
# 18. STRICT DEG SET
# ==============================================================================

strict_sig <-
  meta |>
  
  filter(
    !is.na(rem_padj),
    
    rem_padj <
      PADJ_CUTOFF,
    
    !is.na(randomSummary),
    
    abs(randomSummary) >=
      LFC_CUTOFF,
    
    !is.na(Gene_Symbol),
    
    Gene_Symbol != ""
  ) |>
  
  arrange(
    randomP,
    desc(
      abs(randomSummary)
    )
  ) |>
  
  distinct(
    Gene_Symbol,
    .keep_all = TRUE
  ) |>
  
  transmute(
    Symbol =
      Gene_Symbol,
    
    pooled_log2FC =
      randomSummary,
    
    randomP =
      randomP,
    
    rem_padj =
      rem_padj
  )


message(
  "\nStrict significant DEGs: ",
  nrow(strict_sig)
)


# ==============================================================================
# 19. AUTOMATIC MODE SELECTION
# ==============================================================================

if (
  nrow(strict_sig) >=
  MIN_STRICT_GENES
) {
  
  # --------------------------------------------------------------------------
  # STRICT MODE
  # --------------------------------------------------------------------------
  
  ANALYSIS_MODE <- "strict"
  
  sig <- strict_sig
  
  message(
    "\n>>> STRICT PPI MODE"
  )
  
  message(
    "Using ",
    nrow(sig),
    " statistically significant meta-analysis DEGs."
  )
  
} else {
  
  # --------------------------------------------------------------------------
  # EXPLORATORY MODE
  # --------------------------------------------------------------------------
  
  ANALYSIS_MODE <- "exploratory"
  
  message(
    "\n>>> EXPLORATORY PPI MODE"
  )
  
  message(
    "Only ",
    nrow(strict_sig),
    " genes passed the strict DEG criteria."
  )
  
  message(
    "Switching to top ",
    TOP_EXPLORATORY_GENES,
    " meta-analysis genes."
  )
  
  # Candidate genes:
  #   - must have a gene symbol
  #   - must have valid randomP
  #   - must have |log2FC| >= LFC_CUTOFF
  #
  # Rank primarily by randomP and secondarily by effect size.
  
  exploratory_candidates <-
    meta |>
    
    filter(
      !is.na(randomP),
      randomP > 0,
      
      !is.na(randomSummary),
      
      abs(randomSummary) >=
        LFC_CUTOFF,
      
      !is.na(Gene_Symbol),
      
      Gene_Symbol != ""
    ) |>
    
    arrange(
      randomP,
      desc(
        abs(randomSummary)
      )
    ) |>
    
    distinct(
      Gene_Symbol,
      .keep_all = TRUE
    ) |>
    
    transmute(
      Symbol =
        Gene_Symbol,
      
      pooled_log2FC =
        randomSummary,
      
      randomP =
        randomP,
      
      rem_padj =
        rem_padj
    )
  
  if (
    nrow(exploratory_candidates) < 5
  ) {
    
    stop(
      paste0(
        "\nNot enough genes for even exploratory PPI analysis.\n",
        
        "Genes meeting |log2FC| >= ",
        LFC_CUTOFF,
        ": ",
        nrow(exploratory_candidates),
        "\n\n",
        
        "Consider lowering LFC_CUTOFF or reviewing the meta-analysis."
      )
    )
  }
  
  sig <-
    head(
      exploratory_candidates,
      TOP_EXPLORATORY_GENES
    )
  
  message(
    "Exploratory genes selected: ",
    nrow(sig)
  )
}


# ==============================================================================
# 20. CAP TO STRING LIMIT
# ==============================================================================

if (
  nrow(sig) >
  STRING_MAX_NODES
) {
  
  message(
    "Capping network to ",
    STRING_MAX_NODES,
    " genes."
  )
  
  sig <-
    head(
      sig,
      STRING_MAX_NODES
    )
}


# ==============================================================================
# 21. SAVE THE GENE SET USED FOR PPI
# ==============================================================================

ppi_gene_set <-
  sig |>
  
  mutate(
    Analysis_Mode =
      ANALYSIS_MODE
  )

export(
  ppi_gene_set,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_ppi_gene_set_",
      ANALYSIS_MODE,
      ".csv"
    )
  )
)


# ==============================================================================
# 22. LFC MAP
# ==============================================================================

lfc_map <-
  setNames(
    sig$pooled_log2FC,
    sig$Symbol
  )


# ==============================================================================
# 23. STRING QUERY
# ==============================================================================

message(
  "\nQuerying STRING API..."
)

net_df <-
  string_api_network(
    sig$Symbol
  )


if (
  is.null(net_df)
) {
  
  stop(
    paste0(
      "\nSTRING API returned no usable network.\n",
      "Try lowering STRING_SCORE from ",
      STRING_SCORE,
      " to 500."
    )
  )
}


message(
  "Interactions returned: ",
  nrow(net_df)
)


# ==============================================================================
# 24. SAVE RAW STRING INTERACTIONS
# ==============================================================================

export(
  net_df,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_STRING_interactions_",
      ANALYSIS_MODE,
      ".csv"
    )
  )
)


# ==============================================================================
# 25. BUILD NETWORK
# ==============================================================================

message(
  "\nBuilding PPI network..."
)

g <-
  build_ppi_graph(
    net_df,
    lfc_map
  )


if (
  is.null(g)
) {
  
  stop(
    "\nNo STRING interactions found at the selected confidence threshold."
  )
}


message(
  "Network nodes: ",
  vcount(g)
)

message(
  "Network edges: ",
  ecount(g)
)


# ==============================================================================
# 26. CENTRALITY + HUBS
# ==============================================================================

metrics <-
  compute_metrics(g)

hub_genes <-
  get_hub_genes(metrics)


message(
  "Hub genes: ",
  nrow(hub_genes)
)


# ==============================================================================
# 27. COMMUNITY DETECTION
# ==============================================================================

comm_result <-
  detect_communities(g)

g <-
  comm_result$graph


modules_df <-
  tibble(
    
    Symbol =
      V(g)$name,
    
    Module =
      V(g)$community,
    
    Degree =
      degree(g),
    
    Log2FC =
      V(g)$pooled_log2FC
    
  ) |>
  
  arrange(
    Module,
    desc(Degree)
  )


# ==============================================================================
# 28. EXPORT NETWORK RESULTS
# ==============================================================================

export_network(
  metrics_df = metrics,
  hub_df = hub_genes,
  modules_df = modules_df,
  out_dir = OUT_CSV,
  prefix =
    paste0(
      PREFIX,
      "_",
      ANALYSIS_MODE
    )
)


# ==============================================================================
# 29. ANALYSIS SUMMARY
# ==============================================================================

analysis_summary <-
  tibble(
    
    Cancer =
      PREFIX,
    
    Analysis_Mode =
      ANALYSIS_MODE,
    
    Strict_DEGs =
      nrow(strict_sig),
    
    Genes_Used_for_PPI =
      nrow(sig),
    
    Padj_Cutoff =
      PADJ_CUTOFF,
    
    LFC_Cutoff =
      LFC_CUTOFF,
    
    STRING_Score =
      STRING_SCORE,
    
    Network_Nodes =
      vcount(g),
    
    Network_Edges =
      ecount(g),
    
    Hub_Quantile =
      HUB_QUANTILE,
    
    Number_Hub_Genes =
      nrow(hub_genes)
    
  )


export(
  analysis_summary,
  file.path(
    OUT_CSV,
    paste0(
      PREFIX,
      "_ppi_analysis_summary.csv"
    )
  )
)


# ==============================================================================
# 30. FIGURES
# ==============================================================================

g <-
  set_visual_attrs(
    g,
    metrics
  )


message(
  "\nSaving figures..."
)


save_network_plot(
  
  g,
  
  paste(
    tools::toTitleCase(PREFIX),
    "PPI Network -",
    ANALYSIS_MODE
  ),
  
  file.path(
    OUT_FIG,
    paste0(
      PREFIX,
      "_",
      ANALYSIS_MODE,
      "_ppi_network.png"
    )
  )
)


save_ggraph_plot(
  
  g,
  
  metrics,
  
  paste(
    tools::toTitleCase(PREFIX),
    "PPI Network -",
    ANALYSIS_MODE
  ),
  
  file.path(
    OUT_FIG,
    paste0(
      PREFIX,
      "_",
      ANALYSIS_MODE,
      "_ppi_ggraph.png"
    )
  )
)


save_hub_barplot(
  
  hub_genes,
  
  paste(
    tools::toTitleCase(PREFIX),
    "- Hub Genes -",
    ANALYSIS_MODE
  ),
  
  file.path(
    OUT_FIG,
    paste0(
      PREFIX,
      "_",
      ANALYSIS_MODE,
      "_hub_genes_barplot.png"
    )
  )
)


save_module_plots(
  
  g,
  
  OUT_FIG,
  
  paste0(
    PREFIX,
    "_",
    ANALYSIS_MODE
  )
)


# ==============================================================================
# 31. FINAL MESSAGE
# ==============================================================================

message(
  "\n============================================================"
)

message(
  "PPI ANALYSIS COMPLETE"
)

message(
  "Mode: ",
  ANALYSIS_MODE
)

message(
  "Genes used: ",
  nrow(sig)
)

message(
  "Network nodes: ",
  vcount(g)
)

message(
  "Network edges: ",
  ecount(g)
)

message(
  "============================================================\n"
)