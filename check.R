library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

meta <- import("results/tables/meta-analysis/Parkinson_random_effects_meta_analysis.csv")

ensembl_ids <- unique(
  meta$Gene_ID[!is.na(meta$Gene_ID)]
)

mapping <- bitr(
  ensembl_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

unmapped <- setdiff(
  ensembl_ids,
  mapping$ENSEMBL
)

length(ensembl_ids)
length(mapping$ENSEMBL)
length(unmapped)





test_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unmapped,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME")
)

test_mapping
