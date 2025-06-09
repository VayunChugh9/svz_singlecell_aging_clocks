library(Seurat)
library(tidyverse)
library(here)

# Load subsetted object from step 1
fap_subset <- readRDS(here("data", "fap_subset.rds"))

# Convert counts to dataframe
Convert_to_Dataframe <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if (!is.null(obj[["SCT"]])) obj[["SCT"]] <- NULL
  meta <- obj@meta.data
  meta <- meta[, c("orig.ident", "age", "adjusted_celltype")]
  raw_counts <- t(as.matrix(obj[["RNA"]]@counts))
  raw_counts <- raw_counts[, colSums(raw_counts) > 0]
  df <- as_tibble(cbind(meta, raw_counts))
  return(df)
}

df <- Convert_to_Dataframe(fap_subset) %>%
  group_by(adjusted_celltype, age, orig.ident) %>%
  nest()

# Bootstrapping pseudocells
bootstrap.pseudocells <- function(df, size = 15, n = 100, replace = "dynamic") {
  pseudocells <- c()
  if (replace == "dynamic") {
    if (nrow(df) <= size) { replace <- TRUE } else { replace <- FALSE }
  }
  for (i in seq_len(n)) {
    batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
    pseudocells <- rbind(pseudocells, colSums(batch))
  }
  colnames(pseudocells) <- colnames(df)
  as_tibble(pseudocells)
}

set.seed(42)
df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells))
df2$data <- NULL

saveRDS(df2, here("data", "fap_bootstrap_pseudocell_seed42.rds"))
