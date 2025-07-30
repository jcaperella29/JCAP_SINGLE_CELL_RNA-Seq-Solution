# download_subset_ifnb.R
# Script to download and subset the Kang et al. IFNB Seurat dataset for tutorials/demos
# Requires: Seurat

library(Seurat)
library(dplyr)

# Download IFNB dataset from SeuratData (hosted by Satija Lab)
if (!requireNamespace("SeuratData", quietly = TRUE)) {
    install.packages("SeuratData")
}
library(SeuratData)
SeuratData::InstallData("ifnb")

# Load full dataset
data("ifnb")
seu <- ifnb

# Optionally, normalize if not already
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

# Subset: grab a small random subset of cells (e.g. 200 cells, balanced by group)
meta <- seu@meta.data
table(meta$stim)  # Should be "CTRL" and "STIM"

set.seed(123)
n_per_group <- 100
cells_ctrl <- sample(rownames(meta)[meta$stim == "CTRL"], n_per_group)
cells_stim <- sample(rownames(meta)[meta$stim == "STIM"], n_per_group)
cells_keep <- c(cells_ctrl, cells_stim)

seu_small <- subset(seu, cells = cells_keep)

# Export counts matrix and metadata as CSV
counts_mat <- GetAssayData(seu_small, slot = "counts")
write.csv(as.data.frame(counts_mat), "ifnb_counts_subset.csv")
write.csv(seu_small@meta.data, "ifnb_metadata_subset.csv")

cat("Subsetted IFNB dataset (", length(cells_keep), "cells) written to CSV files.\n")
