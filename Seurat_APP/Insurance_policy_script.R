# insure_policy.R
library(Seurat)
library(ggplot2)
library(dplyr)

# Load raw data
counts <- read.csv("counts.csv", row.names=1)
meta <- read.csv("meta.csv", row.names=1)

# 1. Basic QC filter
seu <- CreateSeuratObject(counts=counts, meta.data=meta)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 10)

# 2. Metadata sanity
stopifnot(length(unique(seu$group)) > 1)
stopifnot(all(!is.na(seu$group)))
print(table(seu$group, seu$batch))

# 3. PCA/UMAP for pre-analysis
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
DimPlot(seu, group.by = "group") + ggtitle("Group Color")
DimPlot(seu, group.by = "batch") + ggtitle("Batch Color")

# 4. Flag if batch == group
if (any(table(seu$group, seu$batch) == ncol(seu))) {
  warning("Batch and group labels are perfectly aligned!")
}

# 5. Save cleaned data
saveRDS(seu, "cleaned_seurat.rds")
write.csv(seu@meta.data, "cleaned_meta.csv")

cat("Insure policy checks complete!\n")
