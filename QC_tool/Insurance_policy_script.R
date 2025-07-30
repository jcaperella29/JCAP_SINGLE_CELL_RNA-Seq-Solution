
# insure_policy.R
library(Seurat)
library(dplyr)

# Load raw data
counts <- read.csv("C:/Users/jcape/Downloads/Seurat_data/ifnb_counts.csv", row.names=1, check.names=FALSE)
meta <- read.csv("C:/Users/jcape/Downloads/Seurat_data/ifnb_metadata.csv", row.names=1, check.names=FALSE)

# Create temp Seurat object for QC/filtering
seu <- CreateSeuratObject(counts=counts, meta.data=meta)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# QC filter (customize as needed)
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 10)

# --- Smart detection of group column ---
meta_cols <- colnames(seu@meta.data)
group_keywords <- c("group", "stim", "condition", "treat", "status", "class", "label", "cluster")
possible_groups <- meta_cols[grepl(paste(group_keywords, collapse="|"), meta_cols, ignore.case = TRUE)]

cat("Possible group columns in your metadata:\n")
print(possible_groups)

if (length(possible_groups) == 0) {
  stop("No obvious group column found! Please edit the script to specify your grouping variable for plotting.")
} else if (length(possible_groups) == 1) {
  selected_group <- possible_groups[1]
  cat("Using column:", selected_group, "\n")
} else {
  # If multiple possible, pick the first, or use interactive select if running interactively
  if (interactive()) {
    selected_group <- select.list(possible_groups, title = "Select a grouping variable for coloring:")
    cat("Using column:", selected_group, "\n")
  } else {
    selected_group <- possible_groups[1]
    cat("Multiple possible group columns. Using first:", selected_group, "\n")
  }
}

# Metadata sanity
stopifnot(length(unique(seu@meta.data[[selected_group]])) > 1)
stopifnot(all(!is.na(seu@meta.data[[selected_group]])))
print(table(seu@meta.data[[selected_group]]))

# Optional: PCA/UMAP check (color by selected group column)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
png("PCA_group.png", width=600, height=500)
print(DimPlot(seu, reduction="pca", group.by=selected_group) + ggtitle(paste("PCA: Colored by", selected_group)))
dev.off()

# Save filtered counts and metadata
cleaned_counts <- GetAssayData(seu, slot = "counts")
write.csv(as.data.frame(cleaned_counts), "cleaned_counts.csv")
write.csv(seu@meta.data, "cleaned_meta.csv")

cat("Insure policy checks complete! Cleaned counts and metadata saved.\n")
