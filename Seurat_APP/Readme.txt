Welcome to JCAP Single Cell RNA-seq solution!
This app allows you to:

Explore your single-cell RNA-seq data interactively

Run differential expression, pathway analysis, classification, and power analysis—all in your browser

Getting Started
Prepare your files:

Counts matrix:

CSV file

Rows = genes/features, Columns = cell/sample names

First column = gene names, first row = sample names

Metadata:

CSV file

Rows = cell/sample names (must match columns in counts)

Columns = sample annotations (e.g. group, batch, cell type, etc.)

First column = sample names

Upload your files:

Click "Upload Counts CSV" and select your counts matrix file

Click "Upload Metadata CSV" and select your metadata file

Create Seurat object:

Click "Create Seurat Object" to load and normalize your data

Analysis Workflow
Run PCA or UMAP:
Visualize your data in reduced dimensions

Differential Expression:

Run DE by condition and ,then cell type

Find condition-only DE genes (not confounded by cell type)

Pathway Analysis:

Enrich differentially expressed gene sets (choose database)

Feature Selection & Classification:

Run Random Forest to identify top marker genes

Classify samples/cells using selected genes

Power Analysis:

Estimate statistical power and minimum sample size for key genes

Visualization:

Generate plots for DE, classification, and power analysis

Create heatmaps of top features

Download Results:

Download tables for all analyses using the buttons at the bottom of the sidebar

Tips & Requirements
Data files must be CSV format and must be pre-filtered for quality (consider using "insurance_policy_script" provided in the repository)

Columns and rows in counts/metadata files must match exactly by sample/cell name

Works best with well-annotated, normalized, and quality-checked datasets

Need Help?
Check the repository README for detailed file format examples and troubleshooting

For support or feature requests, contact the app maintainer or open an issue

Happy analyzing!