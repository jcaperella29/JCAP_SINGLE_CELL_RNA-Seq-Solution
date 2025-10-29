🧬 JCAP Single-Cell RNA-seq Solution

Welcome to JCAP Single-Cell RNA-seq Solution!

This interactive Shiny app lets you explore, analyze, and visualize single-cell RNA-seq datasets entirely in your browser — now with multi-species support (Human, Mouse, Zebrafish, and Fly).

🚀 What You Can Do

Explore your single-cell data interactively (PCA, UMAP)

Run differential expression (DE) by condition and cell type

Find condition-only DE genes (those not confounded by cell type)

Perform pathway enrichment across supported species using g:Profiler

Select top marker genes via Random Forest feature selection

Visualize key genes with the interactive heatmap

Estimate power and required sample size

Download all results in convenient tables and figures

🧩 Supported Species

You can seamlessly analyze data from:

Species	Organism Code	Example IDs
Human	hsapiens	ENSG00000… / Gene symbols
Mouse	mmusculus	ENSMUSG00000… / Gene symbols
Zebrafish	drerio	ENSDARG00000…
Fly	dmelanogaster	FBgn000000…

The app automatically detects whether your dataset uses Ensembl IDs or Gene Symbols, and converts identifiers as needed for enrichment and visualization.

🧱 Getting Started
1️⃣ Prepare Your Files

Counts Matrix (CSV)

Rows = genes/features

Columns = cell/sample names

First column = gene names

First row = sample names

Metadata (CSV)

Rows = cell/sample names (must match count matrix columns)

Columns = annotations (stim, cell_type, batch, etc.)

Must include:

stim (experimental condition)

cell_type (cell classification)

2️⃣ Upload Your Files

Click Upload Counts CSV → select your count matrix

Click Upload Metadata CSV → select your metadata

3️⃣ Create the Seurat Object

Click “Create Seurat Object”

Your data is normalized

stim and cell_type are automatically converted to factors

Metadata is validated for DE analysis

🔬 Analysis Workflow
🧠 PCA / UMAP

Visualize global structure of your dataset

Works for all species automatically

⚖️ Differential Expression

Run DE by Condition → contrasts by stim

Run DE by Cell Type → contrasts by cell_type

Find Condition-Only DE Genes → removes overlaps to isolate genes driven purely by condition

🧬 Pathway Enrichment

Click “Enrich (UP)” or “Enrich (DOWN)” to analyze up/down-regulated DE sets

Automatically maps gene IDs for your selected species

Uses g:Profiler2 with GO, Reactome, KEGG, and MSigDB

🌲 Feature Selection & Classification

Run Random Forest to identify top features by MeanDecreaseGini

Table displays ranked genes and importance values

🔥 Heatmap

Click “Run Heatmap 🔥”

Displays expression for top 10 RF genes, prioritizing those also condition-specific

Works across all supported species with automatic ID translation

Interactive via Plotly (rf_heatmap output)

📊 Power Analysis

Estimate statistical power and required sample size for key DE genes

💾 Download Results

Each major analysis block has a Download button to export:

DE tables (condition, celltype, and condition-only)

Enrichment tables (UP, DOWN, ALL)

Random Forest importance table

Power analysis results

🧠 Tips & Requirements

Files must be CSVs and match exactly by sample/cell name

Data should be pre-filtered and normalized

For best results, use the included insurance_policy_script to QC your input

Works with both Ensembl IDs and Gene Symbols automatically

🛠️ Tech Overview

Built with R + Shiny + Seurat + gprofiler2 + plotly

Modular server observers:

create_seurat — normalization & metadata handling

run_condition_de, run_celltype_de, filter_condition_only — differential expression

enrich_up, enrich_down, enrich_all — pathway enrichment

run_rf — feature selection

run_heatmap — multi-species top-gene visualization

❓Need Help?

See examples and test datasets in the repository

Check logs for notifications (messages appear in the bottom-right corner)

For questions or feature requests, open an issue or contact the maintainer

Happy analyzing!
— JCAP Bioinformatics Team 🧬
