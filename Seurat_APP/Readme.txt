# Welcome to the JCAP Single Cell RNA-seq Solution!

This interactive app lets you:
- Explore your single-cell RNA-seq data visually and statistically
- Run differential expression, pathway enrichment, classification, and power analyses—all in your browser
- Download publication-ready tables and interactive plots for every step

---

## **Getting Started**

### **1. Prepare Your Data**
- **Counts Matrix** (CSV)
  - **Rows:** Genes/features
  - **Columns:** Cell/sample names
  - **First column:** Gene names
  - **First row:** Sample names
- **Metadata** (CSV)
  - **Rows:** Cell/sample names (must match counts matrix columns)
  - **Columns:** Sample annotations (e.g. group, batch, cell type)
  - **First column:** Sample names

### **2. Upload Your Files**
- Click **"Upload Counts CSV"** and select your counts matrix file
- Click **"Upload Metadata CSV"** and select your metadata file

### **3. Create Your Seurat Object**
- Click **"Create Seurat Object"** to load and normalize your data

---

## **Analysis Workflow**

- **PCA & UMAP:**  
  Visualize clusters and patterns in reduced dimensions

- **Differential Expression (DE):**  
  - Run DE by **condition** and by **cell type**
  - Find "condition-only" DE genes (not confounded by cell type)

- **Pathway Analysis:**  
  - Enrich DE gene sets (choose from multiple databases)
  - Results as interactive tables and barplots

- **Volcano Plot:**  
  - Instantly visualize condition-only DE results with a publication-style volcano plot

- **Feature Selection & Classification:**  
  - Identify top marker genes with Random Forest
  - Classify samples/cells and view performance metrics, ROC curves, and variable importance

- **Power Analysis:**  
  - Estimate statistical power and minimum sample size for your key DE genes

- **Heatmaps:**  
  - Visualize top features and group differences

---

## **Download Results**
- **Download all tables and results** directly using the buttons in the sidebar

---

## **Tips & Requirements**

- **Files must be CSV format**
- **Counts and metadata must match by sample/cell name**
- For best results, use **quality-filtered data** (see the included `insurance_policy_script.R` for QC)
- App works best with well-annotated, normalized data

---

## **Need Help?**
- See the repository README for detailed file examples and troubleshooting tips
- For questions or feature requests, contact the maintainer or open a GitHub issue

---

**Happy analyzing!**

