Welcome to JCAP Single Cell RNA-seq solution!
This app allows you to:
•	Explore your single-cell RNA-seq data interactively
•	Run differential expression, pathway analysis, classification, and power analysis—all in your browser

Getting Started
Prepare your files:
Counts matrix:
CSV file
Rows =count, Columns = cells/samples

Metadata:
CSV file
Rows = cell/sample names (must match columns in counts)
Columns = sample annotations (e.g., group, batch, cell type, etc.)

Upload your files:

Click the browse button labeled as "Upload Counts CSV" and select your counts matrix file.

Click the browse button labeled as "Upload Metadata CSV" and select your metadata file.

Create Seurat object:

Click "Create Seurat Object" to load and normalize your data


Analysis Workflow

Step 1:  Exploratory  Visualization
Click Run PCA or  Run UMAP. This will allow you to visualize your data in reduced dimensions. The plots will be shown in the tabs labeled as “PCA Plot” and “UMAP Plot.” 

 Step 2: Differential Expression
Click  “Run DE by Condition” then  “Run DE by Cell Type” and  then click “Find  Condition-only DE Genes.” The last set of hits  are the differentially expressed genes that are not confounded by cell type. Tabulated results from each DE analysis will be  shown in the tabs  labeled “Condition DE Table,” “Cell Type DE table” and “Condition -Only De Table.” Lastly,  you can click “Create Volcano Plot”  to create a volcano plot  using the genes in the Condition-Only DE  Table. Note  that the “Condition-Only De  Table” only takes genes that had adjusted p values of .05 or less by condition  and not p values of .05 or less by cell type.
Step 3 : Pathway Analysis
After performing DE analysis, you can  map the Condition only DE genes to pathways. To do this, you first need to use the “Choose EnrichR Database” drop down menu to pick the database to use when mapping your genes. Next, click “Enrich ALL DE Genes,” “Enrich Upregulated Genes” and “Enrich Downregulated Genes.”  Results will be shown as tables and barplots in the “Pathway Analysis” tab set. The “All DE Table” and “All DE Barplot” show the results of pathway analysis on all of the condition-only  genes and the other tabs in that tab set show the results of pathways analysis for upregulated and downregulated genes and are labeled as “Upregulated Table,” “Unregulated Barplot,” “Downregulated Table” and “Downregulated Barplot”  respectively.

Step 4: Feature Selection & Classification
After  performing DE analysis or optionally performing pathway analysis you can further narrow down your hits using   feature selection via random forest. To do so, click “Run Feature  Selection (Random Forest).” This will create a table ranking the DE genes by variable importance. The table will appear in the tab labeled “Variable Importance.” This is the left most tab in the “Random Forest” tab set. Next you can use the numeric input labeled as “Top N Features to Use for Classification” to determine how many genes you will use as predictors in the classification step. Next you can click “Classify via Random Forest” and  the data will be classified via a  Random Forest model and the following objects will be generated - a table showing the model’s predictions, a table showing the model’s performance metrics, a plot of the ROC curve of the model, and a heatmap  showing the correlation between the selected genes (the ones used as predictors) and the condition.
Step 5 :Power Analysis
Finally, you can click “Run Power Analysis” to perform Power Analysis  on your dataset. This will generate a table showing the statistical power of the dataset for each of the topmost differentially expressed genes (Shown in the  Power Analysis tab set in “Summary Table”). In the same tab you will see a drop down labeled as ‘Select Gene for Power Curve.” As the label suggests, this lets you pick which gene to focus on for the power curve visualization  that will be displayed in the adjacent tab which is labeled as “Power Curve.” Note that the  dashed red line on the plot shows how many samples would be needed to achieve  80% power for  the selected gene and the actual curve is shown in blue.
Tips & Requirements
Data files must be CSV format and must be pre-filtered for quality (consider using "insurance_policy_script" provided in the repository).
Columns and rows in the counts/metadata files must match exactly by sample/cell name.
Works best with well-annotated, normalized, and quality-checked datasets.

Need Help?
Check the repository README for detailed file format examples and troubleshooting.
For support or feature requests, contact the app maintainer or open an issue on GitHub.
Happy analyzing!
