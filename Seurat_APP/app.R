library(shiny)
library(Seurat)
library(ggplot2)
library(DT)
library(plotly)
library(enrichR)
library(randomForest)
library(caret)  # for train/test split
library(markdown)


options(shiny.maxRequestSize = 500*1024^1000)  # 500 MB


ui <- fluidPage(
  titlePanel("JCAP Single Cell RNA-seq Solution"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "flame_theme.css")),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("counts_file", "Upload Counts CSV", accept = ".csv"),
      fileInput("meta_file", "Upload Metadata CSV", accept = ".csv"),
      actionButton("create_seurat", "Create Seurat Object"),
      hr(),
      actionButton("run_pca", "Run PCA"),
      actionButton("run_umap", "Run UMAP"),
      actionButton("run_condition_de", "Run DE by Condition"),
      actionButton("run_celltype_de", "Run DE by Cell Type"),
      actionButton("filter_condition_only", "Find Condition-Only DE Genes"),
      hr(),
      selectInput("gene_to_plot", "Choose Gene to Plot", choices = NULL),
      selectInput("enrichr_db", "Choose EnrichR Database",
                  choices = enrichR::listEnrichrDbs()$libraryName),
      actionButton("enrich_all", "Enrich ALL DE Genes"),
      actionButton("enrich_up", "Enrich UPregulated Genes"),
      actionButton("enrich_down", "Enrich DOWNregulated Genes"),
      actionButton("feature_select", "Run Feature Selection (Random Forest)"),
      numericInput("top_n_features", "Top N Features to Use for Classification:", value = 10, min = 1),
      actionButton("classify_rf", "Classify via Random Forest"),
      actionButton("run_power", "Run Power Analysis"),
      hr(),
      h4("Download Tables"),
      downloadButton("download_condition_de", "Condition DE Table"),
      downloadButton("download_celltype_de", "Cell Type DE Table"),
      downloadButton("download_condition_only", "Condition-Only DE Table"),
      downloadButton("download_enrich_all", "EnrichR All Table"),
      downloadButton("download_enrich_up", "EnrichR Upregulated Table"),
      downloadButton("download_enrich_down", "EnrichR Downregulated Table"),
      downloadButton("download_rf_importance", "Random Forest Importance"),
      downloadButton("download_rf_predictions", "RF Predictions"),
      downloadButton("download_rf_metrics", "RF Metrics"),
      downloadButton("download_power_table", "Power Table")
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("How to Use",
                 includeMarkdown("Readme.txt")
        ),
        
        
        tabPanel("PCA Plot", plotlyOutput("pca_plot")),
        tabPanel("UMAP Plot", plotlyOutput("umap_plot")),
        tabPanel("Condition DE Table", DTOutput("condition_de_table")),
        tabPanel("Condition DE Plot", plotOutput("condition_de_plot")),
        tabPanel("Cell Type DE Table", DTOutput("celltype_de_table")),
        tabPanel("Cell Type DE Plot", plotOutput("celltype_de_plot")),
        tabPanel("Condition-Only DE Table", DTOutput("condition_only_table")),
        tabPanel("Condition-Only DE Plot", plotOutput("condition_only_plot")),
        
        # Pathway Analysis tabset
        tabPanel("Pathway Analysis",
                 tabsetPanel(
                   tabPanel("All DE Table", DTOutput("enrich_all_table")),
                   tabPanel("All DE Barplot", plotlyOutput("enrich_all_plot")),
                   tabPanel("Upregulated Table", DTOutput("enrich_up_table")),
                   tabPanel("Upregulated Barplot", plotlyOutput("enrich_up_plot")),
                   tabPanel("Downregulated Table", DTOutput("enrich_down_table")),
                   tabPanel("Downregulated Barplot", plotlyOutput("enrich_down_plot"))
                 )
        ),
        
        tabPanel("Random Forest",
                 tabsetPanel(
                   tabPanel("Variable Importance", DTOutput("rf_table")),
                   tabPanel("Predictions", DTOutput("rf_predictions")),
                   tabPanel("Metrics", DTOutput("rf_metrics")),
                   tabPanel("ROC Curve", plotlyOutput("rf_roc")),
                   tabPanel("Heatmap", plotlyOutput("rf_heatmap"))
                   
                 )
        ),
        tabPanel("Power Analysis",
                 tabsetPanel(
                   tabPanel("Summary Table",
                            DTOutput("power_table"),
                            selectInput("power_gene", "Select Gene for Power Curve:", choices = NULL)
                   ),
                   tabPanel("Power Curve",
                            plotlyOutput("power_curve")
        )
        )
      )
    )
  )
) )

    



server <- function(input, output, session) {
  rv <- reactiveValues(
    seurat = NULL,
    condition_de = NULL,
    celltype_de = NULL,
    condition_only = NULL
  )
  rv$enrich <- list(all = NULL, up = NULL, down = NULL)
  
  rv$rf_importance <- NULL
  rv$classify <- list()
  
  
  render_enrich_barplot <- function(df) {
    df <- df[order(df$Adjusted.P.value), ]
    df <- head(df, 10)
    
    plot_ly(
      df,
      x = ~-log10(Adjusted.P.value),
      y = ~reorder(Term, -Adjusted.P.value),
      type = "bar",
      orientation = "h",
      marker = list(
        color = 'mediumpurple',
        line = list(color = 'purple4', width = 1)
      )
    ) %>%
      layout(
        title = "Top Enriched Pathways",
        xaxis = list(title = "-log10(Adjusted P-value)"),
        yaxis = list(title = "", automargin = TRUE)
      )
  }
  
  
  observeEvent(input$create_seurat, {
    req(input$counts_file, input$meta_file)
    showNotification("📦 Creating Seurat object...", type = "message")
    
    counts <- read.csv(input$counts_file$datapath, row.names = 1)
    meta <- read.csv(input$meta_file$datapath, row.names = 1)
    
    seu <- CreateSeuratObject(counts = counts, meta.data = meta)
    seu <- NormalizeData(seu)
    rv$seurat <- seu
    
    showNotification("✅ Seurat object created and normalized.", type = "message")
  })
  
  observeEvent(input$run_condition_de, {
    req(rv$seurat)
    showNotification("🧪 Running DE by Condition...", type = "message")
    
    rv$seurat$celltype.stim <- paste(rv$seurat$seurat_annotations, rv$seurat$stim, sep = "_")
    Idents(rv$seurat) <- "celltype.stim"
    de <- FindMarkers(rv$seurat, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
    de$gene <- rownames(de)
    rv$condition_de <- de
    updateSelectInput(session, "gene_to_plot", choices = rownames(de))
    showNotification("✅ Finished DE by Condition.", type = "message")
  
  })
  
  observeEvent(input$run_celltype_de, {
    req(rv$seurat)
    showNotification("🧬 Running DE by Cell Type...", type = "message")
    
    Idents(rv$seurat) <- "seurat_annotations"
    de <- FindMarkers(rv$seurat, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono", verbose = FALSE)
    de$gene <- rownames(de)
    rv$celltype_de <- de
    showNotification("✅ Finished DE by Cell Type.", type = "message")
    
  })
  observeEvent(input$filter_condition_only, {
    req(rv$condition_de, rv$celltype_de)
    showNotification("🔍 Filtering Condition-Only DE Genes...", type = "message")
    
    cond_genes <- rownames(rv$condition_de[rv$condition_de$p_val_adj < 0.05, ])
    celltype_genes <- rownames(rv$celltype_de[rv$celltype_de$p_val_adj < 0.05, ])
    only_cond <- setdiff(cond_genes, celltype_genes)
    filtered <- rv$condition_de[rownames(rv$condition_de) %in% only_cond, ]
    rv$condition_only <- filtered
    
    showNotification("✅ Condition-Only DE Genes filtered.", type = "message")
  })
  
  
  # Tables
  output$condition_de_table <- renderDT({
    req(rv$condition_de)
    datatable(rv$condition_de, options = list(pageLength = 10))
  })
  
  output$celltype_de_table <- renderDT({
    req(rv$celltype_de)
    datatable(rv$celltype_de, options = list(pageLength = 10))
  })
  
  output$condition_only_table <- renderDT({
    req(rv$condition_only)
    datatable(rv$condition_only, options = list(pageLength = 10))
  })
  
  # Plots
  output$condition_de_plot <- renderPlot({
    req(input$gene_to_plot, rv$condition_de)
    VlnPlot(rv$seurat, features = input$gene_to_plot, group.by = "stim") + theme_minimal()
  })
  
  output$celltype_de_plot <- renderPlot({
    req(input$gene_to_plot, rv$celltype_de)
    VlnPlot(rv$seurat, features = input$gene_to_plot, group.by = "seurat_annotations") + theme_minimal()
  })
  
  output$condition_only_plot <- renderPlot({
    req(input$gene_to_plot, rv$condition_only)
    VlnPlot(rv$seurat, features = input$gene_to_plot, group.by = "stim") + theme_minimal()
  })

  observeEvent(input$run_umap, {
    req(rv$seurat)
    showNotification("📉 Running UMAP...", type = "message")
    
    if (!"pca" %in% Reductions(rv$seurat)) {
      rv$seurat <- ScaleData(rv$seurat, verbose = FALSE)
      rv$seurat <- RunPCA(rv$seurat, verbose = FALSE)
    }
    
    rv$seurat <- RunUMAP(rv$seurat, dims = 1:10)
    
    showNotification("✅ UMAP completed.", type = "message")
  })
  
  observeEvent(input$run_pca, {
    req(rv$seurat)
    showNotification("🔎 Running PCA...", type = "message")
    
    rv$seurat <- FindVariableFeatures(rv$seurat, selection.method = "vst", nfeatures = 2000)
    rv$seurat <- ScaleData(rv$seurat, verbose = FALSE)
    rv$seurat <- RunPCA(rv$seurat, verbose = FALSE)
    
    showNotification("✅ PCA completed.", type = "message")
  })
  
  output$umap_plot <- renderPlotly({
    req(rv$seurat)
    req("umap" %in% Reductions(rv$seurat))
    p <- DimPlot(rv$seurat, reduction = "umap", group.by = "seurat_annotations", label = TRUE) + theme_minimal()
    ggplotly(p)
  })
  output$pca_plot <- renderPlotly({
    req(rv$seurat)
    req("pca" %in% Reductions(rv$seurat))
    p <- DimPlot(rv$seurat, reduction = "pca", group.by = "seurat_annotations", label = TRUE) + theme_minimal()
    ggplotly(p)
  })
  observeEvent(input$enrich_all, {
    req(rv$condition_only)
    showNotification("🔬 Running pathway enrichment on ALL DE genes...", type = "message")
    
    genes <- rownames(rv$condition_only)
    db <- input$enrichr_db
    result <- enrichr(genes, db)
    rv$enrich$all <- result[[db]]
    
    showNotification("✅ Enrichment complete (ALL DE genes).", type = "message")
  })
  
  observeEvent(input$enrich_up, {
    req(rv$condition_only)
    showNotification("🔬 Running enrichment on UPregulated genes...", type = "message")
    
    de <- rv$condition_only
    up_genes <- rownames(de[de$avg_log2FC > 0, ])
    db <- input$enrichr_db
    result <- enrichr(up_genes, db)
    rv$enrich$up <- result[[db]]
    
    showNotification("✅ Enrichment complete (UPregulated genes).", type = "message")
  })
  
  observeEvent(input$enrich_down, {
    req(rv$condition_only)
    showNotification("🔬 Running enrichment on DOWNregulated genes...", type = "message")
    
    de <- rv$condition_only
    down_genes <- rownames(de[de$avg_log2FC < 0, ])
    db <- input$enrichr_db
    result <- enrichr(down_genes, db)
    rv$enrich$down <- result[[db]]
    
    showNotification("✅ Enrichment complete (DOWNregulated genes).", type = "message")
  })
  
  output$enrich_all_table <- renderDT({
    req(rv$enrich$all)
    datatable(rv$enrich$all)
  })
  
  output$enrich_up_table <- renderDT({
    req(rv$enrich$up)
    datatable(rv$enrich$up)
  })
  
  output$enrich_down_table <- renderDT({
    req(rv$enrich$down)
    datatable(rv$enrich$down)
  })
  
  output$enrich_all_plot <- renderPlotly({
    req(rv$enrich$all)
    render_enrich_barplot(rv$enrich$all)
  })
  
  output$enrich_up_plot <- renderPlotly({
    req(rv$enrich$up)
    render_enrich_barplot(rv$enrich$up)
  })
  
  output$enrich_down_plot <- renderPlotly({
    req(rv$enrich$down)
    render_enrich_barplot(rv$enrich$down)
  })
  observeEvent(input$feature_select, {
    req(rv$seurat, rv$condition_only)
    showNotification("🌲 Running Random Forest on Condition-Only DE Genes...", type = "message")
    
    # Extract expression matrix of DE genes
    features <- rownames(rv$condition_only)
    expr <- t(as.matrix(GetAssayData(rv$seurat, assay = "RNA", slot = "data")[features, ]))  # log-normalized
    
    # Get labels for classification (stim)
    labels <- rv$seurat$stim
    df <- data.frame(expr)
    df$stim <- as.factor(labels)
    
    # Train/test split
    set.seed(42)
    train_index <- createDataPartition(df$stim, p = 0.7, list = FALSE)
    train_data <- df[train_index, ]
    test_data <- df[-train_index, ]
    
    # Fit Random Forest
    rf_model <- randomForest(stim ~ ., data = train_data, importance = TRUE, ntree = 500)
    
    # Get importance
    importance_df <- as.data.frame(importance(rf_model))
    importance_df$Gene <- rownames(importance_df)
    importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
    
    rv$rf_importance <- importance_df
    
    showNotification("✅ Feature selection complete.", type = "message")
  })
  
  output$rf_table <- renderDT({
    req(rv$rf_importance)
    datatable(rv$rf_importance, options = list(pageLength = 15), rownames = FALSE)
  })
  observeEvent(input$classify_rf, {
    req(rv$seurat, rv$rf_importance, input$top_n_features)
    showNotification("🤖 Classifying using Random Forest...", type = "message")
    
    # ✅ Correct: Use top N most important genes from rf_importance
    top_genes <- head(rv$rf_importance$Gene, input$top_n_features)
    
    # Extract log-normalized expression for top genes
    expr <- t(as.matrix(GetAssayData(rv$seurat, assay = "RNA", slot = "data")[top_genes, ]))
    labels <- rv$seurat$stim
    df <- data.frame(expr)
    df$stim <- factor(labels)
    
    # Train/test split
    set.seed(101)
    split <- createDataPartition(df$stim, p = 0.7, list = FALSE)
    train <- df[split, ]
    test <- df[-split, ]
    
    # Train RF classifier
    model <- randomForest(stim ~ ., data = train, ntree = 500)
    
    # Predict
    pred <- predict(model, newdata = test, type = "response")
    probs <- predict(model, newdata = test, type = "prob")[, 2]  # assuming binary
    
    # Metrics
    actual <- test$stim
    conf <- caret::confusionMatrix(pred, actual, positive = "STIM")
    
    # AUC and ROC
    if (length(levels(actual)) == 2) {
      roc_obj <- pROC::roc(actual, probs)
      auc_val <- pROC::auc(roc_obj)
    } else {
      roc_obj <- NULL
      auc_val <- NA
    }
    
    # Save to reactive values
    rv$classify <- list(
      predictions = data.frame(Actual = actual, Predicted = pred, Probability = probs),
      confusion = conf,
      auc = auc_val,
      roc = roc_obj
    )
    
    showNotification("✅ Classification complete!", type = "message")
  })
  
  
  output$rf_predictions <- renderDT({
    req(rv$classify$predictions)
    datatable(rv$classify$predictions)
  })
  
  output$rf_roc <- renderPlotly({
    req(rv$classify$roc)
    
    roc_df <- data.frame(
      FPR = 1 - rv$classify$roc$specificities,
      TPR = rv$classify$roc$sensitivities
    )
    
    plot_ly(roc_df, x = ~FPR, y = ~TPR, type = 'scatter', mode = 'lines',
            line = list(color = 'purple')) %>%
      layout(
        title = "ROC Curve",
        xaxis = list(title = "False Positive Rate"),
        yaxis = list(title = "True Positive Rate"),
        showlegend = FALSE
      )
  })
  output$rf_metrics <- renderDT({
    req(rv$classify$confusion)
    
    # Extract byClass matrix
    metrics_df <- as.data.frame(t(rv$classify$confusion$byClass))
    metrics_df$Metric <- rownames(metrics_df)
    rownames(metrics_df) <- NULL
    
    # Move Metric column to front
    metrics_df <- metrics_df[, c(ncol(metrics_df), 1:(ncol(metrics_df)-1))]
    
    # Add AUC column (if available)
    if (!is.na(rv$classify$auc)) {
      metrics_df$AUC <- NA
      metrics_df$AUC[1] <- rv$classify$auc  # or repeat for all rows if preferred
    }
    
    datatable(metrics_df, options = list(dom = 't', pageLength = nrow(metrics_df)))
  })
  
  observeEvent(input$run_power, {
    req(rv$seurat, rv$condition_only)
    showNotification("🔬 Running Power Analysis on Top 10 DE Genes...", type = "message")
    
    # Use top 10 most significant DE genes (lowest adj p-value)
    de_tab <- rv$condition_only
    top_genes <- head(rownames(de_tab[order(de_tab$p_val_adj), ]), 10)
    labels <- rv$seurat$stim
    features <- GetAssayData(rv$seurat, assay = "RNA", slot = "data")
    
    power_list <- lapply(top_genes, function(gene) {
      expr <- as.numeric(features[gene, ])
      g1 <- expr[labels == unique(labels)[1]]
      g2 <- expr[labels == unique(labels)[2]]
      if(length(g1) < 2 || length(g2) < 2) return(NULL)
      eff_size <- abs(mean(g1) - mean(g2)) / sqrt((var(g1) + var(g2))/2)
      n1 <- length(g1)
      n2 <- length(g2)
      pwr_res <- pwr::pwr.t.test(n = min(n1, n2), d = eff_size, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
      n_seq <- seq(5, max(100, max(n1, n2) * 2), by = 1)
      powers <- sapply(n_seq, function(n)
        pwr::pwr.t.test(n = n, d = eff_size, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
      )
      min_n_08 <- n_seq[which(powers >= 0.8)[1]]
      list(Gene = gene, EffectSize = round(eff_size, 3), ObservedN1 = n1, ObservedN2 = n2,
           ObservedPower = round(pwr_res$power, 3), N_for_0.8 = min_n_08,
           CurveX = n_seq, CurveY = powers)
    })
    power_list <- Filter(Negate(is.null), power_list)
    
    # Make table
    tab <- data.frame(
      Gene = sapply(power_list, `[[`, "Gene"),
      EffectSize = sapply(power_list, `[[`, "EffectSize"),
      ObservedN1 = sapply(power_list, `[[`, "ObservedN1"),
      ObservedN2 = sapply(power_list, `[[`, "ObservedN2"),
      ObservedPower = sapply(power_list, `[[`, "ObservedPower"),
      N_for_0.8 = sapply(power_list, `[[`, "N_for_0.8")
    )
    rv$power_table <- tab
    
    # Save all curves in a named list
    rv$power_curves <- setNames(
      lapply(power_list, function(pl) data.frame(SampleSizePerGroup = pl$CurveX, Power = pl$CurveY)),
      sapply(power_list, `[[`, "Gene")
    )
    
    # Update gene selector
    updateSelectInput(session, "power_gene", choices = tab$Gene, selected = tab$Gene[1])
  })
  
  output$power_table <- renderDT({
    req(rv$power_table)
    datatable(rv$power_table, options = list(dom = 't'))
  })
  
  output$power_curve <- renderPlotly({
    req(rv$power_curves, input$power_gene)
    curve_df <- rv$power_curves[[input$power_gene]]
    plot_ly(curve_df, x = ~SampleSizePerGroup, y = ~Power, type = "scatter", mode = "lines") %>%
      layout(title = paste("Power Curve:", input$power_gene), 
             xaxis = list(title = "Sample Size per Group"),
             yaxis = list(title = "Power"),
             shapes = list(
               list(type = "line", x0 = min(curve_df$SampleSizePerGroup), x1 = max(curve_df$SampleSizePerGroup),
                    y0 = 0.8, y1 = 0.8, line = list(dash = "dash", color = "red"))
             ))
  })
  
  # Download: Condition DE Table
  output$download_condition_de <- downloadHandler(
    filename = function() {"condition_de_table.csv"},
    content = function(file) {
      write.csv(rv$condition_de, file)
    }
  )
  
  # Download: Cell Type DE Table
  output$download_celltype_de <- downloadHandler(
    filename = function() {"celltype_de_table.csv"},
    content = function(file) {
      write.csv(rv$celltype_de, file)
    }
  )
  
  # Download: Condition-Only DE Table
  output$download_condition_only <- downloadHandler(
    filename = function() {"condition_only_de_table.csv"},
    content = function(file) {
      write.csv(rv$condition_only, file)
    }
  )
  
  # Download: EnrichR All Table
  output$download_enrich_all <- downloadHandler(
    filename = function() {"enrichr_all_table.csv"},
    content = function(file) {
      write.csv(rv$enrich$all, file)
    }
  )
  
  # Download: EnrichR Upregulated Table
  output$download_enrich_up <- downloadHandler(
    filename = function() {"enrichr_upregulated_table.csv"},
    content = function(file) {
      write.csv(rv$enrich$up, file)
    }
  )
  
  # Download: EnrichR Downregulated Table
  output$download_enrich_down <- downloadHandler(
    filename = function() {"enrichr_downregulated_table.csv"},
    content = function(file) {
      write.csv(rv$enrich$down, file)
    }
  )
  
  # Download: Random Forest Importance
  output$download_rf_importance <- downloadHandler(
    filename = function() {"rf_importance_table.csv"},
    content = function(file) {
      write.csv(rv$rf_importance, file)
    }
  )
  
  # Download: Random Forest Predictions
  output$download_rf_predictions <- downloadHandler(
    filename = function() {"rf_predictions_table.csv"},
    content = function(file) {
      write.csv(rv$classify$predictions, file)
    }
  )
  
  # Download: RF Metrics
  output$download_rf_metrics <- downloadHandler(
    filename = function() {"rf_metrics_table.csv"},
    content = function(file) {
      # Use as.data.frame for the confusion/metrics, as rendered in the app
      if (!is.null(rv$classify$confusion)) {
        metrics_df <- as.data.frame(t(rv$classify$confusion$byClass))
        metrics_df$Metric <- rownames(metrics_df)
        rownames(metrics_df) <- NULL
        metrics_df <- metrics_df[, c(ncol(metrics_df), 1:(ncol(metrics_df)-1))]
        if (!is.na(rv$classify$auc)) {
          metrics_df$AUC <- NA
          metrics_df$AUC[1] <- rv$classify$auc
        }
        write.csv(metrics_df, file)
      } else {
        write.csv(data.frame(), file)
      }
    }
  )
  
  # Download: Power Table
  output$download_power_table <- downloadHandler(
    filename = function() {"power_table.csv"},
    content = function(file) {
      write.csv(rv$power_table, file)
    }
  )
  
  output$rf_heatmap <- renderPlotly({
    req(rv$seurat, rv$rf_importance, input$top_n_features)
    # Get top N important genes used for classification
    top_genes <- head(rv$rf_importance$Gene, input$top_n_features)
    
    # Extract expression data
    expr_mat <- as.matrix(GetAssayData(rv$seurat, assay = "RNA", slot = "data")[top_genes, ])
    # Optionally scale genes for visualization
    expr_mat_scaled <- t(scale(t(expr_mat)))
    
    # Sample annotations (e.g. by stim)
    annots <- rv$seurat$stim
    # Assign a color for each group
    annot_colors <- setNames(RColorBrewer::brewer.pal(length(unique(annots)), "Set1"), unique(annots))
    sample_colors <- annot_colors[as.character(annots)]
    
    # For Plotly, row/col names
    heat <- plot_ly(
      z = expr_mat_scaled,
      x = colnames(expr_mat_scaled),
      y = top_genes,
      type = "heatmap",
      colorscale = "Viridis",
      showscale = TRUE
    ) %>%
      layout(
        title = paste("Heatmap of Top", input$top_n_features, "Variable Importance Genes"),
        yaxis = list(title = "Gene", autorange = "reversed"),
        xaxis = list(title = "Sample", tickangle = 45)
      )
    
    heat
  })
  
  
}

shinyApp(ui, server)
