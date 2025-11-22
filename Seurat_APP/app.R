# Multi-species upgrade of your scRNA-seq Shiny app
# - Adds Species selector (Human/Mouse/Fly/Zebrafish)
# - Adds Enrichment backend selector (g:Profiler for cross-species, Enrichr optional)
# - Adds g:Profiler source picker (GO/KEGG/Reactome/MSigDB etc.)
# - Enrichment handlers now route to the selected backend
# - Plot/table logic made backend-agnostic

library(shiny)
library(Seurat)
library(ggplot2)
library(DT)
library(plotly)
library(enrichR)       # optional, works best for human/mouse libs
library(gprofiler2)    # multi-species enrichment
library(randomForest)
library(caret)         # train/test split
library(markdown)
library(pROC)
library(RColorBrewer)

options(shiny.maxRequestSize = 500*1024^1000)  # 500 MB

species_choices <- c(
  "Human" = "hsapiens",
  "Mouse" = "mmusculus",
  "Fly (D. melanogaster)" = "dmelanogaster",
  "Zebrafish (D. rerio)" = "drerio"
)

# g:Profiler sources: keep a sensible default set that exists for most orgs
# (User can toggle more/less; g:Profiler will skip unavailable ones)
gprof_default_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")

ui <- fluidPage(
  titlePanel("JCAP Single Cell RNA-seq Solution (Multi-species)"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "flame_theme.css")),
  
  sidebarLayout(
    sidebarPanel(
      h4("Data"),
      selectInput("species", "Species", choices = species_choices, selected = "hsapiens"),
      selectInput("id_type", "Gene ID type in your matrices (for display only)",
                  choices = c("Symbol", "ENSEMBL", "Entrez"), selected = "Symbol"),
      fileInput("counts_file", "Upload Counts CSV", accept = ".csv"),
      fileInput("meta_file", "Upload Metadata CSV", accept = ".csv"),
      actionButton("create_seurat", "Create Seurat Object"),
      hr(),
      
      h4("Dimensionality Reduction"),
      actionButton("run_pca", "Run PCA"),
      actionButton("run_umap", "Run UMAP"),
      hr(),
      
      h4("Differential Expression"),
      actionButton("run_condition_de", "Run DE by Condition"),
      actionButton("run_celltype_de", "Run DE by Cell Type"),
      actionButton("filter_condition_only", "Find Condition-Only DE Genes"),
      actionButton("make_volcano", "Create Volcano Plot"),
      hr(),
      
      h4("Pathway Enrichment"),
      radioButtons("enrich_backend", "Backend",
                   choices = c("g:Profiler (multi-species)" = "gprof", "Enrichr" = "enrichr"),
                   selected = "gprof"),
      conditionalPanel(
        condition = "input.enrich_backend == 'gprof'",
        checkboxGroupInput("gprof_sources", "g:Profiler Sources",
                           choices = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "MSigDB H", "MSigDB C2"),
                           selected = gprof_default_sources),
        helpText("Tips: GO/Reactome/KEGG are broadly supported across species; MSigDB sets are human-centric and may ortholog-map internally.")
      ),
      conditionalPanel(
        condition = "input.enrich_backend == 'enrichr'",
        selectInput("enrichr_db", "Choose Enrichr Database",
                    choices = enrichR::listEnrichrDbs()$libraryName)
      ),
      actionButton("enrich_all", "Enrich ALL DE Genes"),
      actionButton("enrich_up", "Enrich UPregulated Genes"),
      actionButton("enrich_down", "Enrich DOWNregulated Genes"),
      hr(),
      
      h4("ML: Random Forest"),
      actionButton("feature_select", "Run Feature Selection (Random Forest)"),
      numericInput("top_n_features", "Top N Features to Use for Classification:", value = 10, min = 1),
      actionButton("classify_rf", "Classify via Random Forest"),
      actionButton("run_heatmap", "Run Heatmap 🔥", class = "btn-warning"),
      
      hr(),
      
      h4("Power Analysis"),
      actionButton("run_power", "Run Power Analysis"),
      hr(),
      
      h4("Download Tables"),
      downloadButton("download_condition_de", "Condition DE Table"),
      downloadButton("download_celltype_de", "Cell Type DE Table"),
      downloadButton("download_condition_only", "Condition-Only DE Table"),
      downloadButton("download_enrich_all", "Enrichment: All"),
      downloadButton("download_enrich_up", "Enrichment: Up"),
      downloadButton("download_enrich_down", "Enrichment: Down"),
      downloadButton("download_enrich_all_edges",  "Enrichment: All Term–Gene Edges"),   # 🆕
     downloadButton("download_enrich_up_edges",   "Enrichment: Up Term–Gene Edges"),    # 🆕
      downloadButton("download_enrich_down_edges", "Enrichment: Down Term–Gene Edges"),  # 🆕
      downloadButton("download_rf_importance", "Random Forest Importance"),
      downloadButton("download_rf_predictions", "RF Predictions"),
      downloadButton("download_rf_metrics", "RF Metrics"),
      downloadButton("download_power_table", "Power Table")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("How to Use", includeMarkdown("Readme.txt")),
        tabPanel("PCA Plot", plotlyOutput("pca_plot")),
        tabPanel("UMAP Plot", plotlyOutput("umap_plot")),
        tabPanel("Condition DE Table", DTOutput("condition_de_table")),
        tabPanel("Cell Type DE Table", DTOutput("celltype_de_table")),
        tabPanel("Condition-Only DE Table", DTOutput("condition_only_table")),
        tabPanel("Condition-Only Volcano Plot", plotlyOutput("condition_only_volcano")),
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
                   tabPanel("Power Curve", plotlyOutput("power_curve"))
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
rv <- reactiveValues(
  seurat = NULL,
  condition_de = NULL,
  celltype_de = NULL,
  condition_only = NULL,
  enrich = list(all = NULL, up = NULL, down = NULL),
  enrich_edges = list(all = NULL, up = NULL, down = NULL),  # 🆕 term–gene edges
  rf_importance = NULL,
  classify = list(),
  volcano = NULL,
  power_table = NULL,
  power_curves = NULL
)

  # ---------- Helpers ----------

  # ---- unify to one table shape for BOTH backends ----------------------------
fmt_enrich_for_table <- function(res) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  # g:Profiler result
  if (all(c("term_name","term_size","intersection_size","p_value") %in% colnames(res))) {
    adj_col <- intersect(colnames(res), c("p_adjusted","adjusted_p_value","p_adj"))
    adj <- if (length(adj_col)) res[[adj_col[1]]] else p.adjust(res$p_value, method = "fdr")
    
    # g:Profiler typically has an 'intersection' column with the genes in the term
    genes_col <- if ("intersection" %in% colnames(res)) res$intersection else NA_character_
    
    df <- data.frame(
      Term = paste0(res$term_name, " (", res$source, ")"),
      Overlap = paste0(res$intersection_size, "/", res$term_size),
      P.value = as.numeric(res$p_value),
      Adjusted.P.value = as.numeric(adj),
      Genes = genes_col,                                      # 🆕 keep genes per term
      stringsAsFactors = FALSE, check.names = FALSE
    )
    df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
    rownames(df) <- NULL
    return(df)
  }
  
  # Enrichr result (already has Term/Overlap/P.value/Adjusted.P.value, often Genes)
  if (all(c("Term","Overlap","P.value","Adjusted.P.value") %in% colnames(res))) {
    # If a 'source' column exists, tag the term with it
    if ("source" %in% colnames(res)) {
      res$Term <- paste0(res$Term, " (", res$source, ")")
    } else {
      res$Term <- paste0(res$Term, " (ENRICHR)")
    }
    keep_cols <- intersect(c("Term","Overlap","P.value","Adjusted.P.value","Genes"), colnames(res))
    df <- res[, keep_cols, drop = FALSE]
    df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
    rownames(df) <- NULL
    return(df)
  }
  
  # If Enrichr rbind lost the 'source' attribute but has its own 'source' col
  if ("source" %in% colnames(res) && all(c("Term","Overlap","P.value","Adjusted.P.value") %in% colnames(res))) {
    res$Term <- paste0(res$Term, " (", res$source, ")")
    keep_cols <- intersect(c("Term","Overlap","P.value","Adjusted.P.value","Genes"), colnames(res))
    df <- res[, keep_cols, drop = FALSE]
    df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
    rownames(df) <- NULL
    return(df)
  }
  
  NULL
}
make_edges <- function(enrich_tbl) {
  if (is.null(enrich_tbl) || !"Genes" %in% colnames(enrich_tbl)) return(NULL)
  
  tmp <- enrich_tbl[!is.na(enrich_tbl$Genes) & nchar(enrich_tbl$Genes) > 0, c("Term","Genes"), drop = FALSE]
  if (nrow(tmp) == 0) return(NULL)
  
  edges_list <- lapply(seq_len(nrow(tmp)), function(i) {
    genes <- unlist(strsplit(as.character(tmp$Genes[i]), "[,;]"))
    data.frame(
      Term = tmp$Term[i],
      Gene = trimws(genes),
      stringsAsFactors = FALSE
    )
  })
  
  edges <- do.call(rbind, edges_list)
  edges <- edges[nchar(edges$Gene) > 0, , drop = FALSE]
  rownames(edges) <- NULL
  edges
}

  render_enrich_barplot <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- fmt_enrich_for_table(df)
    if (!all(c("Term", "Adjusted.P.value") %in% names(df))) return(NULL)
    df <- df[order(df$Adjusted.P.value), ]
    df <- head(df, 10)
    
    plot_ly(
      df,
      x = ~-log10(Adjusted.P.value + 1e-300),
      y = ~reorder(Term, -Adjusted.P.value),
      type = "bar",
      orientation = "h"
    ) %>%
      layout(
        title = "Top Enriched Pathways",
        xaxis = list(title = "-log10(FDR)"),
        yaxis = list(title = "", automargin = TRUE)
      )
  }
  
  # Map UI species → g:Profiler organism codes
  gp_org <- function(species_input) {
    switch(tolower(species_input),
           "human"      = "hsapiens",
           "mouse"      = "mmusculus",
           "zebrafish"  = "drerio",
           "fly"        = "dmelanogaster",
           "rat"        = "rnorvegicus",
           "yeast"      = "scerevisiae",
           "mmusculus"  = "mmusculus",    # allow raw codes too
           "drerio"     = "drerio",
           "hsapiens"   = "hsapiens",
           NULL
    )
  }
  
  # Remove synthetic background genes from our sims
  strip_bg <- function(x) x[!grepl("^BG\\d+$", x)]
  

  get_stim_col <- function(seu) {
    cn <- tolower(colnames(seu@meta.data))
    keys <- c("stim","condition","group","treatment")
    hit <- keys[keys %in% cn]
    if (length(hit) == 0) return(NULL)
    colnames(seu@meta.data)[match(hit[1], cn)]
  }
  
  normalize_stim <- function(x) {
    raw <- trimws(as.character(x))
    xl  <- tolower(raw)
    
    # map common synonyms
    ctrl_pat <- c("ctrl","control","vehicle","untreated","baseline","mock","sham","uninfected","con")
    stim_pat <- c("stim","treated","lps","ifn","ifng","ifna","tnf","il6","infection","infected","stimulated")
    
    y <- raw
    y[xl %in% ctrl_pat | raw %in% c("0","FALSE","False","false")] <- "CTRL"
    y[xl %in% stim_pat | raw %in% c("1","TRUE","True","true")]    <- "STIM"
    
    # if still one unique label but raw has ≥2 unique values, split by top-2 values
    if (length(unique(y)) < 2 && length(unique(raw)) >= 2) {
      top2 <- names(sort(table(raw), decreasing = TRUE))[1:2]
      y <- ifelse(raw == top2[1], paste0("L_", top2[1]), paste0("L_", top2[2]))
    }
    
    factor(y)
  }
  
  # optional: quick diagnostics you can keep or remove
  stim_summary <- function(seu, stim_col, ct_col = NULL) {
    ss <- capture.output({
      cat(sprintf("stim column: %s\n", stim_col))
      print(table(seu@meta.data[[stim_col]]))
      if (!is.null(ct_col) && ct_col %in% colnames(seu@meta.data)) {
        cat("\nby cell_type (top 5):\n")
        print(head(table(seu@meta.data[[ct_col]], seu@meta.data[[stim_col]]), 5))
      }
    })
    paste(ss, collapse = "\n")
  }
  
  
  # ---- species code, BG filter ----------------------------------------------
  gp_org <- function(species_input) {
    switch(tolower(species_input),
           "human"="hsapiens","mouse"="mmusculus","zebrafish"="drerio","fly"="dmelanogaster",
           "rat"="rnorvegicus","yeast"="scerevisiae","mmusculus"="mmusculus",
           "drerio"="drerio","hsapiens"="hsapiens", NULL)
  }
  strip_bg <- function(x) x[!grepl("^BG\\d+$", x)]
  
  # ---- g:Profiler with retries (returns data.frame or NULL) ------------------
  safe_gprof <- function(genes, org, retries = 3, sleep_s = c(0.6, 1.2, 2.4)) {
    for (i in seq_len(retries)) {
      try({
        out <- gprofiler2::gost(
          query = unique(genes),
          organism = org,
          sources = c("GO:BP","GO:MF","GO:CC","REAC","KEGG","MSigDBH","MSigDBC2"),
          correction_method = "fdr"
        )
        if (!is.null(out) && !is.null(out$result) && nrow(out$result) > 0) {
          return(out$result)
        }
      }, silent = TRUE)
      if (i < retries) Sys.sleep(sleep_s[i])
    }
    return(NULL)
  }
  
  # ---- Enrichr fallback; choose libs per species -----------------------------
  fallback_enrichr_libs <- function(org) {
    base <- c("GO_Biological_Process_2021","GO_Molecular_Function_2021",
              "GO_Cellular_Component_2021","Reactome_2022")
    if (org %in% c("mmusculus")) c(base, "KEGG_2021_Mouse") else base
  }
  
  run_enrichr_fallback <- function(genes, org) {
    libs <- fallback_enrichr_libs(org)
    res_list <- enrichR::enrichr(unique(genes), libs)
    # Bind and standardize
    out <- do.call(rbind, lapply(names(res_list), function(k) {
      df <- res_list[[k]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      # Enrichr standard columns:
      # Term, Overlap, P.value, Adjusted.P.value, Old.P.value, Odds.Ratio, Combined.Score
      df$source <- k
      df
    }))
    out
  }
  
  # ---- unify to one table shape for BOTH backends ----------------------------
  fmt_enrich_for_table <- function(res) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    # g:Profiler result
    if (all(c("term_name","term_size","intersection_size","p_value") %in% colnames(res))) {
      adj_col <- intersect(colnames(res), c("p_adjusted","adjusted_p_value","p_adj"))
      adj <- if (length(adj_col)) res[[adj_col[1]]] else res$p_value
      df <- data.frame(
        Term = paste0(res$term_name, " (", res$source, ")"),
        Overlap = paste0(res$intersection_size, "/", res$term_size),
        P.value = as.numeric(res$p_value),
        Adjusted.P.value = as.numeric(adj),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
      rownames(df) <- NULL
      return(df)
    }
    
    # Enrichr result (already has Term/Overlap/P.value/Adjusted.P.value)
    if (all(c("Term","Overlap","P.value","Adjusted.P.value") %in% colnames(res))) {
      df <- res[, c("Term","Overlap","P.value","Adjusted.P.value"), drop = FALSE]
      df$Term <- paste0(df$Term, " (", attr(res, "source", exact = TRUE) %||% "ENRICHR", ")")
      df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
      rownames(df) <- NULL
      return(df)
    }
    
    # If Enrichr rbind lost the 'source' attribute, try to append from raw col if present
    if ("source" %in% colnames(res) && all(c("Term","Overlap","P.value","Adjusted.P.value") %in% colnames(res))) {
      res$Term <- paste0(res$Term, " (", res$source, ")")
      df <- res[, c("Term","Overlap","P.value","Adjusted.P.value"), drop = FALSE]
      df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
      rownames(df) <- NULL
      return(df)
    }
    
    NULL
  }
  
  # ---- single entry point to run enrichment with fallback --------------------
  do_enrichment <- function(genes, species_label) {
    genes <- strip_bg(genes)
    if (length(genes) < 3) return(NULL)
    
    org <- gp_org(species_label)
    if (is.null(org)) return(NULL)
    
    gp_res <- safe_gprof(genes, org)
    if (!is.null(gp_res)) return(fmt_enrich_for_table(gp_res))
    
    # fallback to Enrichr
    enr_res <- run_enrichr_fallback(genes, org)
    fmt_enrich_for_table(enr_res)
  }
  
  # ---- barplot used by All/Up/Down ------------------------------------------
  render_enrich_barplot <- function(df, title_txt = "Top Enriched Terms", top_n = 10) {
    req(df, nrow(df) > 0)
    keep <- df[, c("Term","Adjusted.P.value"), drop = FALSE]
    keep <- keep[is.finite(keep$Adjusted.P.value) & keep$Adjusted.P.value > 0, , drop = FALSE]
    if (nrow(keep) == 0) return(NULL)
    keep$mlog10 <- -log10(keep$Adjusted.P.value)
    keep <- keep[order(keep$Adjusted.P.value, -keep$mlog10), , drop = FALSE]
    keep <- head(keep, top_n)
    
    plotly::plot_ly(
      keep,
      x = ~mlog10,
      y = ~reorder(Term, mlog10),
      type = "bar",
      orientation = "h"
    ) %>%
      plotly::layout(
        title = title_txt,
        xaxis = list(title = "-log10(FDR)"),
        yaxis = list(title = "", automargin = TRUE)
      )
  }
  
run_gprof <- function(genes, org) {
  gprofiler2::gost(
    query = unique(genes),
    organism = org,
    sources = c("GO:BP","GO:MF","GO:CC","REAC","KEGG","MSigDBH","MSigDBC2"),
    correction_method = "fdr"
  )$result
}

  
  
  
  run_enrichment <- function(genes, mode = c("all","up","down")) {
    mode <- match.arg(mode)
    backend <- input$enrich_backend
    
    if (backend == "gprof") {
      # g:Profiler (multi-species)
      org <- input$species
      srcs <- input$gprof_sources
      if (is.null(srcs) || length(srcs) == 0) srcs <- gprof_default_sources
      
      gp <- tryCatch({
        gprofiler2::gost(query = genes, organism = org, sources = srcs,
                         correction_method = "fdr", evcodes = FALSE)
      }, error = function(e) NULL)
      
      res <- if (!is.null(gp) && !is.null(gp$result)) gp$result else data.frame()
      res
    } else {
      # Enrichr path (mostly human/mouse; may be sparse for fly/zfish)
      db <- input$enrichr_db
      er <- tryCatch({ enrichr(genes, db) }, error = function(e) NULL)
      if (is.null(er)) return(data.frame())
      er[[db]]
    }
  }
  # ---------- Seurat object ----------
  observeEvent(input$create_seurat, {
    req(input$counts_file, input$meta_file)
    showNotification("📦 Creating Seurat object...", type = "message")
    
    counts <- read.csv(input$counts_file$datapath, row.names = 1, check.names = FALSE)
    meta   <- read.csv(input$meta_file$datapath,   row.names = 1, check.names = FALSE)
    
    # keep your display column name used elsewhere
    if ("cell_type" %in% colnames(meta) && !"seurat_annotations" %in% colnames(meta)) {
      meta$seurat_annotations <- meta$cell_type
    }
    
    seu <- CreateSeuratObject(counts = counts, meta.data = meta)
    seu <- NormalizeData(seu)
    
    # ✅ make sure DE logic sees real factor levels
    if ("cell_type" %in% colnames(seu@meta.data)) {
      seu$cell_type <- factor(trimws(as.character(seu$cell_type)))
    }
    if ("stim" %in% colnames(seu@meta.data)) {
      seu$stim <- factor(trimws(as.character(seu$stim)))
    }
    
    rv$seurat <- seu
    showNotification("✅ Seurat object created and normalized.", type = "message")
  })
  
  

  # ---------- DE & filtering ----------
  observeEvent(input$run_condition_de, {
    req(rv$seurat)
    showNotification("🧪 Running DE by Condition...", type = "message")
    
    if (!all(c("stim","cell_type") %in% colnames(rv$seurat@meta.data))) {
      showNotification("Metadata must contain 'stim' and 'cell_type'.", type = "error")
      return(NULL)
    }
    
    # ensure factors (in case object loaded from cache)
    rv$seurat$stim      <- factor(trimws(as.character(rv$seurat$stim)))
    rv$seurat$cell_type <- factor(trimws(as.character(rv$seurat$cell_type)))
    
    # prefer DE within most abundant cell type
    tab    <- sort(table(rv$seurat$cell_type), decreasing = TRUE)
    top_ct <- names(tab)[1]
    
    rv$seurat$celltype.stim <- paste(rv$seurat$cell_type, rv$seurat$stim, sep = "_")
    Idents(rv$seurat) <- "celltype.stim"
    
    # levels present *within* top_ct (drop unused)
    stim_lv_top <- levels(droplevels(rv$seurat$stim[rv$seurat$cell_type == top_ct]))
    
    if (length(stim_lv_top) >= 2) {
      ident1 <- paste0(top_ct, "_", stim_lv_top[1])
      ident2 <- paste0(top_ct, "_", stim_lv_top[2])
      valid  <- levels(Idents(rv$seurat))
      
      if (ident1 %in% valid && ident2 %in% valid) {
        de <- FindMarkers(rv$seurat, ident.1 = ident1, ident.2 = ident2, verbose = FALSE)
        showNotification(sprintf("DE within '%s': %s vs %s", top_ct, stim_lv_top[1], stim_lv_top[2]), type = "message")
      } else {
        Idents(rv$seurat) <- "stim"
        lv <- levels(Idents(rv$seurat))
        if (length(lv) < 2) { showNotification("Need ≥2 levels in 'stim' to run DE.", type = "error"); return(NULL) }
        de <- FindMarkers(rv$seurat, ident.1 = lv[1], ident.2 = lv[2], verbose = FALSE)
        showNotification(sprintf("Fallback: global DE by 'stim': %s vs %s", lv[1], lv[2]), type = "warning")
      }
    } else {
      # global fallback
      Idents(rv$seurat) <- "stim"
      lv <- levels(Idents(rv$seurat))
      if (length(lv) < 2) { showNotification("Need ≥2 levels in 'stim' to run DE.", type = "error"); return(NULL) }
      de <- FindMarkers(rv$seurat, ident.1 = lv[1], ident.2 = lv[2], verbose = FALSE)
      showNotification(sprintf("Global DE by 'stim': %s vs %s", lv[1], lv[2]), type = "message")
    }
    
    de$gene <- rownames(de)
    rv$condition_de <- de
    showNotification("✅ Finished DE by Condition.", type = "message")
  })
  
  
  observeEvent(input$run_celltype_de, {
    req(rv$seurat)
    showNotification("🧬 Running DE by Cell Type...", type = "message")
    
    if (!"cell_type" %in% colnames(rv$seurat@meta.data)) {
      showNotification("Missing 'cell_type' in metadata.", type = "error"); return(NULL)
    }
    
    Idents(rv$seurat) <- "cell_type"
    lv <- levels(Idents(rv$seurat))
    if (length(lv) < 2) { showNotification("Need ≥2 cell types to run DE.", type = "error"); return(NULL) }
    
    de <- FindMarkers(rv$seurat, ident.1 = lv[1], ident.2 = lv[2], verbose = FALSE)
    de$gene <- rownames(de)
    rv$celltype_de <- de
    showNotification("✅ Finished DE by Cell Type.", type = "message")
  })
  
  
  observeEvent(input$filter_condition_only, {
    req(rv$condition_de, rv$celltype_de)
    showNotification("🔍 Filtering Condition-Only DE Genes...", type = "message")
    
    cond_genes     <- rownames(subset(rv$condition_de, p_val_adj < 0.05))
    celltype_genes <- rownames(subset(rv$celltype_de, p_val_adj < 0.05))
    only_cond      <- setdiff(cond_genes, celltype_genes)
    
    rv$condition_only <- rv$condition_de[rownames(rv$condition_de) %in% only_cond, , drop = FALSE]
    showNotification("✅ Condition-Only DE Genes filtered.", type = "message")
  })
  
  
  
  # ---------- Tables ----------
  output$condition_de_table <- renderDT({ req(rv$condition_de); datatable(rv$condition_de, options = list(pageLength = 10)) })
  output$celltype_de_table  <- renderDT({ req(rv$celltype_de); datatable(rv$celltype_de, options = list(pageLength = 10)) })
  output$condition_only_table <- renderDT({ req(rv$condition_only); datatable(rv$condition_only, options = list(pageLength = 10)) })
  observeEvent(input$run_pca, {
    req(rv$seurat)
    showNotification("🔎 Running PCA...", type = "message")
    
    # pick a sane number of features/PCs based on data size
    nfeat <- min(2000, nrow(rv$seurat))
    rv$seurat <- FindVariableFeatures(rv$seurat, selection.method = "vst", nfeatures = nfeat)
    rv$seurat <- ScaleData(rv$seurat, verbose = FALSE)
    
    npcs <- max(10, min(50, length(VariableFeatures(rv$seurat))))
    rv$seurat <- RunPCA(rv$seurat, features = VariableFeatures(rv$seurat), npcs = npcs, verbose = FALSE)
    
    showNotification("✅ PCA completed.", type = "message")
  })
  
  observeEvent(input$run_umap, {
    req(rv$seurat)
    showNotification("📉 Running UMAP...", type = "message")
    
    # ensure PCA exists
    if (!"pca" %in% Reductions(rv$seurat)) {
      nfeat <- min(2000, nrow(rv$seurat))
      rv$seurat <- FindVariableFeatures(rv$seurat, selection.method = "vst", nfeatures = nfeat)
      rv$seurat <- ScaleData(rv$seurat, verbose = FALSE)
      npcs <- max(10, min(50, length(VariableFeatures(rv$seurat))))
      rv$seurat <- RunPCA(rv$seurat, features = VariableFeatures(rv$seurat), npcs = npcs, verbose = FALSE)
    }
    
    use_dims <- 1:min(30, ncol(Embeddings(rv$seurat, "pca")))
    rv$seurat <- RunUMAP(rv$seurat, reduction = "pca", dims = use_dims, verbose = FALSE)
    
    showNotification("✅ UMAP completed.", type = "message")
  })
  
  
  
  output$umap_plot <- renderPlotly({
    req(rv$seurat); req("umap" %in% Reductions(rv$seurat))
    p <- DimPlot(rv$seurat, reduction = "umap", group.by = "seurat_annotations", label = TRUE) + theme_minimal()
    ggplotly(p)
  })
  output$pca_plot <- renderPlotly({
    req(rv$seurat); req("pca" %in% Reductions(rv$seurat))
    p <- DimPlot(rv$seurat, reduction = "pca", group.by = "seurat_annotations", label = TRUE) + theme_minimal()
    ggplotly(p)
  })
 
# ALL
observeEvent(input$enrich_all, {
  req(rv$condition_only)
  showNotification("🔬 Enriching ALL DE genes...", type = "message")
  genes <- rownames(rv$condition_only)
  tbl <- do_enrichment(genes, input$species)
  rv$enrich$all <- tbl
  rv$enrich_edges$all <- make_edges(tbl)   # 🆕 term–gene edges
  if (is.null(tbl)) showNotification("No terms returned (ALL).", type = "warning")
  showNotification("✅ Enrichment complete (ALL).", type = "message")
})

  
# UP
observeEvent(input$enrich_up, {
  req(rv$condition_only)
  showNotification("🔬 Enriching UP genes...", type = "message")
  de <- rv$condition_only
  genes <- rownames(de[de$avg_log2FC > 0, , drop = FALSE])
  tbl <- do_enrichment(genes, input$species)
  rv$enrich$up <- tbl
  rv$enrich_edges$up <- make_edges(tbl)    # 🆕
  if (is.null(tbl)) showNotification("No terms returned (UP).", type = "warning")
  showNotification("✅ Enrichment complete (UP).", type = "message")
})

 # DOWN
observeEvent(input$enrich_down, {
  req(rv$condition_only)
  showNotification("🔬 Enriching DOWN genes...", type = "message")
  de <- rv$condition_only
  genes <- rownames(de[de$avg_log2FC < 0, , drop = FALSE])
  tbl <- do_enrichment(genes, input$species)
  rv$enrich$down <- tbl
  rv$enrich_edges$down <- make_edges(tbl)  # 🆕
  if (is.null(tbl)) showNotification("No terms returned (DOWN).", type = "warning")
  showNotification("✅ Enrichment complete (DOWN).", type = "message")
})
 

  
  render_enrich_barplot <- function(df, title_txt = "Top Enriched Terms", top_n = 10) {
    req(df, nrow(df) > 0)
    keep <- df[, c("Term","Adjusted.P.value"), drop = FALSE]
    keep <- keep[is.finite(keep$Adjusted.P.value) & keep$Adjusted.P.value > 0, , drop = FALSE]
    if (nrow(keep) == 0) return(NULL)
    keep$mlog10 <- -log10(keep$Adjusted.P.value)
    keep <- keep[order(keep$Adjusted.P.value, -keep$mlog10), , drop = FALSE]
    keep <- head(keep, top_n)
    
    plotly::plot_ly(
      keep,
      x = ~mlog10,
      y = ~reorder(Term, mlog10),
      type = "bar",
      orientation = "h"
    ) %>%
      plotly::layout(
        title = title_txt,
        xaxis = list(title = "-log10(FDR)"),
        yaxis = list(title = "", automargin = TRUE)
      )
  }
  
  
  output$enrich_all_table  <- renderDT({ req(rv$enrich$all);  datatable(rv$enrich$all) })
  output$enrich_up_table   <- renderDT({ req(rv$enrich$up);   datatable(rv$enrich$up) })
  output$enrich_down_table <- renderDT({ req(rv$enrich$down); datatable(rv$enrich$down) })
  
  output$enrich_all_plot  <- renderPlotly({ req(rv$enrich$all);  render_enrich_barplot(rv$enrich$all) })
  output$enrich_up_plot   <- renderPlotly({ req(rv$enrich$up);   render_enrich_barplot(rv$enrich$up) })
  output$enrich_down_plot <- renderPlotly({ req(rv$enrich$down); render_enrich_barplot(rv$enrich$down) })
  
  # ---------- Random Forest ----------
  observeEvent(input$feature_select, {
    req(rv$seurat, rv$condition_only)
    showNotification("🌲 Running Random Forest on Condition-Only DE Genes...", type = "message")
    
    features <- rownames(rv$condition_only)
    expr <- t(as.matrix(GetAssayData(rv$seurat, assay = "RNA", slot = "data")[features, , drop = FALSE]))
    labels <- rv$seurat$stim
    df <- data.frame(expr)
    df$stim <- as.factor(labels)
    
    set.seed(42)
    train_index <- createDataPartition(df$stim, p = 0.7, list = FALSE)
    train_data <- df[train_index, , drop = FALSE]
    test_data  <- df[-train_index, , drop = FALSE]
    
    rf_model <- randomForest(stim ~ ., data = train_data, importance = TRUE, ntree = 500)
    
    importance_df <- as.data.frame(importance(rf_model))
    importance_df$Gene <- rownames(importance_df)
    importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
    rv$rf_importance <- importance_df
    
    showNotification("✅ Feature selection complete.", type = "message")
  })
  
  output$rf_table <- renderDT({ req(rv$rf_importance); datatable(rv$rf_importance, options = list(pageLength = 15), rownames = FALSE) })
  
  observeEvent(input$classify_rf, {
    req(rv$seurat, rv$rf_importance, input$top_n_features)
    showNotification("🤖 Classifying using Random Forest...", type = "message")
    
    top_genes <- head(rv$rf_importance$Gene, input$top_n_features)
    expr <- t(as.matrix(GetAssayData(rv$seurat, assay = "RNA", slot = "data")[top_genes, , drop = FALSE]))
    labels <- rv$seurat$stim
    df <- data.frame(expr)
    df$stim <- factor(labels)
    
    set.seed(101)
    split <- createDataPartition(df$stim, p = 0.7, list = FALSE)
    train <- df[split, , drop = FALSE]
    test  <- df[-split, , drop = FALSE]
    
    model <- randomForest(stim ~ ., data = train, ntree = 500)
    
    pred  <- predict(model, newdata = test, type = "response")
    probs <- tryCatch({ predict(model, newdata = test, type = "prob")[, 2] }, error = function(e) rep(NA_real_, nrow(test)))
    
    actual <- test$stim
    conf <- caret::confusionMatrix(pred, actual)
    
    if (length(levels(actual)) == 2 && !all(is.na(probs))) {
      roc_obj <- pROC::roc(actual, probs)
      auc_val <- pROC::auc(roc_obj)
    } else {
      roc_obj <- NULL
      auc_val <- NA
    }
    
    rv$classify <- list(
      predictions = data.frame(Actual = actual, Predicted = pred, Probability = probs),
      confusion = conf,
      auc = auc_val,
      roc = roc_obj
    )
    
    showNotification("✅ Classification complete!", type = "message")
  })
  
  output$rf_predictions <- renderDT({ req(rv$classify$predictions); datatable(rv$classify$predictions) })
  
  output$rf_roc <- renderPlotly({
    req(rv$classify$roc)
    roc_df <- data.frame(
      FPR = 1 - rv$classify$roc$specificities,
      TPR = rv$classify$roc$sensitivities
    )
    plot_ly(roc_df, x = ~FPR, y = ~TPR, type = 'scatter', mode = 'lines') %>%
      layout(title = "ROC Curve", xaxis = list(title = "False Positive Rate"), yaxis = list(title = "True Positive Rate"), showlegend = FALSE)
  })
  
  output$rf_metrics <- renderDT({
    req(rv$classify$confusion)
    metrics_df <- as.data.frame(t(rv$classify$confusion$byClass))
    metrics_df$Metric <- rownames(metrics_df)
    rownames(metrics_df) <- NULL
    metrics_df <- metrics_df[, c(ncol(metrics_df), 1:(ncol(metrics_df)-1))]
    if (!is.na(rv$classify$auc)) {
      metrics_df$AUC <- NA
      metrics_df$AUC[1] <- rv$classify$auc
    }
    datatable(metrics_df, options = list(dom = 't', pageLength = nrow(metrics_df)))
  })
  
  # ---------- Power analysis ----------
  observeEvent(input$run_power, {
    req(rv$seurat, rv$condition_only)
    showNotification("🔬 Running Power Analysis on Top 10 DE Genes...", type = "message")
    
    de_tab <- rv$condition_only
    top_genes <- head(rownames(de_tab[order(de_tab$p_val_adj), , drop = FALSE]), 10)
    labels <- rv$seurat$stim
    features <- GetAssayData(rv$seurat, assay = "RNA", slot = "data")
    
    power_list <- lapply(top_genes, function(gene) {
      expr <- as.numeric(features[gene, ])
      lv <- unique(labels)
      if (length(lv) < 2) return(NULL)
      g1 <- expr[labels == lv[1]]
      g2 <- expr[labels == lv[2]]
      if (length(g1) < 2 || length(g2) < 2) return(NULL)
      eff_size <- abs(mean(g1) - mean(g2)) / sqrt((var(g1) + var(g2))/2)
      n1 <- length(g1); n2 <- length(g2)
      pwr_res <- pwr::pwr.t.test(n = min(n1, n2), d = eff_size, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
      n_seq <- seq(5, max(100, max(n1, n2) * 2), by = 1)
      powers <- sapply(n_seq, function(n) pwr::pwr.t.test(n = n, d = eff_size, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power)
      min_n_08 <- n_seq[which(powers >= 0.8)[1]]
      list(Gene = gene, EffectSize = round(eff_size, 3), ObservedN1 = n1, ObservedN2 = n2,
           ObservedPower = round(pwr_res$power, 3), N_for_0.8 = min_n_08,
           CurveX = n_seq, CurveY = powers)
    })
    
    power_list <- Filter(Negate(is.null), power_list)
    tab <- data.frame(
      Gene = sapply(power_list, `[[`, "Gene"),
      EffectSize = sapply(power_list, `[[`, "EffectSize"),
      ObservedN1 = sapply(power_list, `[[`, "ObservedN1"),
      ObservedN2 = sapply(power_list, `[[`, "ObservedN2"),
      ObservedPower = sapply(power_list, `[[`, "ObservedPower"),
      N_for_0.8 = sapply(power_list, `[[`, "N_for_0.8")
    )
    rv$power_table <- tab
    rv$power_curves <- setNames(lapply(power_list, function(pl) data.frame(SampleSizePerGroup = pl$CurveX, Power = pl$CurveY)), sapply(power_list, `[[`, "Gene"))
    updateSelectInput(session, "power_gene", choices = tab$Gene, selected = tab$Gene[1])
  })
  
  output$power_table <- renderDT({ req(rv$power_table); datatable(rv$power_table, options = list(dom = 't')) })
  
  output$power_curve <- renderPlotly({
    req(rv$power_curves, input$power_gene)
    curve_df <- rv$power_curves[[input$power_gene]]
    plot_ly(curve_df, x = ~SampleSizePerGroup, y = ~Power, type = "scatter", mode = "lines") %>%
      layout(title = paste("Power Curve:", input$power_gene),
             xaxis = list(title = "Sample Size per Group"),
             yaxis = list(title = "Power"),
             shapes = list(list(type = "line", x0 = min(curve_df$SampleSizePerGroup), x1 = max(curve_df$SampleSizePerGroup), y0 = 0.8, y1 = 0.8, line = list(dash = "dash"))))
  })
  
  # ---------- Downloads ----------
  output$download_condition_de <- downloadHandler(filename = function() {"condition_de_table.csv"}, content = function(file) { write.csv(rv$condition_de, file) })
  output$download_celltype_de  <- downloadHandler(filename = function() {"celltype_de_table.csv"},  content = function(file) { write.csv(rv$celltype_de, file) })
  output$download_condition_only <- downloadHandler(filename = function() {"condition_only_de_table.csv"}, content = function(file) { write.csv(rv$condition_only, file) })
  
  output$download_enrich_all  <- downloadHandler(filename = function() {"enrichment_all.csv"},  content = function(file) { write.csv(rv$enrich$all, file) })
  output$download_enrich_up   <- downloadHandler(filename = function() {"enrichment_up.csv"},   content = function(file) { write.csv(rv$enrich$up, file) })
  output$download_enrich_down <- downloadHandler(filename = function() {"enrichment_down.csv"}, content = function(file) { write.csv(rv$enrich$down, file) })
  
  output$download_rf_importance <- downloadHandler(filename = function() {"rf_importance_table.csv"}, content = function(file) { write.csv(rv$rf_importance, file) })
  output$download_rf_predictions <- downloadHandler(filename = function() {"rf_predictions_table.csv"}, content = function(file) { write.csv(rv$classify$predictions, file) })
  output$download_rf_metrics <- downloadHandler(filename = function() {"rf_metrics_table.csv"}, content = function(file) {
    if (!is.null(rv$classify$confusion)) {
      metrics_df <- as.data.frame(t(rv$classify$confusion$byClass))
      metrics_df$Metric <- rownames(metrics_df)
      rownames(metrics_df) <- NULL
      metrics_df <- metrics_df[, c(ncol(metrics_df), 1:(ncol(metrics_df)-1))]
      if (!is.na(rv$classify$auc)) { metrics_df$AUC <- NA; metrics_df$AUC[1] <- rv$classify$auc }
      write.csv(metrics_df, file)
    } else {
      write.csv(data.frame(), file)
    }
  })
  
  output$download_power_table <- downloadHandler(filename = function() {"power_table.csv"}, content = function(file) { write.csv(rv$power_table, file) })
  
  # ---------- Volcano ----------
  observeEvent(input$make_volcano, {
    req(rv$condition_only)
    df <- rv$condition_only
    df$log10p <- -log10(df$p_val_adj + 1e-300)  # Avoid -Inf
    df$gene <- rownames(df)
    rv$volcano <- df
    showNotification("Volcano plot generated!", type = "message")
  })
  
  output$condition_only_volcano <- renderPlotly({
    req(rv$volcano)
    plot_ly(rv$volcano, x = ~avg_log2FC, y = ~log10p, text = ~gene, type = "scatter", mode = "markers",
            marker = list(size = 8, opacity = 0.6)) %>%
      layout(title = "Volcano Plot: Condition-Only DE Genes",
             xaxis = list(title = "Log2 Fold Change"),
             yaxis = list(title = "-log10 Adjusted P-value"))
  })
  
  observeEvent(input$run_heatmap, {
    req(rv$seurat, rv$rf_importance)
    showNotification("🔥 Generating heatmap from top 10 selected genes...", type = "message")
    
    # 1) Top 10 features by importance
    imp <- rv$rf_importance
    stopifnot(all(c("Gene","MeanDecreaseGini") %in% colnames(imp)))
    imp <- imp[order(imp$MeanDecreaseGini, decreasing = TRUE), , drop = FALSE]
    top_fs <- head(as.character(imp$Gene), 10)
    
    # 2) Prefer genes that are also in the Condition-Only table (if present)
    cond_only_genes <- if (!is.null(rv$condition_only)) rownames(rv$condition_only) else character(0)
    
    # case-insensitive overlap to avoid Cxcl10/CXCL10 mismatches
    lc <- function(x) tolower(trimws(as.character(x)))
    preferred <- top_fs[lc(top_fs) %in% lc(cond_only_genes)]
    
    # 3) If not enough preferred genes, fill with remaining top features present in the Seurat object
    present_in_obj <- top_fs[lc(top_fs) %in% lc(rownames(rv$seurat))]
    genes <- unique(c(preferred, present_in_obj))
    genes <- genes[genes %in% rownames(rv$seurat)]
    
    if (length(genes) < 3) {
      showNotification("Not enough overlapping top features to plot heatmap.", type = "error")
      return(NULL)
    }
    
    # 4) Scale data for these genes (quietly)
    DefaultAssay(rv$seurat) <- "RNA"
    try({ rv$seurat <- ScaleData(rv$seurat, features = genes, verbose = FALSE) }, silent = TRUE)
    
    # 5) Group by stim (condition)
    cond_col <- if ("stim" %in% colnames(rv$seurat@meta.data)) "stim" else NULL
    if (is.null(cond_col)) { 
      showNotification("Missing 'stim' column to group heatmap.", type = "error")
      return(NULL)
    }
    Idents(rv$seurat) <- cond_col
    
    # 6) Plot (Plotly in your RF » Heatmap tab)
    p <- DoHeatmap(rv$seurat, features = genes, group.by = cond_col) +
      ggtitle(sprintf("Heatmap • Top features (preferring condition-only) • n=%d", length(genes))) +
      theme(axis.text.y = element_text(size = 7))
    
    output$rf_heatmap <- renderPlotly({ ggplotly(p) })
    showNotification("✅ Heatmap generated.", type = "message")
  })
  

output$download_enrich_all_edges <- downloadHandler(
  filename = function() {"enrichment_all_edges.csv"},
  content = function(file) {
    write.csv(rv$enrich_edges$all, file, row.names = FALSE)
  }
)

output$download_enrich_up_edges <- downloadHandler(
  filename = function() {"enrichment_up_edges.csv"},
  content = function(file) {
    write.csv(rv$enrich_edges$up, file, row.names = FALSE)
  }
)

output$download_enrich_down_edges <- downloadHandler(
  filename = function() {"enrichment_down_edges.csv"},
  content = function(file) {
    write.csv(rv$enrich_edges$down, file, row.names = FALSE)
  }
)

output$download_enrich_all_edges <- downloadHandler(
  filename = function() {"enrichment_all_edges.csv"},
  content = function(file) {
    write.csv(rv$enrich_edges$all, file, row.names = FALSE)
  }
)

output$download_enrich_up_edges <- downloadHandler(
  filename = function() {"enrichment_up_edges.csv"},
  content = function(file) {
    write.csv(rv$enrich_edges$up, file, row.names = FALSE)
  }
)

output$download_enrich_down_edges <- downloadHandler(
  filename = function() {"enrichment_down_edges.csv"},
  content = function(file) {
    write.csv(rv$enrich_edges$down, file, row.names = FALSE)
  }
)


  
}



shinyApp(ui, server)



