library(shiny)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(plotly)
library(RColorBrewer)
library(grid)

options(shiny.maxRequestSize=60*1024^2)

# Define UI
ui <- fluidPage(
  titlePanel("HeraNorm Endogenous Control & DEG Analysis Tool"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("countsFile", "Upload Count Matrix (CSV format)", accept = c("text/csv", ".csv")),
      fileInput("groupFile", "Upload Group Information (TXT format)", accept = c("text/plain", ".txt")),
      actionButton("analyze", "Run Analysis"),
      hr(),
      
      # Filtering parameters
      h4("Filtering Parameters: Endogenous Control"),
      numericInput("endoPvalue", "Minimum p-value threshold (default: 0.7)", value = 0.7, min = 0, max = 1, step = 0.01),
      numericInput("endoBaseMean", "Minimum baseMean value (default: 100)", value = 100, min = 0),
      numericInput("endoLog2FCMin", "Minimum log2FoldChange (default: -0.2)", value = -0.2, step = 0.01),
      numericInput("endoLog2FCMax", "Maximum log2FoldChange (default: 0.2)", value = 0.2, step = 0.01),
      hr(),
      
      h4("Filtering Parameters: DEG"),
      numericInput("degPvalue", "Maximum p-value threshold (default: 0.05)", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("degBaseMean", "Minimum baseMean value (default: 100)", value = 100, min = 0),
      numericInput("degLog2FC", "Minimum absolute log2FoldChange (default: 0.2)", value = 0.2, min = 0, step = 0.01),
      hr(),
      
      # Plotting options
      h4("Plotting Options"),
      selectInput("plotGeneType", "Select Gene Type", choices = c("DEG", "Endogenous Control")),
      uiOutput("geneSelector"),
      actionButton("plotGene", "Plot Gene Expression"),
      hr(),
      
      h4("in silico normalization"),
      uiOutput("endoControlSelector"),
      uiOutput("degSelector"),
      actionButton("plotRelative", "Plot in silico normalization"),
      hr(),
      
      # Heatmap parameters
      h4("Heatmap Settings"),
      numericInput("heatmapNumGenes", "Number of DEGs to Display", value = 50, min = 10, max = 500),
      selectInput("heatmapColor", "Color Scheme", 
                  choices = c("Blue-White-Red", "Purple-White-Green", "Blue-Yellow-Red"), 
                  selected = "Blue-White-Red"),
      hr(),
      
      # Volcano plot parameters
      h4("Volcano Plot Settings"),
      checkboxInput("showLabels", "Show Top Gene Labels", value = TRUE),
      numericInput("topGeneNum", "Number of Top Genes to Label", value = 10, min = 1, max = 50, step = 1),
      hr(),
      
      # Download buttons
      downloadButton("downloadEndogenous", "Download Endogenous Control Table"),
      downloadButton("downloadDeg", "Download DEG Table"),
      downloadButton("downloadGenePlot", "Download Gene Expression Plot (PDF)"),
      downloadButton("downloadRelativePlot", "Download In Silico normalization (PDF)"),
      downloadButton("downloadHeatmap", "Download Heatmap (PDF)"),
      downloadButton("downloadVolcano", "Download Volcano Plot (PDF)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Endogenous Control Genes", tableOutput("endogenousTable")),
        tabPanel("DEG Genes", tableOutput("degTable")),
        tabPanel("Gene Expression Plot", plotOutput("genePlot")),
        tabPanel("In Silico Normalization Plot", plotOutput("relativePlot")),
        tabPanel("DEG Heatmap", plotOutput("heatmap", height = "600px")),
        tabPanel("Volcano Plot", plotlyOutput("volcanoPlot", height = "600px"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  
  # Create reactive values to store analysis results
  analysisResults <- reactiveValues(data = NULL)
  
  # Replace the eventReactive with progress tracking
  observeEvent(input$analyze, {
    # Only proceed if files are uploaded
    req(input$countsFile, input$groupFile)
    
    # Start progress bar
    withProgress(message = 'Analyzing data...', value = 0, {
      # reading files
      setProgress(message = "Reading input files...")
      incProgress(0.1)
      
      # Read data
      counts_matrix <- read.csv(input$countsFile$datapath, row.names = 1)
      sample_info <- read.table(input$groupFile$datapath, header = TRUE, sep = "\t", row.names = 1)
      
      # Check sample names
      if (!all(rownames(sample_info) %in% colnames(counts_matrix))) {
        showNotification("Sample names in group file don't match count matrix columns!", type = "error")
        return(NULL)
      }
      counts_matrix <- counts_matrix[, rownames(sample_info)]
      
      # creating dataset
      setProgress(message = "Creating dataset...")
      incProgress(0.2)
      
      # Create DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = counts_matrix, 
                                    colData = sample_info, 
                                    design = ~ condition)
      dds <- dds[rowSums(counts(dds)) > 1, ]
      
      # running DESeq analysis (most time-consuming step)
      setProgress(message = "Running DESeq analysis...")
      incProgress(0.3)
      
      # Run DESeq analysis
      dds <- DESeq(dds)
      
      # processing results
      setProgress(message = "Processing results...")
      incProgress(0.2)
      
      # Get normalized counts and results
      normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
      raw_counts <- as.data.frame(counts(dds))
      res <- results(dds)
      
      # Process results
      res_df <- as.data.frame(res) %>%
        mutate(
          gene_name = rownames(res),
          endogenous_control = ifelse(
            pvalue > input$endoPvalue &
              baseMean > input$endoBaseMean &
              log2FoldChange > input$endoLog2FCMin &
              log2FoldChange < input$endoLog2FCMax,
            TRUE, FALSE),
          DEG = ifelse(
            pvalue < input$degPvalue &
              baseMean > input$degBaseMean &
              abs(log2FoldChange) > input$degLog2FC,
            TRUE, FALSE)
        )
      vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
      vst_data <- assay(vsd) 
      # Store results
      analysisResults$data <- list(
        vst_data = vst_data,
        raw_counts = raw_counts,
        normalized_counts = normalized_counts,
        sample_info = sample_info,
        res_df = res_df,
        endogenous_control_data = res_df %>% filter(endogenous_control == TRUE),
        deg_data = res_df %>% filter(DEG == TRUE)
      )
      
      # Set progress to 100%
      setProgress(1, message = 'Analysis complete!')
      showNotification("Analysis completed successfully!", type = "message")
    })
  })
  
  # Create a reactive expression to access analysis results
  dataInput <- reactive({
    req(analysisResults$data)
    analysisResults$data
  })
  
  # Dynamic UI elements
  output$geneSelector <- renderUI({
    data <- dataInput()
    req(data)
    gene_list <- if (input$plotGeneType == "DEG") {
      data$deg_data$gene_name
    } else {
      data$endogenous_control_data$gene_name
    }
    selectInput("selectedGene", "Select Gene", choices = gene_list)
  })
  output$endoControlSelector <- renderUI({
    data <- dataInput()
    req(data)
    selectInput("selectedEndoControl", "Select Endogenous Control Gene", 
                choices = data$endogenous_control_data$gene_name)
  })
  
  output$degSelector <- renderUI({
    data <- dataInput()
    req(data)
    selectInput("selectedDeg", "Select Target Gene (DEG)", 
                choices = data$deg_data$gene_name)
  })
  
  # Table outputs
  output$endogenousTable <- renderTable({
    data <- dataInput()
    req(data)
    data$endogenous_control_data %>%
      select(gene_name, baseMean, log2FoldChange, pvalue, padj)
  })
  
  output$degTable <- renderTable({
    data <- dataInput()
    req(data)
    data$deg_data %>%
      select(gene_name, baseMean, log2FoldChange, pvalue, padj)
  })
  
  # Download handlers
  output$downloadEndogenous <- downloadHandler(
    filename = function() { paste("endogenous_controls_", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      data <- dataInput()
      req(data)
      write.csv(data$endogenous_control_data, file, row.names = FALSE)
    }
  )
  
  output$downloadDeg <- downloadHandler(
    filename = function() { paste("DEGs_", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      data <- dataInput()
      req(data)
      write.csv(data$deg_data, file, row.names = FALSE)
    }
  )
  
  # Plot data preparation
  genePlotData <- reactive({
    data <- dataInput()
    req(data, input$selectedGene)
    gene <- input$selectedGene
    df <- data.frame(
      Sample = colnames(data$normalized_counts),
      Expression = as.numeric(data$normalized_counts[gene, ]),
      Group = data$sample_info$condition
    )
    df[complete.cases(df), ]
  })
  
  relativePlotData <- reactive({
    data <- dataInput()
    req(data, input$selectedEndoControl, input$selectedDeg)
    df <- data.frame(
      Sample = colnames(data$normalized_counts),
      Endo = as.numeric(data$normalized_counts[input$selectedEndoControl, ]),
      DEG = as.numeric(data$normalized_counts[input$selectedDeg, ]),
      Group = data$sample_info$condition
    ) %>%
      mutate(Relative = DEG / Endo)
    df[complete.cases(df), ]
  })
  
  # Plot rendering
  output$genePlot <- renderPlot({
    df <- genePlotData()
    wilcox_p <- wilcox.test(Expression ~ Group, data = df)$p.value
    
    ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = paste("Expression of", input$selectedGene),
           subtitle = paste("Wilcoxon p-value:", format.pval(wilcox_p, digits = 3)),
           x = "Group", y = "Normalized Expression") +
      theme_minimal() +
      theme(legend.position = "none",
            aspect.ratio = 1,
            text = element_text(size = 16))
  })
  
  output$relativePlot <- renderPlot({
    df <- relativePlotData()
    wilcox_p <- wilcox.test(Relative ~ Group, data = df)$p.value
    
    ggplot(df, aes(x = Group, y = Relative, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = paste(input$selectedDeg, "/", input$selectedEndoControl),
           subtitle = paste("Wilcoxon p-value:", format.pval(wilcox_p, digits = 3)),
           x = "Group", y = "Relative Expression") +
      theme_minimal() +
      theme(legend.position = "none",
            aspect.ratio = 1,
            text = element_text(size = 16))
  })
  
  # Plot download handlers
  output$downloadGenePlot <- downloadHandler(
    filename = function() {
      paste("gene_plot_", input$selectedGene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      df <- genePlotData()
      wilcox_p <- wilcox.test(Expression ~ Group, data = df)$p.value
      
      p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        labs(title = paste("Expression of", input$selectedGene),
             subtitle = paste("Wilcoxon p-value:", format.pval(wilcox_p, digits = 3)),
             x = "Group", y = "Normalized Expression") +
        theme_minimal() +
        theme(legend.position = "none",
              aspect.ratio = 1,
              text = element_text(size = 16))
      
      ggsave(file, plot = p, device = "pdf", width = 8, height = 6)
    }
  )
  
  output$downloadRelativePlot <- downloadHandler(
    filename = function() {
      paste("relative_expression_", input$selectedDeg, "_", 
            input$selectedEndoControl, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      df <- relativePlotData()
      wilcox_p <- wilcox.test(Relative ~ Group, data = df)$p.value
      
      p <- ggplot(df, aes(x = Group, y = Relative, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        labs(title = paste(input$selectedDeg, "/", input$selectedEndoControl),
             subtitle = paste("Wilcoxon p-value:", format.pval(wilcox_p, digits = 3)),
             x = "Group", y = "Relative Expression") +
        theme_minimal() +
        theme(legend.position = "none",
              aspect.ratio = 1,
              text = element_text(size = 16))
      
      ggsave(file, plot = p, device = "pdf", width = 8, height = 6)
    }
  )
  
  # Heatmap data preparation
  # Heatmap data preparation
  heatmapData <- reactive({
    req(input$heatmapNumGenes, input$heatmapColor)
    data <- dataInput()
    req(data, nrow(data$deg_data) > 0)
    # Select top DEGs by significance
    deg_genes <- data$deg_data %>%
      arrange(pvalue) %>%
      head(input$heatmapNumGenes) %>%
      pull(gene_name)
    
    norm_data <- data$vst_data
    
    # Subset DEG data
    deg_data <- norm_data[deg_genes, , drop = FALSE]
    # Apply row Z-scoring
    scaled_data <- t(scale(t(deg_data)))
    # Get condition info
    condition_info <- data$sample_info$condition
    names(condition_info) <- rownames(data$sample_info)
    # Group
    grouped_samples <- split(colnames(deg_data), condition_info[colnames(deg_data)])
    # clustering
    clustered_samples <- lapply(grouped_samples, function(group_samples) {
      if(length(group_samples) > 2) {
        # clustering inside
        hc <- hclust(dist(t(deg_data[, group_samples, drop = FALSE])), method = "average")
        group_samples[hc$order]
      } else {
        group_samples
      }
    })
    # order_samples
    ordered_samples <- unlist(clustered_samples)
    # reorder
    scaled_data <- scaled_data[, ordered_samples, drop = FALSE]
    # group annotation
    annotation_df <- data.frame(Condition = condition_info[colnames(scaled_data)],
                                row.names = colnames(scaled_data))
    # color map
    color_map <- switch(input$heatmapColor,
                        "Blue-White-Red" = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                        "Purple-White-Green" = colorRampPalette(brewer.pal(n = 11, name = "PiYG"))(100),
                        "Blue-Yellow-Red" = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(100))
    # gaps col
    gaps_col <- cumsum(sapply(clustered_samples, length))
    if (length(gaps_col) > 1) {
      gaps_col <- gaps_col[-length(gaps_col)]
      gaps_col <- gaps_col[gaps_col < ncol(deg_data)]
    } else {
      gaps_col <- integer(0)
    }
    list(
      matrix = scaled_data,
      annotation = annotation_df,
      color_map = color_map,
      gaps_col = gaps_col
    )
  })
  # Render heatmap
  output$heatmap <- renderPlot({
    hm_data <- heatmapData()
    
    pheatmap(hm_data$matrix,
             annotation_col = hm_data$annotation,
             color = hm_data$color_map,
             show_rownames = ifelse(input$heatmapNumGenes <= 50, TRUE, FALSE),
             show_colnames = TRUE,
             cluster_rows = TRUE,
             cluster_cols = FALSE, 
             gaps_col = hm_data$gaps_col, 
             fontsize_row = ifelse(input$heatmapNumGenes <= 30, 10, 8),
             main = paste("Top", nrow(hm_data$matrix), "Differentially Expressed Genes"),
             border_color = NA)
  })
  outputOptions(output, "heatmap", suspendWhenHidden = FALSE)
  
  # Download heatmap
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste("DEG_heatmap_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      hm_data <- heatmapData()
      height_val <- max(8, 0.5 * nrow(hm_data$matrix) + 2)
      
      pdf(file, width = 10, height = height_val)
      
      hm <- pheatmap(hm_data$matrix,
                     annotation_col = hm_data$annotation,
                     color = hm_data$color_map,
                     show_rownames = ifelse(input$heatmapNumGenes <= 50, TRUE, FALSE),
                     show_colnames = TRUE,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     gaps_col = hm_data$gaps_col,
                     fontsize_row = ifelse(input$heatmapNumGenes <= 30, 10, 8),
                     main = paste("Top", nrow(hm_data$matrix), "Differentially Expressed Genes"),
                     border_color = NA)
      
      if (requireNamespace("grid", quietly = TRUE)) {
        grid::grid.newpage()
        grid::grid.draw(hm$gtable)
      } else {
        print(hm$gtable)
      }
      
      dev.off()
    }
  )
  
  
  # Volcano plot data preparation
  volcanoData <- reactive({
    data <- dataInput()
    req(data)
    
    data$res_df %>%
      mutate(
        significance = case_when(
          
          pvalue < input$degPvalue &
            baseMean > input$degBaseMean &
            abs(log2FoldChange) > input$degLog2FC ~ ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"),
          TRUE ~ "Not Significant"
        ),
        log10p = -log10(pvalue + 1e-300)  # Handle p-value = 0
      )
  })
  
  # Render interactive volcano plot
  output$volcanoPlot <- renderPlotly({
    vdata <- volcanoData()
    req(vdata)
    p <- ggplot(vdata, aes(x = log2FoldChange, y = log10p, 
                           color = significance, 
                           text = paste("Gene:", gene_name,
                                        "<br>log2FC:", round(log2FoldChange, 2),
                                        "<br>p-value:", format.pval(pvalue, digits = 3),
                                        "<br>baseMean:", round(baseMean, 1)))) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_vline(xintercept = c(-input$degLog2FC, input$degLog2FC), 
                 linetype = "dashed", color = "black", alpha = 0.5) +
      geom_hline(yintercept = -log10(input$degPvalue + 1e-300), 
                 linetype = "dashed", color = "black", alpha = 0.5) +
      scale_color_manual(values = c(
        "Up-regulated" = "red",
        "Down-regulated" = "blue",
        "Not Significant" = "gray"
      )) +
      labs(x = "log2 Fold Change", y = "-log10(p-value)",
           title = "Volcano Plot of Differential Expression",
           color = "Gene Type") +
      theme_minimal() +
      theme(legend.position = "bottom")
    # Add labels for top DEGs
    if (input$showLabels && nrow(vdata) > 0) {
      top_genes <- vdata %>%
        filter(significance %in% c("Up-regulated", "Down-regulated")) %>%
        arrange(pvalue) %>%
        head(input$topGeneNum)
      p <- p + 
        geom_text(data = top_genes, 
                  aes(label = gene_name), 
                  color = "black", size = 3, vjust = -0.5, show.legend = FALSE)
    }
    # Convert to plotly
    ggplotly(p, tooltip = "text") %>%
      layout(hoverlabel = list(bgcolor = "white", font = list(size = 12)),
             legend = list(orientation = "h", y = -0.2))
  })
  
  # Download volcano plot as PDF
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste("volcano_plot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      vdata <- volcanoData()
      p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue + 1e-300), 
                             color = significance)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_vline(xintercept = c(-input$degLog2FC, input$degLog2FC), 
                   linetype = "dashed", color = "black", alpha = 0.5) +
        geom_hline(yintercept = -log10(input$degPvalue + 1e-300), 
                   linetype = "dashed", color = "black", alpha = 0.5) +
        scale_color_manual(values = c(
          "Up-regulated" = "red",
          "Down-regulated" = "blue",
          "Not Significant" = "gray"
        )) +
        labs(x = "log2 Fold Change", y = "-log10(p-value)",
             title = "Volcano Plot of Differential Expression",
             color = "Gene Type") +
        theme_minimal() +
        theme(legend.position = "bottom")
      # Add labels for top DEGs in PDF version
      if (input$showLabels && nrow(vdata) > 0) {
        top_genes <- vdata %>%
          filter(significance %in% c("Up-regulated", "Down-regulated")) %>%
          arrange(pvalue) %>%
          head(input$topGeneNum)
        p <- p + 
          ggrepel::geom_text_repel(data = top_genes, 
                                   aes(label = gene_name), 
                                   color = "black", size = 3, max.overlaps = 20)
      }
      ggsave(file, plot = p, device = "pdf", width = 10, height = 8)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
