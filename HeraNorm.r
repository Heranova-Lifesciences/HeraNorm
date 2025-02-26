library(shiny)
library(DESeq2)
library(dplyr)
library(ggplot2)

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
      numericInput("endoPvalue", "Minimum p-value threshold (default: 0.8)", value = 0.8, min = 0, max = 1, step = 0.01),
      numericInput("endoBaseMean", "Minimum baseMean value (default: 1000)", value = 1000, min = 0),
      numericInput("endoLog2FCMin", "Minimum log2FoldChange (default: -0.02)", value = -0.02, step = 0.01),
      numericInput("endoLog2FCMax", "Maximum log2FoldChange (default: 0.02)", value = 0.02, step = 0.01),
      hr(),

      h4("Filtering Parameters: DEG"),
      numericInput("degPvalue", "Maximum p-value threshold (default: 0.05)", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("degBaseMean", "Minimum baseMean value (default: 500)", value = 500, min = 0),
      numericInput("degLog2FC", "Minimum absolute log2FoldChange (default: 0.2)", value = 0.2, min = 0, step = 0.01),
      hr(),

      # Plotting options
      h4("Plotting Options"),
      selectInput("plotGeneType", "Select Gene Type", choices = c("DEG", "Endogenous Control")),
      uiOutput("geneSelector"),  # Dynamically generated gene selector
      actionButton("plotGene", "Plot Gene Expression"),
      hr(),

      h4("Relative Gene Expression"),
      uiOutput("endoControlSelector"),
      uiOutput("degSelector"),
      actionButton("plotRelative", "Plot Relative Expression"),
      hr(),

      # Download buttons
      downloadButton("downloadEndogenous", "Download Endogenous Control Table"),
      downloadButton("downloadDeg", "Download DEG Table"),
      hr(),
      downloadButton("downloadGenePlot", "Download Gene Expression Plot (PDF)"),
      downloadButton("downloadRelativePlot", "Download Relative Expression Plot (PDF)")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Endogenous Control Genes", tableOutput("endogenousTable")),
        tabPanel("DEG Genes", tableOutput("degTable")),
        tabPanel("Gene Expression Plot", plotOutput("genePlot")),
        tabPanel("Relative Expression Plot", plotOutput("relativePlot"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  # Data processing logic
  dataInput <- eventReactive(input$analyze, {
    req(input$countsFile, input$groupFile)

    counts_matrix <- read.csv(input$countsFile$datapath, row.names = 1)
    sample_info <- read.table(input$groupFile$datapath, header = TRUE, sep = "\t", row.names = 1)

    if (!all(rownames(sample_info) %in% colnames(counts_matrix))) {
      showNotification("The sample order in the group information file does not match the column order in the count matrix. Please check your files.", type = "error")
      return(NULL)
    }

    counts_matrix <- counts_matrix[, rownames(sample_info)]
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = sample_info, design = ~ condition)
    dds <- dds[rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds)

    normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
    raw_counts <- as.data.frame(counts(dds))

    res <- results(dds)
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

    endogenous_control_data <- res_df %>% filter(endogenous_control == TRUE)
    deg_data <- res_df %>% filter(DEG == TRUE)

    list(
      raw_counts = raw_counts,
      normalized_counts = normalized_counts,
      sample_info = sample_info,
      endogenous_control_data = endogenous_control_data,
      deg_data = deg_data
    )
  })

  # Dynamic gene selectors
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
    selectInput("selectedEndoControl", "Select Endogenous Control Gene", choices = data$endogenous_control_data$gene_name)
  })

  output$degSelector <- renderUI({
    data <- dataInput()
    req(data)
    selectInput("selectedDeg", "Select Target Gene (DEG)", choices = data$deg_data$gene_name)
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

  # Download tables
  output$downloadEndogenous <- downloadHandler(
    filename = function() { paste("endogenous_control_table_", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      data <- dataInput()
      req(data)
      write.csv(data$endogenous_control_data %>% select(gene_name, baseMean, log2FoldChange, pvalue, padj), file, row.names = FALSE)
    }
  )

  output$downloadDeg <- downloadHandler(
    filename = function() { paste("deg_table_", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      data <- dataInput()
      req(data)
      write.csv(data$deg_data %>% select(gene_name, baseMean, log2FoldChange, pvalue, padj), file, row.names = FALSE)
    }
  )

  # Plotting logic: Gene Expression Plot
  genePlotData <- reactive({
    data <- dataInput()
    req(data)
    normalized_counts <- data$normalized_counts
    sample_info <- data$sample_info
    selected_gene <- input$selectedGene
    req(selected_gene %in% rownames(normalized_counts))

    data.frame(
      Sample = colnames(normalized_counts),
      Expression = as.numeric(normalized_counts[selected_gene, ]),
      Group = sample_info$condition
    )
  })

  output$genePlot <- renderPlot({
    plot_data <- genePlotData()

    # Compute Wilcoxon p-value
    wilcox_p <- wilcox.test(Expression ~ Group, data = plot_data)$p.value

    # Generate the plot
    ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(
        title = paste0("Gene Expression: ", input$selectedGene),
        subtitle = paste("Wilcoxon p-value:", signif(wilcox_p, 3)),
        x = "Group",
        y = "Normalized Expression"
      ) +
      theme_minimal()
  })

  output$downloadGenePlot <- downloadHandler(
    filename = function() {
      paste("gene_expression_plot_", input$selectedGene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 7)
      print(ggplot(genePlotData(), aes(x = Group, y = Expression, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        labs(
          title = paste0("Gene Expression: ", input$selectedGene),
          x = "Group",
          y = "Normalized Expression"
        ) +
        theme_minimal())
      dev.off()
    }
  )

  # Plotting logic: Relative Expression Plot
  relativePlotData <- reactive({
    data <- dataInput()
    req(data)
    raw_counts <- data$raw_counts
    sample_info <- data$sample_info
    endo_gene <- input$selectedEndoControl
    deg_gene <- input$selectedDeg
    req(all(c(endo_gene, deg_gene) %in% rownames(raw_counts)))

    relative_expression <- raw_counts[deg_gene, ] / raw_counts[endo_gene, ]
    data.frame(
      Sample = colnames(raw_counts),
      RelativeExpression = as.numeric(relative_expression),
      Group = sample_info$condition
    )
  })

  output$relativePlot <- renderPlot({
    plot_data <- relativePlotData()

    # Compute Wilcoxon p-value
    wilcox_p <- wilcox.test(RelativeExpression ~ Group, data = plot_data)$p.value

    # Generate the plot
    ggplot(plot_data, aes(x = Group, y = RelativeExpression, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      labs(
        title = paste0("Relative Expression: ", input$selectedDeg, " / ", input$selectedEndoControl),
        subtitle = paste("Wilcoxon p-value:", signif(wilcox_p, 3)),
        x = "Group",
        y = "Relative Expression"
      ) +
      theme_minimal()
  })

  output$downloadRelativePlot <- downloadHandler(
    filename = function() {
      paste("relative_expression_plot_", input$selectedDeg, "_", input$selectedEndoControl, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 7)
      print(ggplot(relativePlotData(), aes(x = Group, y = RelativeExpression, fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        labs(
          title = paste0("Relative Expression: ", input$selectedDeg, " / ", input$selectedEndoControl),
          x = "Group",
          y = "Relative Expression"
        ) +
        theme_minimal())
      dev.off()
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
