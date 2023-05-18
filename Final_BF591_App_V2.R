# Importing required libraries/packages
library(shiny)  # for building interactive web apps
library(shinydashboard)  # for creating dashboard layouts
library(DT)  # for interactive tables
library(ggplot2)  # for creating plots
library(dplyr)  # for data manipulation
library(colourpicker)  # for interactive color picking
library(pheatmap)  # for creating heatmap visualizations
library(gridExtra)  # for arranging multiple grid-based plots on a page
library(tidyr)  # for data tidying
library(forcats)  # for working with categorical variables
library(tibble)  # for creating and manipulating tibble data frames
# Increasing the size limit for file uploads to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)
# App 1 - Sample Information Exploration
# App 1 UI function
app1_ui <- function(id) {
  # The NS function generates a namespace function for input and output
  ns <- NS(id)
  
  # fluidPage creates a page with automatically adjusting layout
  fluidPage(
    titlePanel("Sample Information Exploration"),  # title of the page
    
    # sidebarLayout divides the UI into a sidebar and a main panel
    sidebarLayout(
      # sidebarPanel contains the input controls
      sidebarPanel(
        # fileInput creates a file upload control
        fileInput(ns("data_file"), "Upload CSV File",
                  multiple = FALSE,  # only one file can be uploaded
                  accept = c("text/csv",  # the types of files that can be uploaded
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        # uiOutput renders a reactive UI element
        uiOutput(ns("columnSelect")),
        # actionButton creates a clickable button
        actionButton(ns("plotButton"), "Plot")
      ),
      
      # mainPanel contains the output elements
      mainPanel(
        # tabsetPanel creates a tabbed panel to switch between outputs
        tabsetPanel(
          id = ns("tabs"),
          # tabPanel creates an individual tab
          tabPanel("Summary", tableOutput(ns("summaryTable"))),  # tableOutput renders a reactive table
          tabPanel("Data", DTOutput(ns("dataTable"))),  # DTOutput renders a DataTable element
          tabPanel("Plots", plotOutput(ns("histogram")))  # plotOutput renders a reactive plot
        )
      )
    )
  )
}

# App 2
app2_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    # App title
    titlePanel("Counts Matrix Exploration"),
    
    # Sidebar layout
    sidebarLayout(
      sidebarPanel(
        # File input
        fileInput(ns("file"), "Choose a CSV file", accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv")),
        
        # Slider input controls
        sliderInput(ns("varianceThreshold"), "Variance percentile threshold:", min = 0, max = 100, value = 10, step = 1),
        sliderInput(ns("nonZeroThreshold"), "Non-zero samples threshold:", min = 0, max = 100, value = 10, step = 1),
        checkboxInput(ns("logTransform"), "Log-transform counts for visualization", value = FALSE)
        
      ),
      
      # Main panel
      mainPanel(
        tabsetPanel(
          tabPanel("Summary", tableOutput(ns("summaryTable"))),
          tabPanel("Scatter Plots", plotOutput(ns("scatterPlots"))),
          tabPanel("Clustered Heatmap", plotOutput(ns("heatmap"))),
          tabPanel("PCA", plotOutput(ns("pca")))
        )
      )
    )
  )
}


# App 3 - DESeq2 Results Visualization
app3_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("DESeq2 Results Visualization"),
    helpText("Upload your CSV file containing DESeq2 results to visualize the data as a volcano plot and a table."),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("file"), "Choose CSV File", accept = ".csv"),
        radioButtons(ns("x_axis"), "X-Axis",
                     choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                     selected = "log2FoldChange"),
        radioButtons(ns("y_axis"), "Y-Axis",
                     choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                     selected = "padj"),
        colourInput(ns("color1"), "Base point color", "green"),
        colourInput(ns("color2"), "Highlight point color", "blue"),
        sliderInput(ns("p_slider"), "P-Value Filter Magnitude", min = -300, max = 0, value = -150),
        actionButton(ns("plot_btn"), "Plot")
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Plot", plotOutput(ns("volcano"))),
                    tabPanel("Table", tableOutput(ns("table")))
        ))
    )
  )
}




# App 4 - Individual Gene Expressions Visualization
app4_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Visualization of Individual Gene Expressions"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("countsFile"), "Upload Normalized Counts Matrix CSV",
                  multiple = FALSE,
                  accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv")),
        fileInput(ns("sampleInfoFile"), "Upload Sample Information Matrix CSV",
                  multiple = FALSE,
                  accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv")),
        uiOutput(ns("categoricalFieldSelect")),
        uiOutput(ns("geneSelect")),
        selectInput(ns("plotType"), "Select Plot Type", choices = c("Bar Plot", "Box Plot", "Violin Plot", "Beeswarm Plot")),
        actionButton(ns("plotButton"), "Plot")
      ),
      mainPanel(
        plotOutput(ns("geneExpressionPlot"))
      )
    )
  )
}


# Combine the UIs into a single Shiny app
# The navbarPage function creates a page with a top level navigation bar
ui <- navbarPage("Bioinformatics Processes App",
                 # tabPanel creates an individual tab
                 tabPanel("Sample Info", app1_ui("app1")),
                 tabPanel("Counts Matrix", app2_ui("app2")),
                 tabPanel("DESeq2 Results", app3_ui("app3")),
                 tabPanel("Gene Expressions", app4_ui("app4")),
                 # footer adds a footnote to the page
                 footer = "This app allows users to explore and visualize different bioinformatics processes, including sample information exploration, counts matrix exploration, DESeq2 results visualization, and individual gene expressions visualization.")
                 
# Server function for App 1
app1_server <- function(id) {
  # Using moduleServer to isolate input and output
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive expression to read the uploaded CSV file
    data <- reactive({
      # Ensuring file is uploaded before attempting to read
      req(input$data_file)
      
      # Reading the CSV file
      df <- read.csv(input$data_file$datapath, stringsAsFactors = FALSE, header = TRUE)
      
      # Renaming the first column as "gene"
      colnames(df)[1] <- "gene"
      
      # Return the data frame
      df
    })
    
    # Render the summary table
    output$summaryTable <- renderTable({
      # Accessing the reactive 'data'
      df <- data()
      
      # Creating the summary data frame
      summary_df <- data.frame(Column = colnames(df),
                               Type = sapply(df, class),
                               Mean_SD = sapply(df, function(x) if (is.numeric(x)) paste0(round(mean(x), 2), " (+/- ", round(sd(x), 2), ")") else "NA"),
                               Distinct_Values = sapply(df, function(x) if (is.factor(x)) toString(unique(x)) else "NA"))
      
      # Return the summary data frame
      summary_df
    }, rownames = FALSE)
    
    # Render the data table
    output$dataTable <- renderDT({
      # Create a datatable from the reactive 'data'
      datatable(data())
    })
    
    # Reactive expression for numeric columns
    columns <- reactive({
      colnames(data()[, sapply(data(), is.numeric)])
    })
    
    # UI output for column select input
    output$columnSelect <- renderUI({
      selectInput(ns("column"), "Select a Column", choices = columns())
    })
    
    # Render the histogram plot
    output$histogram <- renderPlot({
      # Ensure the column input is selected
      req(input$column)
      
      # Accessing the reactive 'data'
      df <- data()
      
      # Selected column for the plot
      plot_column <- input$column
      
      # Check if a column is selected
      if (!is.null(plot_column) && plot_column != "") {
        # Plotting with selected column
        ggplot(df, aes(x = !!sym(plot_column))) +
          geom_histogram(bins = 30, fill = "steelblue", color = "black") +
          labs(x = plot_column, y = "Count", title = "Histogram")
      } else {
        # Default plot when no column is selected
        ggplot() +
          geom_histogram(fill = "steelblue", color = "black") +
          labs(x = "Data", y = "Count", title = "Histogram")
      }
    })
  })
}
# Define app2_server module
app2_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive expression to read the uploaded CSV file
    data <- reactive({
      req(input$file)
      read.csv(input$file$datapath, header = TRUE, row.names = 1)
    })
    
    # Reactive expression for filtered data based on variance and non-zero counts
    filteredData <- reactive({
      req(data())
      df <- data()
      
      # Calculate variances, medians and non-zero counts for each gene (row)
      variances <- apply(df, 1, var)
      medians <- apply(df, 1, median)
      non_zero_counts <- rowSums(df > 0)
      
      # Calculate thresholds for variance and non-zero counts
      variance_threshold <- quantile(variances, input$varianceThreshold / 100)
      non_zero_threshold <- ceiling(input$nonZeroThreshold * ncol(df) / 100)
      
      # Filter genes based on the thresholds
      filtered_df <- df[variances >= variance_threshold & non_zero_counts >= non_zero_threshold, ]
      
      filtered_df
    })
    
    # Render the summary table
    output$summaryTable <- renderTable({
      req(data())
      req(filteredData())
      df <- data()
      filtered_df <- filteredData()
      
      # Calculate summaries
      total_samples <- ncol(df)
      total_genes <- nrow(df)
      passing_genes <- nrow(filtered_df)
      not_passing_genes <- total_genes - passing_genes
      
      passing_percentage <- round((passing_genes / total_genes) * 100, 2)
      not_passing_percentage <- round((not_passing_genes / total_genes) * 100, 2)
      
      # Create summary table
      summary_df <- data.frame(
        Category = c("Total Samples", "Total Genes", "Passing Genes", "Not Passing Genes"),
        Count = c(total_samples, total_genes, passing_genes, not_passing_genes),
        Percentage = c("", "100%", paste0(passing_percentage, "%"), paste0(not_passing_percentage, "%"))
      )
      
      summary_df
    }, rownames = FALSE)
    
    # Render scatter plots
    output$scatterPlots <- renderPlot({
      req(data())
      req(filteredData())
      df <- data()
      filtered_df <- filteredData()
      
      # Calculate variances, medians, and zero counts
      variances <- apply(df, 1, var)
      medians <- apply(df, 1, median)
      zero_counts <- rowSums(df == 0)
      
      # Create a data frame for plotting
      plot_df <- data.frame(
        Gene = rownames(df),
        Variance = variances,
        Median = medians,
        ZeroCount = zero_counts,
        Passing = rownames(df) %in% rownames(filtered_df)
      )
      
      # Create scatter plots
      scatter_plot1 <- ggplot(plot_df) +
        geom_point(aes(x = log10(Median), y = log10(Variance), color = Passing), alpha = 0.6) +
        scale_color_manual(values = c("lightblue", "blue")) +
        labs(x = "Log10(Median Count)", y = "Log10(Variance)", title = "Median Count vs Variance") +
        theme_minimal()
      
      scatter_plot2 <- ggplot(plot_df) +
        geom_point(aes(x = log10(Median), y = ZeroCount,color = Passing), alpha = 0.6) +
        scale_color_manual(values = c("lightblue", "blue")) +
        labs(x = "Log10(Median Count)", y = "Number of Zeros", title = "Median Count vs Number of Zeros") +
        theme_minimal()

      # Combine the two scatter plots in a grid
      grid.arrange(scatter_plot1, scatter_plot2, ncol = 2)
    })

    # Render heatmap of the filtered data
    output$heatmap <- renderPlot({
      req(filteredData())
      filtered_df <- filteredData()

      # If the logTransform option is selected, apply the log2 transformation to the data
      if (input$logTransform) {
        heatmap_data <- log2(filtered_df + 1)
      } else {
        heatmap_data <- filtered_df
      }

      # Create the heatmap using the pheatmap function
      pheatmap(
        heatmap_data,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        show_colnames = TRUE,
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        legend = TRUE
      )
    })

    # Render PCA plot of the filtered data
    output$pca <- renderPlot({
      req(filteredData())
      filtered_df <- filteredData()

      # Perform PCA on the transposed filtered data
      pca <- prcomp(t(filtered_df), scale = TRUE)

      # Create a data frame from the PCA results for plotting
      pca_df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        Sample = rownames(pca$x),
        Color = factor(rownames(pca$x))
      )

      # Calculate the proportion of variance explained by the first two principal components
      explained_variance <- round((pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100, 2)

      # Create the PCA scatter plot
      pca_plot <- ggplot(pca_df) +
        geom_point(aes(x = PC1, y = PC2, color = Color), alpha = 0.6) +
        labs(
          x = paste0("PC1 (", explained_variance[1], "%)"),
          y = paste0("PC2 (", explained_variance[2], "%)"),
          title = "PCA Scatter Plot"
        ) +
        theme_minimal() +
        theme(legend.title = element_blank())

      pca_plot
    })
  })
}     

# Define app3_server module
# Define app3_server module
app3_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    load_data <- reactive({
      req(input$file)
      data <- read.csv(input$file$datapath)
      colnames(data)[1] <- "gene"
      return(data)
    })
    
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
      ggplot(dataf, aes_string(x = sym(x_name), y = paste("neg_log10_", y_name, sep = ""))) +
        geom_point(data = dataf, aes_string(color = paste("neg_log10_", y_name, sep = "")), size = 0.5, na.rm = TRUE) +
        scale_color_gradientn(colors = c(color1, color2)) +
        labs(title = "Volcano Plot", x = x_name, y = paste("-log10(", y_name,")")) +
        theme_minimal()
    }
    
    draw_table <- function(dataf, slider) {
      filtered_dataf <- dataf[which(dataf$padj < 10^slider), ]
      filtered_dataf$pvalue <- formatC(filtered_dataf$pvalue, digits = 3, format = "e")
      filtered_dataf$padj <- formatC(filtered_dataf$padj, digits = 3, format = "e")
      return(filtered_dataf)
    }
    
    observeEvent(input$plot_btn, {
      data <- load_data()
      if (!is.null(data)) {
        output$volcano <- renderPlot({
          data <- within(data, {neg_log10_padj <- -log10(padj)})
          if (is.null(data)) {
            return()
          }
          volcano_plot(data, input$x_axis, input$y_axis, input$p_slider, input$color1, input$color2)
        })
        
        output$table <- renderTable({
          data <- load_data()
          if (is.null(data)) {
            return()
          }
          draw_table(data, input$p_slider)
        }, digits = 3)
      }
    })
  })
}


# Define app4_server module
app4_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Read data from files
    countsData <- reactive({
      req(input$countsFile)
      inFile <- input$countsFile
      read.csv(inFile$datapath, row.names = 1)
    })
    
    sampleInfoData <- reactive({
      req(input$sampleInfoFile)
      inFile <- input$sampleInfoFile
      read.csv(inFile$datapath, header=TRUE)
      
    })
    
    # Update categorical field selection
    output$categoricalFieldSelect <- renderUI({
      req(sampleInfoData())
      selectInput(ns("categoricalField"), "Select Categorical Field", choices = colnames(sampleInfoData()))
    })
    
    # Update gene selection
    output$geneSelect <- renderUI({
      req(countsData())
      selectizeInput(ns("selectedGene"), "Select a Gene", choices = rownames(countsData()), multiple = FALSE, options = list(maxOptions = length(rownames(countsData()))))
    })
    
    # Create the plot
    output$geneExpressionPlot <- renderPlot({
      req(input$plotButton)
      input$plotButton
      
      req(countsData(), sampleInfoData(), input$selectedGene, input$categoricalField)
      
      gene_counts <- countsData()[input$selectedGene, , drop = FALSE]
      sample_info <- sampleInfoData()[, input$categoricalField, drop = FALSE]
      
      
      cat("Debugging gene_counts:\n")
      print(gene_counts)
      cat("Debugging sample_info:\n")
      print(sample_info)
      
      df <- data.frame(Sample = colnames(gene_counts), GeneCount = as.numeric(gene_counts), Group = sample_info[colnames(gene_counts), , drop = FALSE])
      
      
      cat("Debugging df:\n")
      print(df)
      
      plotType <- input$plotType
      print(head(df))
      
      gg <- ggplot(df, aes(x = Group, y = GeneCount)) +
        theme_bw() +
        labs(title = paste("Gene Expression for", input$selectedGene), x = "Group", y = "Normalized Gene Counts")
      print(gg)
      cat("Debugging data frame used for plot:\n")
      print(df)
      if (plotType == "boxplot") {
        gg <- gg + geom_boxplot()
      } else if (plotType == "violinplot") {
        gg <- gg + geom_violin() + geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.5)
      } else if (plotType == "barplot") {
        gg <- gg + stat_summary(fun.data = mean_se, geom = "bar", fill = "lightblue") +
          stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
      }
      
      plot(gg)
    })
    
    # Helper function for error bars in bar plot
    mean_se <- function(x) {
      data.frame(y = mean(x, na.rm = TRUE), ymin = mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE) / sqrt(length(x)), ymax = mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE) / sqrt(length(x)))
    }
  })
}


server <- function(input, output, session) {
  # Call server modules for each app
  app1_server("app1")
  app2_server("app2")
  app3_server("app3")
  app4_server("app4")
}

shinyApp(ui, server)
