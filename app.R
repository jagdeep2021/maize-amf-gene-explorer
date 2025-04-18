# Load required packages
library(shiny)
library(shinyWidgets)
library(shinyjs)
library(tidyverse)
library(ggplot2)
library(DT)
library(bslib)

# Source the plot function
source("plot_function.R")

# Read the table data once at startup
table_data <- read_csv("data/DE_ASE_table.csv") %>%
  select(geneID, Working.Symbol, ASE_type, ASE_sign, Gen_type, Myc_gene)

# UI
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2c3e50",
    secondary = "#18bc9c",
    success = "#2ecc71",
    info = "#3498db",
    warning = "#f39c12",
    danger = "#e74c3c"
  ),
  useShinyjs(),  # Initialize shinyjs
  titlePanel(
    div(
      h1("Maize-AMF Gene Explorer", style = "color: #2c3e50; margin-bottom: 20px;"),
      p("Explore gene expression patterns and allele-specific expression in maize with AMF treatment", 
        style = "color: #7f8c8d; font-size: 16px;")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "background-color: #f8f9fa; padding: 20px; border-radius: 5px;",
      div(
        style = "margin-bottom: 20px;",
        h4("Input Parameters", style = "color: #2c3e50; margin-bottom: 15px;"),
        textAreaInput("gene_list", "Enter Gene List (one per line):", 
                     height = "200px",
                     placeholder = "Example:\nZm00001eb164860\nbx13"),
        div(
          style = "margin-top: 20px;",
          h5("Filters", style = "color: #2c3e50; margin-bottom: 10px;"),
          checkboxGroupInput("trt_filter", "Select Treatments:",
                           choices = c("M", "C"),
                           selected = c("M", "C"),
                           inline = TRUE),
          checkboxGroupInput("time_filter", "Select Time Points:",
                           choices = c("T1", "T2"),
                           selected = c("T1", "T2"),
                           inline = TRUE)
        ),
        div(
          style = "margin-top: 20px;",
          actionButton("plot_button", "Generate Plots", 
                      class = "btn-primary",
                      style = "width: 100%; margin-bottom: 10px;"),
          downloadButton("download_plots", "Download Plots",
                        class = "btn-success",
                        style = "width: 100%;")
        ),
        progressBar(id = "progress", value = 0, display_pct = TRUE)
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Results", 
                 div(
                   style = "margin-bottom: 20px;",
                   h4("Gene Information", style = "color: #2c3e50;"),
                   dataTableOutput("gene_table")
                 ),
                 div(
                   style = "margin-bottom: 20px;",
                   h4("Expression Patterns by Genotype", style = "color: #2c3e50;"),
                   plotOutput("expression_plot", height = "400px")
                 ),
                 div(
                   style = "margin-bottom: 20px;",
                   h4("Allele-Specific Expression Results", style = "color: #2c3e50;"),
                   plotOutput("ase_plot", height = "400px")
                 )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive value to store plot objects
  plots <- reactiveValues(
    de_expr = NULL,
    ase = NULL,
    table_data = NULL
  )
  
  # Reactive expression for genes
  genes <- reactive({
    req(input$gene_list)
    genes <- strsplit(input$gene_list, "\n")[[1]]
    genes[genes != ""]  # Remove empty lines
  })
  
  # Generate plots when button is clicked
  observeEvent(input$plot_button, {
    # Show loading message
    show("loading")
    
    # Update progress bar
    updateProgressBar(session, "progress", value = 0, title = "Starting plot generation...")
    
    # Get genes from reactive expression
    genes_list <- genes()
    
    if (length(genes_list) == 0) {
      showNotification("Please enter at least one gene", type = "error")
      hide("loading")
      return()
    }
    
    # Generate plots with error handling
    tryCatch({
      # Read and store table data
      plots$table_data <- read_csv("data/DE_ASE_table.csv") %>%
        select(geneID, Working.Symbol, ASE_type, ASE_sign, Gen_type, Myc_gene)
      
      # Update progress
      updateProgressBar(session, "progress", value = 20, title = "Generating expression plot...")
      
      # Generate expression plot
      expr_result <- tryCatch({
        plot_expression_and_amreads(
          query_genes = genes_list,
          title_prefix = "Expression Patterns by Genotype",
          trt_filter = input$trt_filter,
          time_filter = input$time_filter
        )
        TRUE
      }, error = function(e) {
        showNotification(paste("Error generating expression plot:", e$message), type = "error")
        FALSE
      })
      
      if (!expr_result) {
        hide("loading")
        return()
      }
      
      # Store expression plot path
      plots$de_expr <- file.path("plots", "Expression Patterns by Genotype_expression_plot.png")
      
      # Update progress
      updateProgressBar(session, "progress", value = 60, title = "Generating ASE plot...")
      
      # Generate ASE plot
      ase_result <- tryCatch({
        plot_ase_expression(
          query_genes = genes_list,
          title_prefix = "Allele-Specific Expression Results",
          trt_filter = input$trt_filter,
          time_filter = input$time_filter
        )
        TRUE
      }, error = function(e) {
        showNotification(paste("Error generating ASE plot:", e$message), type = "error")
        FALSE
      })
      
      if (!ase_result) {
        hide("loading")
        return()
      }
      
      # Store ASE plot path
      plots$ase <- file.path("plots", "Allele-Specific Expression Results_ASE_expression_plot.png")
      
      # Update progress
      updateProgressBar(session, "progress", value = 100, title = "Done!")
      
      # Hide loading message
      hide("loading")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      updateProgressBar(session, "progress", value = 0, title = "Error occurred")
      hide("loading")
    })
  })
  
  # Render plots with error handling
  output$expression_plot <- renderPlot({
    req(plots$de_expr)
    tryCatch({
      if (file.exists(plots$de_expr)) {
        img <- png::readPNG(plots$de_expr)
        grid::grid.raster(img)
      } else {
        print(paste("Plot file not found:", plots$de_expr))
        NULL
      }
    }, error = function(e) {
      print(paste("Error loading expression plot:", e$message))
      NULL
    })
  })
  
  output$ase_plot <- renderPlot({
    req(plots$ase)
    tryCatch({
      if (file.exists(plots$ase)) {
        img <- png::readPNG(plots$ase)
        grid::grid.raster(img)
      } else {
        print(paste("Plot file not found:", plots$ase))
        NULL
      }
    }, error = function(e) {
      print(paste("Error loading ASE plot:", e$message))
      NULL
    })
  })
  
  # Add table functionality
  output$gene_table <- renderDataTable({
    # Get the gene list
    selected_genes <- unlist(strsplit(input$gene_list, "\n"))
    selected_genes <- selected_genes[selected_genes != ""]  # Remove empty lines
    
    if (length(selected_genes) == 0) {
      return(NULL)
    }
    
    # Filter table for selected genes
    table_data %>%
      filter(geneID %in% selected_genes | Working.Symbol %in% selected_genes)
  }, options = list(
    pageLength = 10,
    scrollX = TRUE,
    autoWidth = TRUE,
    searching = TRUE
  ))
  
  # Download handler for all plots
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("plots.zip")
    },
    content = function(file) {
      tryCatch({
        # Create a temporary directory
        temp_dir <- tempdir()
        
        # Copy all plot files to the temporary directory
        plot_files <- c(
          plots$de_expr,
          plots$ase
        )
        
        # Check if all files exist
        if (!all(file.exists(plot_files))) {
          stop("Some plot files are missing. Please generate the plots first.")
        }
        
        # Copy files to temp directory
        file.copy(plot_files, temp_dir)
        
        # Create zip file
        zip_file <- file.path(temp_dir, "plots.zip")
        zip::zip(zip_file, 
                 files = basename(plot_files),
                 root = temp_dir)
        
        # Copy zip file to download
        file.copy(zip_file, file)
      }, error = function(e) {
        showNotification(paste("Error creating zip file:", e$message), type = "error")
      })
    }
  )
  
  # Clear plots when inputs change
  observe({
    input$gene_list
    input$trt_filter
    input$time_filter
    
    plots$de_expr <- NULL
    plots$ase <- NULL
  })
}

# Run the app
shinyApp(ui = ui, server = server)