# This example is adapted from
# Grolemund, G. (2015). Learn Shiny - Video Tutorials. URL:https://shiny.rstudio.com/tutorial/

library(shiny)
library(shinyalert)

# set max file upload limit
options(shiny.maxRequestSize = 10 * 1024^2)


# Define UI
ui <- fluidPage(

  # title page
  titlePanel("CellSearchR: automated annotation of single cell RNAseq data
            with reference mapping"),
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      br(),

      tags$b("Description: CellSearchR is an R package designed to automate the
      annotation of single cell RNA seq data by leveraging the Cellxgene
      database. CellSearchR is unique in that it is the first R package to use
      the Cellxgene census database for reference mapping and automatic cell
      annotation. "),

      br(),
      br(),

      # input
      tags$p("Instruction: Select the provided reference dataset and
      upload your query sc-RNAseq count matrix for annotation. After
      uploading, click run at the bottom to start reference mapping!"),

      br(),

      # input
      shinyalert::useShinyalert(force = TRUE),  # Set up shinyalert

      # selectable list for available  cell embedding method
      selectInput("embed", "Choose your cell embedding method",
                  list("",
                       "Seurat")),

      # selectable list for available reference dataset
      selectInput("reference", "Choose your reference dataset",
                  list("",
                    "Cellxgene3K")),

      # file input for query dataset
      fileInput("query", "Upload your query dataset: scRNAseq count
                matrix (Seurat object in .rda format)"),

      # instruction for demo file download
      tags$div(
        "You may find the demo query dataset at `\\data\\covid.pbmc3k.rda` or
      download from",
        tags$a(href="https://drive.google.com/file/d/1q7a9DtWTQhlNT1sKl8WyExiepQMEL4qF/view?usp=sharing",
               "here.", style = "color:blue")
      ),

      br(),
      br(),


      # actionButton
      actionButton(inputId = "actionButton",
                   label = "Run"),
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot annotation results",
                           div("Instructions: Enter values and click 'Run' at the bottom left side."),
                           br(),
                           plotOutput("RNAseqPlot")),
                  tabPanel("Summary of annotation results",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Summary of Information Criteria Values:"),
                           br(),
                           h4("Bayesian information criterion (BIC)"),
                           verbatimTextOutput("textOutBIC"),
                           h4("Integrated Complete Likelihood (ICL)"),
                           verbatimTextOutput("textOutICL"),
                           h4("Akaike Information Criterion (AIC)"),
                           verbatimTextOutput("textOutAIC")),
        )
      )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # selected embeddings
  embedMethod <- reactive({
    data <- switch(input$embed,
           "Seurat" = "Seurat")
    return(data)
  })

  # selected reference dataset
  refData <- reactive({
    data <- switch(input$reference,
           "Cellxgene3K" = CellSearchR::Cellxgene3K)
    return(data)
  })

  # upload query dataset
  queryData <- reactive({
    req(input$query)
    data <- load(input$query$datapath)
    return(data)
  })

  # run bottom logic
  observeEvent(input$actionButton, {
    withProgress(message = "Running reference mapping", {
      mapping_results <- CellSearchR::searchCell(queryData(),
                                                 refData(),
                                                 embedMethod())
      plot <- CellSearchR::plotCellSearchR(mapping_results)
    })
  })

  # output plot
  output$RNAseqPlot <- renderPlot({
    return(plot)
  })

}

# Create Shiny app ----
shiny::shinyApp(ui, server)

# [END]
