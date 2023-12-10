# This example is adapted from
# Grolemund, G. (2015). Learn Shiny - Video Tutorials. URL:https://shiny.rstudio.com/tutorial/

library(shiny)
library(shinyalert)

# set max file upload limit
options(shiny.maxRequestSize = 100 * 1024^2)


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

      # selectable list for available reference dataset
      selectInput("reference", "Choose your reference dataset",
                  list("",
                    "Cellxgene3K")),

      # file input for query dataset
      fileInput("query", "Upload your query dataset: scRNAseq count
                matrix (.rda format)"),
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
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Pairs Plot of Log-transformed RNAseq Count Dataset:"),
                           br(),
                           # textOutput("refData"),


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

  # selected reference dataset
  refData <- reactive({
    switch(input$reference,
           "Cellxgene3K" = CellSearchR::Cellxgene3K)
  })

  # upload query dataset
  queryData <- reactive({
    req(input$query)
    load(input$query$datapath)
  })

  # run bottom logic

  mappingResults <- eventReactive(eventExpr = input$actionButton, {
    # CellSearchR::CellSearchR(
    #   refData = refData(),
    #   queryData = queryData()
    # )

    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
        incProgress(0.1, detail = paste("Doing part", "1"))
        Sys.sleep(1)
        incProgress(0.5, detail = paste("Doing part", "1"))
        Sys.sleep(100)

      })

  })


}

# Create Shiny app ----
shiny::shinyApp(ui, server)

# [END]
