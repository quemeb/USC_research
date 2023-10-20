#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File"),
      selectInput("version", "Choose Genome Version", choices = c("37", "38")),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      # Define tabs
      tabsetPanel(
        tabPanel("rsIDs", 
                 fluidRow(
                   column(8),  # Empty column to push the button to the right
                   column(4, align = "right", downloadButton("downloadRsTable", "Download rsIDs"))
                 ),
                 tableOutput("rsTable")
        ),
        tabPanel("Results", 
                 fluidRow(
                   column(8),  
                   column(4, align = "right", downloadButton("downloadResults", "Download Results"))
                 ),
                 tableOutput("resultsTable")
        ),
        tabPanel("Ensembl IDs", 
                 fluidRow(
                   column(8),  
                   column(4, align = "right", downloadButton("downloadEnsembl", "Download Ensembl IDs"))
                 ),
                 tableOutput("ensemblTable")
        )
      )
    )
  )
)

server <- function(input, output) {
  
  # This will store the rs_ID_clean dataframe and annotations
  rsData <- reactiveVal()
  resultsData <- reactiveVal()
  ensemblData <- reactiveVal()
  
  observeEvent(input$submit, {
    inFile <- input$file
    genomeVersion <- input$version
    
    if (!is.null(inFile)) {
      # Set the variables file_path and genome_version 
      file_path <<- inFile$datapath
      genome_version <<- genomeVersion
      
      # Display a modal informing the user that processing has started
      showModal(modalDialog(
        title = "Processing",
        "File received. Loading results... Please wait.",
        footer = NULL
      ))
      
      # Now source the rsid_Ensembl.R script
      source("C:/Users/bryan/Desktop/USC_research/Huaiyu/Snps-to-panther/rsid_Ensembl.R")
      rsData(rs_ID_clean$refsnp_id)
      
      # Get ANNOQR annotations
      source("C:/Users/bryan/Desktop/USC_research/Huaiyu/Snps-to-panther/geneIDs_AnnoQ.R")
      df_results <- get_annoq_annotations(rs_ID_clean$refsnp_id)
      resultsData(df_results) # Assuming df_results is what the geneIDs_AnnoQ.R script outputs
      
      # Fetch unique Ensembl IDs
      source("C:/Users/bryan/Desktop/USC_research/Huaiyu/Snps-to-panther/extractEnsemblIDs.R")
      unique_ensembl_matches <- extract_unique_ensembl_ids(df_results)
      ensemblData(unique_ensembl_matches)
      
      # Remove the modal after processing
      removeModal()
    }
  })
  
  # Render the rsID table
  output$rsTable <- renderTable({data.frame(rsID = rsData())})
  
  # Render the annotations table
  output$resultsTable <- renderTable({resultsData()})
  
  # Render the Ensembl ID table
  output$ensemblTable <- renderTable({data.frame(Ensembl_ID = ensemblData())})
  
######### DOWNLOAD ##########
  
  # Download handler for rsTable
  output$downloadRsTable <- downloadHandler(
    filename = function() {
      paste("rsTable-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(data.frame(rsID = rsData()), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  # Download handler for resultsTable
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("Results-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(resultsData(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  # Download handler for ensemblTable
  output$downloadEnsembl <- downloadHandler(
    filename = function() {
      paste("Ensembl_IDs-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(data.frame(Ensembl_ID = ensemblData()), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
