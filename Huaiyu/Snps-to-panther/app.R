#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(shiny)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File"),
      selectInput("version", "Choose Genome Version", choices = c("37", "38")),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      tableOutput("rsTable")
    )
  )
)
server <- function(input, output) {
  # This will store the rs_ID_clean dataframe
  rsData <- reactiveVal()
  
  observeEvent(input$submit, {
    inFile <- input$file
    genomeVersion <- input$version
    
    if (!is.null(inFile)) {
      # Set the variables file_path and genome_version 
      file_path <<- inFile$datapath
      genome_version <<- genomeVersion
      
      # Now source the rsid_Ensembl.R script
      source("C:/Users/bryan/Desktop/USC_research/Huaiyu/Snps-to-panther/rsid_Ensembl.R")
      
      # Assign the rs_ID_clean$refsnp_id to rsData
      rsData(rs_ID_clean$refsnp_id)
    }
  })
  
  # Render the table using the data stored in rsData
  output$rsTable <- renderTable({
    data.frame(rsID = rsData())
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
