library(shiny)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File"),
      selectInput("version", "Choose Genome Version", choices = c("37", "38")),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      fluidRow(
        column(8),  # Empty column to push the button to the right
        column(4, align = "right", downloadButton("downloadModifiedDf", "Download Modified Data"))
      ),
      tableOutput("finalTable")
    )
  )
)

server <- function(input, output) {
  
  finalData <- reactiveVal()
  
  observeEvent(input$submit, {
    inFile <- input$file
    genomeVersion <- input$version
    
    if (!is.null(inFile)) {
      file_path <<- inFile$datapath
      genome_version <<- genomeVersion
      
      showModal(modalDialog(
        title = "Processing",
        "File received. Loading results... Please wait.",
        footer = NULL
      ))
      
      source("C:\\Users\\bryan\\Desktop\\USC_research\\Huaiyu\\Snps-to-panther\\rsid_Ensembl.R")
      source("C:\\Users\\bryan\\Desktop\\USC_research\\Huaiyu\\Snps-to-panther\\geneIDs_AnnoQ.R")
      
      df <- read.csv(file_path)
      
      snp_pos <- apply(df, 1, function(row) {
        paste(row['chr'], row['start'], row['end'], sep = ":")
      })
      
      rs_ID <- rs_ID_get(snp_pos, genome_version)
      rs_ID <- as.data.frame(rs_ID)
      
      df <- match_rsID_to_df(df, rs_ID)
      
      annotations_list <- lapply(df$rsID, process_annotations)
      annotations_df <- do.call(rbind, annotations_list)
      
      final_df <- match_annotations_to_df(df, annotations_df)
      
      # Select and rename columns for display
      final_df <- final_df[, c("chr", "start", "end", "rsID", "Ensembl_Genes", "RefSeq_Genes")]
      colnames(final_df) <- c("Chromosome", "Start", "End", "rsID", "Ensembl Gene", "RefSeq Gene")
      
      finalData(final_df)
      
      removeModal()
    }
  })
  
  output$finalTable <- renderTable({
    finalData()
  })
  
  output$downloadModifiedDf <- downloadHandler(
    filename = function() {
      paste("Modified_Data-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      temp_file <- tempfile(fileext = ".txt")
      write.table(finalData(), temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
      file.copy(temp_file, file)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
