# Ensure necessary packages are installed
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

# Load required libraries
library(httr)
library(jsonlite)

# Function to retrieve SNP information from Ensembl API for the specified genome version
rs_ID_get <- function(snp_positions, genome_version) {
  base_url <- ifelse(genome_version == "37", 
                     "https://grch37.rest.ensembl.org", 
                     "https://rest.ensembl.org") # Static URL for GRCh38
  endpoint <- "/overlap/region/human/"
  
  results <- lapply(snp_positions, function(pos) {
    url <- paste0(base_url, endpoint, pos, "?feature=variation")
    response <- GET(url, accept("application/json"))
    stop_for_status(response)
    content(response, as = "parsed", type = "application/json")
  })
  
  # Combine results into a single data frame and ensure unique entries
  results <- do.call(rbind, lapply(results, function(x) {
    if (length(x) > 0) {
      data.frame(do.call(rbind, x), stringsAsFactors = FALSE)
    } else {
      NULL
    }
  }))
  
  return(results)
}

# Function to match rsID back to df
match_rsID_to_df <- function(df, rs_ID) {
  df$rsID <- NA # Initialize the new rsID column
  
  # Loop through each row in df
  for (i in 1:nrow(df)) {
    start_pos <- df$start[i]
    
    # Loop through each row in rs_ID
    for (j in 1:nrow(rs_ID)) {
      row_values <- unlist(rs_ID[j, ])
      
      # Check if start_pos appears twice in the row
      if (sum(row_values == start_pos) == 2) {
        # Look for the cell that starts with "rs"
        rsid <- row_values[grep("^rs[0-9]{7}", row_values)]
        
        if (length(rsid) > 0) {
          df$rsID[i] <- rsid[1]
        }
      }
    }
  }
  
  return(df)
}
