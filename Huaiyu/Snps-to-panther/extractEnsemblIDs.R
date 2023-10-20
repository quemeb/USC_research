# Load necessary libraries
library(stringr)

# Function to extract unique Ensembl Gene IDs
extract_unique_ensembl_ids <- function(df) {
  # Create a regex pattern to capture Ensembl gene IDs
  pattern <- "ENSG[0-9]{10,12}" # Assuming there are 10-12 digits following "ENSG"
  
  # Filter columns with "ensembl" in their name
  ensembl_cols <- df[, grepl("ensembl", colnames(df), ignore.case = TRUE)]
  
  # Extract all matches
  matches <- unlist(lapply(ensembl_cols, function(col) {
    str_extract_all(col, pattern)
  }))
  
  # Remove any NA values and find unique matches
  unique_ensembl_matches <- unique(matches[!is.na(matches)])
  
  return(unique_ensembl_matches)
}