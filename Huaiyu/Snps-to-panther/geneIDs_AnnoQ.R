if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("AnnoQR", quietly = TRUE))
  devtools::install_github("USCbiostats/AnnoQR")

# Assuming rs_ID_clean$refsnp_id is a list/vector of rsid values

# Desired fields
desired_fields <- c(
  "ANNOVAR_ensembl_Gene_ID",
  "ANNOVAR_ensembl_Closest_gene(intergenic_only)",
  "SnpEff_ensembl_Gene_ID",
  "VEP_ensembl_Gene_ID",
  "ANNOVAR_refseq_Gene_ID",
  "ANNOVAR_refseq_Closest_gene(intergenic_only)",
  "SnpEff_refseq_Gene_ID",
  "VEP_refseq_Gene_ID"
)

# Function to retrieve desired fields from result
extract_desired_fields <- function(result) {
  sapply(desired_fields, function(field) {
    if (!is.null(result[[field]])) {
      return(result[[field]])
    } else {
      return(NA)
    }
  })
}

# Function to query ANNOQR and retrieve results
get_annoq_annotations <- function(rsid_list) {
  results <- lapply(rsid_list, function(rsid) {
    query_result <- rsidQuery(rsid)
    
    # Check if there are hits
    if (query_result$hits$total$value > 0) {
      # Take the first hit
      hit <- query_result$hits$hits[[1]]$`_source`
      return(extract_desired_fields(hit))
    } else {
      return(rep(NA, length(desired_fields)))
    }
  })
  
  # Convert results to a data.frame
  df_results <- as.data.frame(do.call(rbind, results), stringsAsFactors = FALSE)
  colnames(df_results) <- desired_fields
  
  return(df_results)
}



