if (!require("httr", quietly = TRUE))
  install.packages("httr")

if (!require("jsonlite", quietly = TRUE))
  install.packages("jsonlite")

library(httr)
library(jsonlite)


# Define the Base URL
Annotations_URL <- "http://annoq.org/api-v2/graphql"

# Convert a list of fields into a GraphQL query format
create_annotations_query_string <- function(annotations) {
  paste(annotations, collapse = "\n")
}

# Perform a GraphQL query
perform_graphql_query <- function(query) {
  response <- POST(
    Annotations_URL, 
    content_type_json(), 
    body = list(query = query),
    encode = "json"
  )
  stop_for_status(response)
  content(response, "text", encoding = "UTF-8")
}


# Desired fields
annotations_to_retrieve <- c(
  "ANNOVAR_ensembl_Gene_ID",
  "ANNOVAR_ensembl_Closest_gene",
  "SnpEff_ensembl_Gene_ID",
  "VEP_ensembl_Gene_ID",
  "ANNOVAR_refseq_Gene_ID",
  "ANNOVAR_refseq_Closest_gene",
  "SnpEff_refseq_Gene_ID",
  "VEP_refseq_Gene_ID",
  "enhancer_linked_genes"
)


# Get SNP by rsID
rsidQuery <- function(rsID, annotations_to_retrieve) {
  annotations_query_string <- create_annotations_query_string(annotations_to_retrieve)
  query <- sprintf('
  query {
    get_SNPs_by_RsID(rsID: "%s", query_type_option: SNPS, filter_args: {exists: ["rs_dbSNP151"]}) {
      snps {
        %s
      }
    }
  }', rsID, annotations_query_string)
  
  response_content <- perform_graphql_query(query)
  data <- fromJSON(response_content, flatten = TRUE)
  data$data$get_SNPs_by_RsID$snps
}

# Function to process and extract gene annotations for a single rsID
process_annotations <- function(rsID) {
  query_result <- rsidQuery(rsID, annotations_to_retrieve)
  
  if (nrow(query_result) > 0) {
    # Combine Ensembl genes, excluding NAs
    ensembl_genes <- unlist(query_result[grep("ensembl", names(query_result))])
    ensembl_genes <- ensembl_genes[!is.na(ensembl_genes)]
    
    # Separate Ensembl genes by "-" and discard parts after ":"
    ensembl_genes <- unlist(strsplit(ensembl_genes, "[-:]"))
    ensembl_genes <- ensembl_genes[grepl("^ENSG", ensembl_genes)]
    ensembl_genes <- unique(ensembl_genes)
    
    # Combine RefSeq genes, excluding NAs
    refseq_genes <- unlist(query_result[grep("refseq", names(query_result))])
    refseq_genes <- refseq_genes[!is.na(refseq_genes)]
    
    # Separate RefSeq genes by "-" and discard parts after ":"
    refseq_genes <- unlist(strsplit(refseq_genes, "[-:]"))
    refseq_genes <- unique(refseq_genes)
    
    return(data.frame(rsID = rsID,
                      Ensembl_Genes = paste(ensembl_genes, collapse = ","),
                      RefSeq_Genes = paste(refseq_genes, collapse = ","),
                      stringsAsFactors = FALSE))
  } else {
    return(data.frame(rsID = rsID, Ensembl_Genes = NA, RefSeq_Genes = NA, stringsAsFactors = FALSE))
  }
}


# Function to match gene annotations back to the original df by rsID
match_annotations_to_df <- function(df, annotations_df) {
  # Merge the annotations with the original df on rsID
  merged_df <- merge(df, annotations_df, by = "rsID", all.x = TRUE)
  return(merged_df)
}





