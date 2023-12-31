AnnoQR2.0
================
Bryan Queme
2023-09-08

``` r
# Try to load AnnoQR
if (!requireNamespace("AnnoQR", quietly = TRUE)) {
  # Check for devtools and install it if necessary
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  # Load devtools
  library(devtools)
  # Install AnnoQR from GitHub
  install_github("USCbiostats/AnnoQR")
}

# Load AnnoQR
library(AnnoQR)
```

``` r
# Only using the basic annotation types
basic_annotations = c("chr", "pos", "ref", "alt",
                      "ANNOVAR_ensembl_Effect",
                      "ANNOVAR_ensembl_Closest_gene(intergenic_only)",
                      "ANNOVAR_ensembl_summary",
                      "SnpEff_ensembl_Effect",
                      "SnpEff_ensembl_Effect_impact",
                      "SnpEff_ensembl_Gene_name",
                      "SnpEff_ensembl_Gene_ID",
                      "SnpEff_ensembl_HGVSc",
                      "SnpEff_ensembl_summary",
                      "VEP_ensembl_Consequence",
                      "VEP_ensembl_summary",
                      "ANNOVAR_refseq_Effect",
                      "ANNOVAR_refseq_Closest_gene(intergenic_only)",
                      "ANNOVAR_refseq_summary",
                      "SnpEff_refseq_Effect",
                      "SnpEff_refseq_Effect_impact",
                      "SnpEff_refseq_Gene_name",
                      "SnpEff_refseq_Gene_ID",
                      "SnpEff_refseq_HGVSc",
                      "SnpEff_refseq_summary",
                      "VEP_refseq_Consequence",
                      "VEP_refseq_summary",
                      "ANNOVAR_ucsc_Effect",
                      "ANNOVAR_ucsc_Closest_gene(intergenic_only)",
                      "ANNOVAR_ucsc_summary",
                      "rs_dbSNP151")
```

``` r
hits <- function(variants) {
  # Finding number of hits
  len <- length(variants$hits$hits)
  
  # Early return if no annotations
  if (len == 0) {
    return(data.frame())
  }
  
  # Creating an empty data frame with pre-allocated size
  n <- length(basic_annotations)
  df <- data.frame(matrix(ncol = n, nrow = len))
  colnames(df) <- basic_annotations
  
  # Loop through hits and annotations
  for (i in seq_len(len)) {
    hit <- variants$hits$hits[[i]]$'_source'
    for (j in seq_along(basic_annotations)) {
      annotation <- basic_annotations[j]
      # Fill the data frame with annotation info
      df[i, j] <- ifelse(is.null(hit[[annotation]]), " ", hit[[annotation]])
    }
  }
  
  return(df)
}
```

``` r
#url <- choose.files()
url <- "https://raw.githubusercontent.com/quemeb/USC_research/main/Kelly/ex_peaks_example.csv"
cdata <- read.csv(url)
```

``` r
# Function to retrieve hits for each region
get_hits_for_region <- function(chr, start, end) {
  variants = regionQuery(contig = chr, start = start, end = end)
  hits(variants)
}

# Function to fetch all hits
get_all_hits <- function(cdata, basic_annotations) {
  len_regions <- nrow(cdata)
  
  # List to store individual data frames
  list_df <- vector("list", len_regions)
  
  # Loop through each region
  for(i in seq_len(len_regions)) {
    list_df[[i]] <- get_hits_for_region(cdata$chr[i], cdata$start[i], cdata$end[i])
  }
  
  # Combine all data frames into one
  df_all <- do.call(rbind, list_df)
  
  # Add column names
  colnames(df_all) <- basic_annotations
  
  return(df_all)
}

df_all <- get_all_hits(cdata, basic_annotations)
```

``` r
write.table(df_all, "ex_peaks.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```
