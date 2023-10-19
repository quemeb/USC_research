if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(biomaRt)
library(stringr)
library(dplyr)

# Set up Ensembl
ensembl = useEnsembl(biomart="snps", dataset = "hsapiens_snp", version = "GRCh37")

# Function to get information from Ensembl
rs_ID_get <- function(snp_positions){
  getBM(attributes = c("refsnp_id", "allele", "chr_name", "chrom_start", 
                       "allele_1", "minor_allele", "minor_allele_freq", 
                       "clinical_significance"),
        filters = "chromosomal_region",
        values = snp_positions,
        mart = ensembl)
}

# Load file
#file_path <- "/cloud/project/SNPs_Seq37_example.csv"
df <- read.csv(file_path)

# Prepare data
snp_pos <- apply(df, 1, paste, collapse = ":")
rs_ID <- rs_ID_get(snp_pos)

# Filter data directly (replaces the loop)
rs_ID_clean <- rs_ID %>% 
  filter(str_starts(refsnp_id, "rs"), !is.na(minor_allele_freq))

# Check for missing chrom_start
missing_start <- setdiff(df$chr_start, rs_ID_clean$chrom_start)

# If any are missing, add them to rs_ID_clean
if (length(missing_start) > 0) {
  missing_data <- rs_ID %>% 
    filter(chrom_start %in% missing_start, str_starts(allele, "[:upper:]/[:upper:]"), str_starts(refsnp_id, "rs"))
  
  rs_ID_clean <- rbind(rs_ID_clean, missing_data)
}

# The rs_ID_clean should now be ready for use.
