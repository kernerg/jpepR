#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Read command line arguments
Trait1 <- as.character(args[1])  # Primary trait
aux_traits <- strsplit(args[2], ",")[[1]]  # Vector of auxiliary trait names
pip_threshold <- as.numeric(args[3])  # PIP threshold
aux_pip_thres <- as.numeric(args[4])
aux_nbr_thres <- as.numeric(args[5])
dir <- as.character(args[6]) # directory

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(glue)
})

# Set directories
dir_sqlite <- sprintf("%s/SQLite", dir)
dir_output <- sprintf("%s/V_traits", dir)
dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

# Read primary trait data
primary_file <- sprintf("%s/%s_SQLite.db", dir_sqlite, Trait1)
con <- dbConnect(SQLite(), primary_file)
primary_data <- as.data.table(dbGetQuery(con, "SELECT CHR, BP, RSID, A1, A2, PIP, BETA FROM trait_data"))
dbDisconnect(con)

# Filter by PIP threshold and remove duplicates
primary_data <- primary_data[primary_data$PIP > pip_threshold, ]
primary_data <- unique(primary_data, by = "RSID")

# Initialize the V_trait matrix
V_traits <- primary_data[, .(CHR, BP, RSID, A1, A2, PIP, BETA)]
setorder(V_traits, RSID)

# Process each auxiliary trait
aux_traits <- gsub("_SQLite.db", "", aux_traits)
aux_traits <- setdiff(aux_traits, Trait1)  # Exclude the primary trait

for (Trait2 in aux_traits) {
  cat("Processing:", Trait2, "\n")
  aux_file <- sprintf("%s/%s_SQLite.db", dir_sqlite, Trait2)
  con <- dbConnect(SQLite(), aux_file)
  aux_data <- as.data.table(dbGetQuery(con, "SELECT RSID, A1, A2, PIP, BETA FROM trait_data"))
  dbDisconnect(con)

  # Merge with primary data on RSID
  merged_data <- merge(V_traits, aux_data, by = "RSID", all.x = TRUE, suffixes = c("", "_aux"))

  # Align effect alleles (flip if necessary)
  flip_index <- !is.na(merged_data$A1_aux) & (merged_data$A1 != merged_data$A1_aux)
  merged_data[flip_index, BETA_aux := -BETA_aux]
  merged_data[flip_index, A1_aux := A1]
  merged_data[flip_index, A2_aux := A2]

  # Calculate signed PIP product
  merged_data[, paste0(Trait2, "_pos") := pmax(PIP * PIP_aux * sign(BETA * BETA_aux), 0)]
  merged_data[, paste0(Trait2, "_neg") := pmax(-PIP * PIP_aux * sign(BETA * BETA_aux), 0)]
  
  # Complete NA with zeros
  merged_data[is.na(get(paste0(Trait2, "_pos"))), paste0(Trait2, "_pos") := 0]
  merged_data[is.na(get(paste0(Trait2, "_neg"))), paste0(Trait2, "_neg") := 0]

  # Add to V_trait
  V_traits <- cbind(V_traits, as.matrix(merged_data[, .SD, .SDcols = patterns(paste0(Trait2, "_"))]))
}

# Remove rows with all zero columns (excluding RSID, A1, A2, PIP)
aux_trait_columns <- setdiff(colnames(V_traits), c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA"))
non_zero_rows <- rowSums(V_traits[, ..aux_trait_columns, with = FALSE]) > 0
V_traits <- V_traits[non_zero_rows, ]

# Filter auxiliary traits
# Keep auxiliary traits with aux_nbr_thres SNPs with PIP > aux_pip_thres
filt1 <- V_traits[, lapply(.SD, function(x) sum(x > aux_pip_thres)), .SDcols = aux_trait_columns]
aux_trait_columns_filtered <- names(filt1)[which(filt1 > aux_nbr_thres)]
new_columns <- c(c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA"), aux_trait_columns_filtered)
V_traits <- V_traits[, ..new_columns]

# Save the final V_trait matrix
output_file <- sprintf("%s/%s_V_traits_thres%s.tsv.gz", dir_output, Trait1, pip_threshold)
fwrite(V_traits, file = output_file, sep = "\t")

cat("V_trait matrix saved to:", output_file, "\n")


