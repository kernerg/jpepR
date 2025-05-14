#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Read command line arguments
Trait1 <- as.character(args[1])  # GWAS trait name
epi_arg <- args[2] # Directory for epi bed files
pip_threshold <- as.numeric(args[3])  # Minimum PIP to keep in final output
output_dir <- as.character(args[4])  # Output directory

dir_epi_default <- sprintf("%s/../epimap", output_dir)
dir_epi <- if (epi_arg == "default") dir_epi_default else as.character(epi_arg) 
use_default_epimap <- (epi_arg == "default")

# Paths to default precomputed files
default_epi_file <- sprintf("%s/Epimap_tracks_overlapping_focal.tsv.gz", dir_epi_default)
default_bckg_epi_file <- sprintf("%s/Epimap_bckg_expectation.rds", dir_epi_default)

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(glue)
})

# Set of functions
# Function to load epigenomic tracks, with outlier removal option
load_epi_tracks <- function(bed_files, pattern, outlier_length = 500000) {
  bed_list <- lapply(bed_files, function(bed) {
    bedname <- gsub(pattern, "", basename(bed))
    bed_dt <- fread(bed, col.names = c("CHR", "START", "STOP"))
    bed_dt[, `:=`(CHR = gsub("chr", "", CHR), START = as.numeric(START), STOP = as.numeric(STOP))]
    bed_dt <- bed_dt[CHR %in% 1:22]
    bed_dt[, track := bedname]
    return(bed_dt)
  })
    
  # Combine all BED files into one large data.table
  combined_bed <- rbindlist(bed_list)
  setkey(combined_bed, CHR, START, STOP)

  # Remove outliers (default: peak length >500 kb)
  combined_bed[, LENGTH := STOP-START]
  outlier_tracks <- combined_bed[LENGTH > outlier_length, unique(track)]
  combined_bed <- combined_bed[!track %in% outlier_tracks]

  return(list(combined = combined_bed, outlier = outlier_tracks))
}

overlap_assignment <- function(query_dt, overlaps) {
  # Extract unique (RSID, track) pairs
  overlap_pairs <- unique(overlaps[, .(RSID, track)])
  
  # Reshape to a wide format to get one row per SNP
  wide_tracks <- dcast(overlap_pairs, RSID ~ track, value.var = "track", fun.aggregate = length)
  wide_tracks[, paste0("NA", "") := NULL]

  # Merge back with the original data
  setkey(query_dt, RSID)
  setkey(wide_tracks, RSID)
  query_dt <- wide_tracks[query_dt]
  
  return(query_dt)
}

# Function to merge summary statistics with epigenomic tracks
compute_overlaps_bulk <- function(query_dt, combined_bed) {
  query_dt[, `:=`(CHR = as.numeric(CHR), START = as.numeric(BP), STOP = as.numeric(BP))]
  setkey(query_dt, CHR, START, STOP)
  
  # Perform foverlaps call
  overlaps <- foverlaps(query_dt, combined_bed, type = "any", nomatch = NA)

  # Annotate overlaping tracks with 1
  query_dt_overlap <-  overlap_assignment(query_dt, overlaps)
  
  return(query_dt_overlap)
}

# Function to compute background expectation for epigenomic tracks (using chr1 of 1KG)
compute_bckg <- function(mat) {

  # track names
  non_tracks <- c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA", "START", "STOP", "CM", "REF", "ALT")
  tracks <- setdiff(colnames(mat), non_tracks)

  # Define chunk size (adjust based on available memory)
  chunk_size <- 10000  # Number of rows per chunk
  n_rows <- nrow(mat)

  # Initialize storage for column sums
  track_sums <- rep(0, length(tracks))
  total_denom <- 0

  # Process data in chunks
  for (i in seq(1, n_rows, by = chunk_size)) {
    # Define the row range for this chunk
    row_range <- i:min(i + chunk_size - 1, n_rows)
    
    # Extract the chunk and convert to matrix
    sub_chunk <- as.matrix(mat[row_range, ..tracks])
    
    # Compute row sums for normalization
    row_sums <- rowSums(sub_chunk)
    row_sums[row_sums == 0] <- 1  # Handle division by zero
    
    # Normalize the chunk
    sub_chunk <- sub_chunk / row_sums
    
    # Accumulate column sums
    track_sums <- track_sums + colSums(sub_chunk)
    
    # Accumulate total denominator
    total_denom <- total_denom + sum(colSums(sub_chunk))
  }

  # Final normalized expectation
  bckg_expectation <- track_sums / total_denom

  return(bckg_expectation)
}

EM_algo <- function(binary_mat, priors, bckg_expectation,
                    tol = 1e-8, maxiter = 200) {
  ## Algorithm starts
  apt <- binary_mat
  theta <- priors
  diff <- 2
  iter <- 0
  while (diff > tol && iter < maxiter) {
    
    ## save previous instance to compute distance between consec iterations
    apt_prev <- apt
    
    ## update prior values
    denominator <- 0
    for (t in seq_len(ncol(apt))) {
      weight_pip <- pip ## use pip to upweight SNP
      weight_size <- bckg_expectation[[colnames(apt)[t]]] # use size to downweight track
      numerator <- sum(apt[, t] * weight_pip / weight_size)
      denominator <- denominator + numerator
      theta[t] <- numerator
    }
    theta <- theta / denominator
    
    ## update posteriors
    for (t in seq_len(ncol(apt))) {
      apt[, t] <- cpt[, t] * theta[t]
    }
    apt <- as.matrix(apt / rowSums(apt))
    apt[is.nan(apt)] <- 0
    
    diff <- mean((apt - apt_prev)^2)
    iter <- iter + 1
  }
  return(list(apt = apt, iter = iter))
}

# Set directories
db_path <- sprintf("%s/SQLite/%s_SQLite.db", output_dir, Trait1)
dir_output <- sprintf("%s/V_tissues", output_dir)
dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

# Load summary statistics for the focal trait
con <- dbConnect(SQLite(), db_path)
sumstats <- as.data.table(dbGetQuery(con, "SELECT CHR, BP, RSID, A1, A2, PIP, BETA FROM trait_data"))
dbDisconnect(con)
sumstats <- unique(sumstats)

# exclude exonic regions
data(exons)
colnames(exons) <- c("CHR", "START", "STOP", "gene", "b", "strand")
exons[, CHR := as.numeric(gsub("chr", "", CHR))]
exons <- exons[CHR %in% 1:22]
setkey(exons, CHR, START, STOP)

# Load epigenomic annotations
if (use_default_epimap) {
  cat("Reading default epigenomic files from: ", default_epi_file, "\n")
  combined_bed <- fread(default_epi_file)
  setkey(combined_bed, CHR, START, STOP)
  bckg_expectation <- readRDS(default_bckg_epi_file)
} else {
  folder <- dir_epi
  pattern <- ".bed"
  #pattern <- "_liftedOverTo37.bed3"
  # List all BED files
  bed_files <- list.files(folder, pattern = pattern, full.names = TRUE)

  # Bulk load all BED files
  out_epi <- load_epi_tracks(bed_files, pattern)
  combined_bed <- out_epi$combined
  setkey(combined_bed, CHR, START, STOP)
  outlier_tracks <- out_epi$outlier

  cat("Read epigenomic files from: ", folder, "\n")

  # read bckg files
  cat("Reading 1KG chr1 to compute background expectations", "\n")

  bckg <- fread(sprintf("%s/1000G.EUR.QC.1.bim", dir_external))
  colnames(bckg) <- c("CHR", "RSID", "CM", "BP", "REF", "ALT")

  # Obtain raw overlaps between 1KG chr1 and epigenomic tracks for bckg expectation
  bckg_epi <- compute_overlaps_bulk(bckg, combined_bed)
  setkey(bckg_epi, CHR, START, STOP)
  no_exon_bckg <- foverlaps(bckg_epi, exons, type = "any", nomatch = 0L)
  bckg_epi_noexon <- bckg_epi[!no_exon_bckg[, .(CHR, i.START, i.STOP)], ]
  setkey(bckg_epi_noexon, CHR, START, STOP)
}

# Obtain raw overlaps between summary statistics and epigenomic tracks
cat("Merging epigenomic files with sumstats", "\n")
sumstats_epi <- compute_overlaps_bulk(sumstats, combined_bed)

# Remove duplicate SNPs, keep highest PIP
setorder(sumstats_epi, CHR, START, STOP, -PIP)

# Filter by PIP threshold
sumstats_epi <- sumstats_epi[PIP > pip_threshold]

# Set keys for efficient overlapping using data.table
setkey(sumstats_epi, CHR, START, STOP)

# Perform the anti-join to find regions in sumstats that do NOT overlap with exons
no_exon <- foverlaps(sumstats_epi, exons, type = "any", nomatch = 0L)

# Subset sumstats to keep only regions that have no match with exons
sumstats_epi_noexon <- sumstats_epi[!no_exon[, .(CHR, i.START, i.STOP)], ]
setkey(sumstats_epi_noexon, CHR, START, STOP)

# Load biosample metadata
data(metadata)
data(sampleID)
mark_matrix_names <- sampleID
mark_matrix_names[, id := 1:.N]
setnames(mark_matrix_names, "V1", "SampleID")

# merge tracks with metadata
non_tracks <- c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA", "START", "STOP")
tracks <- setdiff(colnames(sumstats_epi_noexon), non_tracks)
track_id <- as.numeric(sapply(strsplit(tracks, split = "[.]"), `[[`, 3))
track_info <- data.table(track = tracks, id = track_id)
track_info <- merge(track_info, mark_matrix_names, by = "id", all.x = TRUE)

# Load track metadata
metadata <- metadata[id %in% track_info$SampleID, ]
metadata[, SampleID := as.character(id)]
metadata[is.na(SampleID), SampleID := Annot]

# Remove cancer samples
track_info <- merge(track_info, metadata[, .(SampleID, GROUP, infoline)], by = "SampleID", all.x = TRUE)
cancer_samples <- track_info[GROUP == "Cancer", track]
tracks <- setdiff(track_info[, track], cancer_samples)
track_info <- track_info[track %in% tracks]
sumstats_epi_noexon[, paste0(cancer_samples, "") := NULL]

# Measured tissues
tissues <- track_info[, unique(GROUP)]

# Compute background expectations (input to EM algorithm)
if (!use_default_epimap) {
  bckg_epi_noexon[, paste0(cancer_samples, "") := NULL]
  bckg_expectation <- compute_bckg(bckg_epi_noexon)
  cat("Background expectations computed", "\n")
}

# Subset to valid annotations
valid_tracks <- intersect(tracks, colnames(sumstats_epi_noexon))
mat_EM <- sumstats_epi_noexon[, ..valid_tracks]
pip <- sumstats_epi_noexon[, PIP]

# Initialize priors (uninformative by default) (input to EM algorithm)
theta0 <- rep(1 / ncol(mat_EM), ncol(mat_EM))

# Normalize the SNP-to-tissue matrix (convert to conditional probability table)
cpt <- mat_EM / rowSums(mat_EM)
cpt <- as.matrix(cpt)
cpt[is.nan(cpt)] <- 0  # Replace NaNs with 0

cat("SNP-to-tissue matrix (cpt) prepared with", ncol(cpt), "annotations and", nrow(cpt), "SNPs.\n")

# Run EM algorithm to obtain SNP-to-tracks posteriors
out_em <- EM_algo(cpt, priors = theta0, bckg_expectation = bckg_expectation)
apt <- out_em$apt
iter_nbr <- out_em$iter

# Group results by tissue category
info_columns <- intersect(non_tracks, colnames(sumstats_epi_noexon))
v_tissue <- sumstats_epi_noexon[, ..info_columns]

# sum posteriors across tissue groups
for (t in tissues) {
  tissue_tracks <- track_info[GROUP == t, track]
  if(length(which(colnames(apt) %in% tissue_tracks)) > 1){
    snp_count <- rowSums(apt[,which(colnames(apt) %in% tissue_tracks)])
  } else{
    snp_count <- apt[,which(colnames(apt) %in% tissue_tracks)]
  }
  v_tissue[, paste0(t, "") := snp_count]
}
v_tissue[,`:=`(START = NULL, STOP = NULL)]

# Write to gzipped file
output_file <- sprintf("%s/%s_V_tissues_thres%s.tsv.gz", dir_output, Trait1, pip_threshold)
fwrite(v_tissue, file = output_file, sep = "\t")

cat("V_tissue matrix saved to:", output_file, "\n")

# # Read the file as a data.table
# V_tissue <- fread(output_file, sep = "\t", header = TRUE)
