#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# ─────────────────────────────────────────────────────
# Parse Arguments
focal       <- args[1]                      # Trait name
epi_arg      <- args[2]                      # "default" or path to custom epi BEDs
pip_threshold <- as.numeric(args[3])         # Minimum PIP threshold
output_dir   <- args[4]                      # Output directory

# ─────────────────────────────────────────────────────
# Setup
suppressPackageStartupMessages({
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(jpepR)
})

# Directories
dir_epi_default <- file.path(output_dir, "../epimap")
dir_epi         <- if (epi_arg == "default") dir_epi_default else epi_arg
use_default     <- epi_arg == "default"

# Output location
dir_output <- file.path(output_dir, "V_tissues")
dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────
# Step 1: Load Summary Statistics
sumstats <- load_summary_stats(focal, output_dir)

# ─────────────────────────────────────────────────────
# Step 2: Exclude Exonic Regions
exon_dat <- data.table::as.data.table(jpepR::exons)
colnames(exon_dat) <- c("CHR", "START", "STOP", "gene", "b", "strand")
suppressWarnings(exon_dat[, CHR := as.numeric(gsub("chr", "", CHR))])
exon_dat <- exon_dat[CHR %in% 1:22]
sumstats <- exclude_exonic_regions(sumstats, exons = exon_dat)

# ─────────────────────────────────────────────────────
# Step 3: Load Tracks (default or custom)
if (use_default) {
  message("Using default EpiMap tracks from jpepR package")
  message("\nTracks have been filtered to retain only those overlapping the example SNPs")
  data("tracks", package = "jpepR")
  combined_bed <- tracks
  data("bckg_expectation", package = "jpepR")
  metadata     <- get("metadata", envir = asNamespace("jpepR"))
  sampleID     <- get("sampleID", envir = asNamespace("jpepR"))
} else {
  message("Using custom epigenomic tracks from ", dir_epi)
  combined_bed <- load_and_filter_tracks(dir_epi)
  metadata     <- load_custom_metadata(dir_epi)
  sampleID     <- load_custom_sample_ids(dir_epi)
}

# ─────────────────────────────────────────────────────
# Step 4: Compute Overlap with Summary Stats
sumstats_epi <- compute_track_overlaps(sumstats, combined_bed)

# ─────────────────────────────────────────────────────
# Step 5: Compute Background (if using custom epi)
if (!use_default) {
  bckg_expectation <- compute_background_matrix(dir_epi, combined_bed, exons)
}

# ─────────────────────────────────────────────────────
# Step 6: Prepare Metadata
track_info <- build_track_metadata(sumstats_epi, metadata, sampleID)

# ─────────────────────────────────────────────────────
# Step 7: Run EM Algorithm
mat_EM <- prepare_em_input(sumstats_epi, track_info)
pip    <- sumstats_epi$PIP
theta0 <- rep(1 / ncol(mat_EM), ncol(mat_EM))
cpt    <- normalize_cpt(mat_EM)
em_out <- run_em_algorithm(cpt, pip, theta0, bckg_expectation)

# ─────────────────────────────────────────────────────
# Step 8: Aggregate by Tissue
v_tissue <- aggregate_by_tissue(
  em_matrix = em_out$apt,
  snp_info = sumstats_epi,
  track_info = track_info
)

# ─────────────────────────────────────────────────────
# Step 9: Write Output
output_file <- file.path(dir_output, sprintf("%s_V_tissues_thres%s.tsv.gz", focal, pip_threshold))
fwrite(v_tissue, output_file, sep = "\t")

message("✅ V_tissue matrix saved to: ", output_file)
