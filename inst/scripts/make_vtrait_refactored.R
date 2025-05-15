#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# ─────────────────────────────────────────────────────
# Parse Arguments
focal          <- args[1]                      # Focal trait name
aux_traits     <- strsplit(args[2], ",")[[1]]  # Vector of auxiliary trait names
pip_threshold  <- as.numeric(args[3])          # PIP threshold
aux_pip_thres  <- as.numeric(args[4])          # PIP cutoff for auxiliary traits
aux_nbr_thres  <- as.numeric(args[5])          # Minimum number of non-zero SNPs
output_dir     <- args[6]                      # Output directory

# ─────────────────────────────────────────────────────
# Setup
suppressPackageStartupMessages({
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(jpepR)
})

# Directories
dir_sqlite <- file.path(output_dir, "SQLite")
dir_output <- file.path(output_dir, "V_traits")
dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────
# Step 1: Load focal trait summary stats
focal_data <- load_summary_stats(focal, output_dir)
focal_data <- focal_data[PIP > pip_threshold]
focal_data <- unique(focal_data, by = "RSID")
setorder(focal_data, RSID)

# ─────────────────────────────────────────────────────
# Step 2: Initialize V_traits matrix
V_traits <- focal_data[, .(CHR, BP, RSID, A1, A2, PIP, BETA)]

# ─────────────────────────────────────────────────────
# Step 3: Process auxiliary traits
aux_traits <- setdiff(gsub("_SQLite.db", "", aux_traits), focal)
for (aux in aux_traits) {
  message("Processing auxiliary trait: ", aux)
  aux_data <- load_summary_stats(aux, output_dir)[, .(RSID, A1, A2, PIP, BETA)]

  merged <- merge(V_traits, aux_data, by = "RSID", all.x = TRUE, suffixes = c("", "_aux"))
  flip <- !is.na(merged$A1_aux) & (merged$A1 != merged$A1_aux)
  merged[flip, `:=`(BETA_aux = -BETA_aux, A1_aux = A1, A2_aux = A2)]

  pos_col <- paste0(aux, "_pos")
  neg_col <- paste0(aux, "_neg")

  merged[, (pos_col) := pmax(PIP * PIP_aux * sign(BETA * BETA_aux), 0)]
  merged[, (neg_col) := pmax(-PIP * PIP_aux * sign(BETA * BETA_aux), 0)]
  merged[is.na(get(pos_col)), (pos_col) := 0]
  merged[is.na(get(neg_col)), (neg_col) := 0]

  new_cols <- merged[, .SD, .SDcols = c(pos_col, neg_col)]
  V_traits <- cbind(V_traits, as.matrix(new_cols))
}

# ─────────────────────────────────────────────────────
# Step 4: Filter by signal strength
aux_cols <- setdiff(colnames(V_traits), c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA"))
keep_rows <- rowSums(V_traits[, ..aux_cols]) > 0
V_traits <- V_traits[keep_rows]

filt_aux <- V_traits[, lapply(.SD, function(x) sum(x > aux_pip_thres)), .SDcols = aux_cols]
keep_cols <- names(filt_aux)[filt_aux > aux_nbr_thres]
keep_all <- c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA", keep_cols)
V_traits <- V_traits[, ..keep_all]

# ─────────────────────────────────────────────────────
# Step 5: Save output
output_file <- file.path(dir_output, sprintf("%s_V_traits_thres%s.tsv.gz", focal, pip_threshold))
fwrite(V_traits, file = output_file, sep = "\t")
message("✅ V_trait matrix saved to: ", output_file)
