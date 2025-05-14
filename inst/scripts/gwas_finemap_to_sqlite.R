#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Read command line arguments
Trait1 <- as.character(args[1])  # GWAS trait name
block_file <- as.character(args[2])  # LD block file
sumstat_file <- as.character(args[3])  # Summary statistics file
output_dir <- as.character(args[4])  # Output directory
pip_threshold <- as.numeric(args[5])  # Minimum PIP to keep in final output
force_flag <- ifelse(length(args) >= 6, as.logical(args[6]), FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(RSQLite)
  library(Rmpfr)
  library(stringr)
})

# Create output directory if it doesn't exist
dir.create(sprintf("%s/SQLite", output_dir), showWarnings = FALSE, recursive = TRUE)

# # Directory containing the sumstats files
# dir_sumstats <- "/n/groups/price/ldsc/sumstats_formatted_2021/new"

# # List all files with pattern
# files <- list.files(dir_sumstats, pattern = "\\.sumstats_ext\\.gz$", full.names = TRUE)
# files <- files[-which(files == "/n/groups/price/ldsc/sumstats_formatted_2021/new/PASS_Autism_Grove2019.sumstats_ext.gz")]

# # Initialize a list to collect column names
# all_headers <- list()

# # Efficient header collection using readLines
# for (f in files) {
#   con <- gzfile(f, open = "rt")  # open connection in text mode
#   header_line <- tryCatch(readLines(con, n = 1), error = function(e) NULL)
#   close(con)  # close it explicitly
  
#   if (!is.null(header_line)) {
#     header <- strsplit(header_line, "\t|,")[[1]]
#     all_headers[[f]] <- header
#   }
# }

# # Flatten into a vector of column names
# colnames_all <- unlist(all_headers)
# unique_cols <- unique(colnames_all)

# # Helper function to search candidates
# find_matches <- function(patterns, choices) {
#   matches <- unlist(lapply(patterns, function(p) grep(p, choices, value = TRUE, ignore.case = TRUE)))
#   unique(matches)
# }

# # Define pattern lists (exclude indicates patterns that should not accompany the include pattern)
# patterns <- list(
#   CHR  = list(include = c("^chr", "^chrom", "^#chr", "^#CHR", "^#Chr", "^Chr", "^CHR"),
#               exclude = c("chromosome_length", "chromEnd")),
#   BP   = list(include = c("^bp$", "^pos", "^Pos", "^POS", "^base", "^BP"),
#               exclude = c("position_sd")),
#   A1   = list(include = c("^a1$", "^effect_allele$", "^A1$"),
#               exclude = c("a1_freq")),
#   A2   = list(include = c("^a2$", "^non_effect_allele", "^A2$"),
#               exclude = c("a2_freq")),
#   RSID = list(include = c("rsid", "^snp$", "markername", "variant_id"),
#               exclude = c("imputed")),
#   BETA = list(include = c("^beta", "^BETA", "^Beta", "^effect", "log_odds", "^Effect", "^or", "^OR", "^EFFECT", "^Z"),
#               exclude = c("frequency", "CI", "healthspan", "longevity", "lifespan", "EFFECTIVE")),
#   P    = list(include = c("^p$", "^pval", "^p_value", "^P", "^p"),
#               exclude = c("POS", "pos", "Pos"))
# )

# match_with_exclusion <- function(cols, include, exclude = NULL) {
#   any(
#     sapply(include, function(inc_pat) {
#       # Which columns match the include pattern?
#       inc_matches <- grep(inc_pat, cols, value = TRUE, perl = TRUE)
      
#       # Of those, filter out any that match an exclude pattern
#       if (!is.null(exclude)) {
#         inc_matches <- inc_matches[
#           !sapply(inc_matches, function(col)
#             any(sapply(exclude, function(exc_pat)
#               grepl(exc_pat, col, perl = TRUE)))
#           )
#         ]
#       }
      
#       length(inc_matches) > 0
#     })
#   )
# }

# # List of files missing some categories
# missing_by_file <- list()

# patterns2 <- patterns[-which(names(patterns) == "RSID")]
# for (f in names(all_headers)) {
#   cols <- all_headers[[f]]
#   missing_categories <- names(patterns2)[
#     sapply(patterns2, function(p) !match_with_exclusion(cols, p$include, p$exclude))
#   ]
#   if (length(missing_categories) > 0) {
#     missing_by_file[[f]] <- missing_categories
#   }
# }

# # View files with issues
# missing_by_file

# Prepare output file
db_path <- sprintf("%s/SQLite/%s_SQLite.db", output_dir, Trait1)

if (file.exists(db_path) && !force_flag) {
  cat("Output already exists at:", db_path, "\n")
  cat("Skipping computation. Use --force TRUE to overwrite.\n")
  quit(save = "no", status = 0)
} else if (file.exists(db_path) && force_flag) {
  cat("Overwriting existing file at:", db_path, "\n")
}

con <- dbConnect(SQLite(), db_path)

# Load summary statistics
sumstats <- fread(sumstat_file)

# Step 1: Compound columns
compound_patterns <- list(
  chr_pos_ = list(pattern = "chr_pos_", chr_fun = function(x) gsub("chr", "", x), bp_fun = identity),
  MarkerName = list(pattern = "MarkerName", chr_fun = identity, bp_fun = identity)
)

# Apply only the first matching transformation
for (entry in compound_patterns) {
  matches <- grep(entry$pattern, colnames(sumstats), value = TRUE)
  if (length(matches) > 0) {
    colname <- matches[1]
    split_vals <- tstrsplit(sumstats[[colname]], split = ":", fixed = TRUE)
    
    sumstats[, CHR := entry$chr_fun(split_vals[[1]])]
    sumstats[, BP := entry$bp_fun(split_vals[[2]])]
    sumstats[, (colname) := NULL]
    
    break  # Stop after first match
  }
}

# Step 2: SNP in chr:pos_ref_alt
if ("SNP" %in% colnames(sumstats)) {
  # Check for pattern like 1:12345_A_T or chr1:12345_A_T
  if (all(grepl("^[^:]+:[0-9]+_[ACGT]+_[ACGT]+$", sumstats$SNP))) {
    
    # First split on colon (chr:pos_ref_alt → c("chr", "pos_ref_alt"))
    split_chr_pos <- tstrsplit(sumstats$SNP, ":", fixed = TRUE)
    
    # Then split second part on underscore to get pos, ref, alt
    split_rest <- tstrsplit(split_chr_pos[[2]], "_", fixed = TRUE)
    
    sumstats[, CHR := gsub("^chr", "", split_chr_pos[[1]])]
    sumstats[, BP  := as.integer(split_rest[[1]])]
  }
}

# Step 3: Harmonize column names
patterns <- list(
  CHR  = list(include = c("^chr", "^chrom", "^#chr", "^#CHR", "^#Chr", "^Chr", "^CHR"),
              exclude = c("chromosome_length", "chromEnd")),
  BP   = list(include = c("^bp$", "^pos", "^Pos", "^POS", "^base", "^BP"),
              exclude = c("position_sd")),
  A1   = list(include = c("^a1$", "^effect_allele$", "^A1$"),
              exclude = c("a1_freq")),
  A2   = list(include = c("^a2$", "^non_effect_allele", "^A2$"),
              exclude = c("a2_freq")),
  RSID = list(include = c("rsid", "^snp$", "markername", "variant_id"),
              exclude = c("imputed")),
  BETA = list(include = c("^beta", "^BETA", "^Beta", "^effect", "log_odds", "^Effect", "^or", "^OR", "^EFFECT", "^Z"),
              exclude = c("frequency", "CI", "healthspan", "longevity", "lifespan", "EFFECTIVE")),
  P    = list(include = c("^p$", "^pval", "^p_value", "^P", "^p"),
              exclude = c("POS", "pos", "Pos"))
)

match_with_exclusion <- function(cols, include, exclude = NULL) {
  for (p in include) {
    matches <- grep(p, cols, value = TRUE, perl = TRUE)
    if (length(matches) > 0) {
      valid <- matches[!sapply(matches, function(m) any(sapply(exclude, function(e) grepl(e, m, perl = TRUE))))]
      if (length(valid) > 0) return(valid[1])
    }
  }
  return(NULL)
}

# Rename columns
for (target in names(patterns)) {
  matched_col <- match_with_exclusion(colnames(sumstats), patterns[[target]]$include, patterns[[target]]$exclude)
  if (!is.null(matched_col) && !(target %in% colnames(sumstats))) {
    setnames(sumstats, matched_col, target)
  }
}

# Step 4: Fallback for RSID via content check
if (!"RSID" %in% colnames(sumstats)) {
  for (col in colnames(sumstats)) {
    values <- sumstats[[col]]
    if (is.character(values) || is.factor(values)) {
      rs_fraction <- mean(grepl("^rs[0-9]+$", values))
      if (rs_fraction > 0.5) {
        setnames(sumstats, col, "RSID")
        break
      }
    }
  }
}

# Create RSID as CHR:BP:A1:A2 if not already present
if (!"RSID" %in% colnames(sumstats)) {
  if (all(c("CHR", "BP", "A1", "A2") %in% colnames(sumstats))) {
    sumstats[, RSID := paste0(CHR, ":", BP, ":", A1, ":", A2)]
  } else {
    warning("Cannot create RSID: missing one of CHR, BP, A1, or A2")
  }
}

# Final sumstat format
sumstats <- sumstats[, .(CHR, BP, RSID, A1, A2, BETA, P)]

# Function to perform single-causal variant fine-mapping from p-value
pval_to_posterior <- function(pval) {
  # Convert p-value to z-score
  # zsc <- qnorm(pval / 2, lower.tail = FALSE)
  zsc <- qnorm(as.numeric(log(pval) - log(2)), lower.tail = FALSE, log.p = TRUE)

  # Compute log probability from the normal distribution's PDF
  logprob <- -dnorm(zsc, log = TRUE)

  # Compute posterior probabilities using log-sum-exp normalization
  logsumexp <- function(log_values) {
    max_log <- max(log_values)
    max_log + log(sum(exp(log_values - max_log)))
  }

  posterior <- exp(logprob - logsumexp(logprob))
  return(posterior)
}

# Load LD blocks
blocks <- fread(block_file)[V2 %in% unique(sumstats$CHR), ]
setnames(blocks, c("LOCUS", "CHR", "START", "STOP"))
setkey(blocks, CHR, START, STOP)

# Fine-mapping across LD blocks
out_all <- data.table()
for (loc in 1:nrow(blocks)) {
  chr <- blocks[loc, CHR]
  start <- blocks[loc, START]
  stop <- blocks[loc, STOP]

  # Focus on SNPs within this block
  block_snps <- sumstats[CHR == chr & BP >= start & BP <= stop, ]
  if (nrow(block_snps) > 1) {
    # Convert p-values to high precision for stability
    pval_mpfr <- mpfr(block_snps$P, precBits = 100)
    pips <- pval_to_posterior(pval_mpfr)
    block_snps[, PIP := pips]
    block_snps <- block_snps[PIP > pip_threshold, .(CHR, BP, RSID, A1, A2, BETA, PIP)]
    out_all <- rbind(out_all, block_snps)
  }
}

# Save results to SQLite
dbWriteTable(con, "trait_data", out_all, overwrite = TRUE)
dbDisconnect(con)

cat("Data saved to:", db_path, "\n")

