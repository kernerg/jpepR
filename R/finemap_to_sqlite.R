#' Fine-map summary statistics and save to SQLite
#'
#' @param trait Trait name
#' @param sumstat_file Path to .sumstats.gz file
#' @param ld_blocks Path to LD block file, or NULL to use internal default
#' @param output_dir Output directory
#' @param pip_threshold Minimum PIP threshold (default = 0.01)
#' @param force Logical. Whether to overwrite existing SQLite file (default = FALSE)
#'
#' @import data.table
#' @import DBI
#' @import RSQLite
#' @importFrom Rmpfr mpfr
#' @export
finemap_to_sqlite <- function(trait, sumstat_file, ld_blocks, output_dir,
                              pip_threshold = 0.01, force = FALSE) {

  # Create output path
  db_path <- file.path(output_dir, "SQLite", sprintf("%s_SQLite.db", trait))
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)

  # ─────────────────────────────────────────────────────
  # Step 1: Skip or overwrite existing output
  if (file.exists(db_path) && !force) {
    message("SQLite already exists: ", db_path, ". Use `force = TRUE` to overwrite.")
    return(invisible(NULL))
  }

  # ─────────────────────────────────────────────────────
  # Step 2: Load summary statistics and handle column names
  sumstats <- data.table::fread(sumstat_file)

  # Handle compound columns
  compound_patterns <- list(
    chr_pos_ = list(pattern = "chr_pos_", chr_fun = function(x) gsub("chr", "", x), bp_fun = identity),
    MarkerName = list(pattern = "MarkerName", chr_fun = identity, bp_fun = identity)
  )

  for (entry in compound_patterns) {
    matches <- grep(entry$pattern, colnames(sumstats), value = TRUE)
    if (length(matches) > 0) {
      colname <- matches[1]
      split_vals <- tstrsplit(sumstats[[colname]], ":", fixed = TRUE)
      sumstats[, CHR := entry$chr_fun(split_vals[[1]])]
      sumstats[, BP := entry$bp_fun(split_vals[[2]])]
      sumstats[, (colname) := NULL]
      break
    }
  }

  # Parse chr:pos_ref_alt format
  if ("SNP" %in% colnames(sumstats)) {
    if (all(grepl("^[^:]+:[0-9]+_[ACGT]+_[ACGT]+$", sumstats$SNP))) {
      split_chr_pos <- tstrsplit(sumstats$SNP, ":", fixed = TRUE)
      split_rest <- tstrsplit(split_chr_pos[[2]], "_", fixed = TRUE)
      sumstats[, CHR := gsub("^chr", "", split_chr_pos[[1]])]
      sumstats[, BP  := as.integer(split_rest[[1]])]
    }
  }

  # Harmonize column names
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

  for (target in names(patterns)) {
    matched_col <- match_with_exclusion(colnames(sumstats), patterns[[target]]$include, patterns[[target]]$exclude)
    if (!is.null(matched_col) && !(target %in% colnames(sumstats))) {
      setnames(sumstats, matched_col, target)
    }
  }

  # Fallback for RSID
  if (!"RSID" %in% colnames(sumstats)) {
    for (col in colnames(sumstats)) {
      values <- sumstats[[col]]
      if (is.character(values) || is.factor(values)) {
        if (mean(grepl("^rs[0-9]+$", values)) > 0.5) {
          setnames(sumstats, col, "RSID")
          break
        }
      }
    }
  }

  # Create RSID if missing
  if (!"RSID" %in% colnames(sumstats)) {
    if (all(c("CHR", "BP", "A1", "A2") %in% colnames(sumstats))) {
      sumstats[, RSID := paste0(CHR, ":", BP, ":", A1, ":", A2)]
    } else {
      warning("Cannot create RSID: missing one of CHR, BP, A1, or A2")
    }
  }

  # Validate
  jpepR::validate_sumstats(sumstats)
  
  # Final formatting
  sumstats <- sumstats[, .(CHR, BP, RSID, A1, A2, BETA, P)]

  # ─────────────────────────────────────────────────────
  # Step 3: LD block fine-mapping
  if (is.null(ld_blocks)) {
    message("Using built-in LD blocks")
    blocks <- jpepR::ld_blocks
    data.table::setnames(blocks, c("LOCUS", "CHR", "START", "STOP"))
  } else if (is.character(ld_blocks) && file.exists(ld_blocks)) {
    blocks <- data.table::fread(ld_blocks)
  } else {
    stop("`ld_blocks` must be NULL or a valid file path.")
  }
 
  # Validate
  jpepR::validate_ld_blocks(blocks)
  data.table::setkey(blocks, CHR, START, STOP)

  # PIP calculation
  pval_to_posterior <- function(pval) {
    zsc <- qnorm(as.numeric(log(pval) - log(2)), lower.tail = FALSE, log.p = TRUE)
    logprob <- -dnorm(zsc, log = TRUE)
    logsumexp <- function(logv) max(logv) + log(sum(exp(logv - max(logv))))
    exp(logprob - logsumexp(logprob))
  }

  out_all <- data.table()
  for (i in 1:nrow(blocks)) {
    chr <- blocks[i, CHR]
    start <- blocks[i, START]
    stop <- blocks[i, STOP]

    block_snps <- sumstats[CHR == chr & BP >= start & BP <= stop]
    if (nrow(block_snps) > 1) {
      pval_mpfr <- mpfr(block_snps$P, precBits = 100)
      pips <- pval_to_posterior(pval_mpfr)
      block_snps[, PIP := pips]
      block_snps <- block_snps[PIP > pip_threshold, .(CHR, BP, RSID, A1, A2, BETA, PIP)]
      out_all <- rbind(out_all, block_snps)
    }
  }

  # ─────────────────────────────────────────────────────
  # Step 4: Save results to SQLite
  con <- dbConnect(RSQLite::SQLite(), db_path)
  dbWriteTable(con, "trait_data", out_all, overwrite = TRUE)
  dbDisconnect(con)

  message("✅ Fine-mapped data saved to: ", db_path)
}
