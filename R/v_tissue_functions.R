#' @importFrom data.table := setkey copy foverlaps as.data.table

#' @title Load focal summary statistics
#'
#' @description
#' Loads summary statistics from an SQLite database for a given trait.
#'
#' @param trait Name of the focal trait
#' @param output_dir Directory containing SQLite databases
#'
#' @return A data.table of summary statistics
#' @export
load_summary_stats <- function(trait, output_dir) {
  db_path <- file.path(output_dir, "SQLite", sprintf("%s_SQLite.db", trait))
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  sumstats <- data.table::as.data.table(DBI::dbGetQuery(con, "SELECT CHR, BP, RSID, A1, A2, PIP, BETA FROM trait_data"))
  DBI::dbDisconnect(con)
  data.table::as.data.table(unique(sumstats))
}

#' Exclude SNPs overlapping exonic regions
#'
#' @param sumstats A data.table of summary statistics
#' @param exons A data.table of exon coordinates
#'
#' @return Filtered summary statistics data.table
#' @export
exclude_exonic_regions <- function(sumstats, exons) {
  sumstats_copy <- data.table::as.data.table(copy(sumstats))
  exons_copy <- data.table::as.data.table(copy(exons))
  sumstats_copy[, `:=`(CHR = as.numeric(CHR), START = BP, STOP = BP)]
  data.table::setkey(sumstats_copy, CHR, START, STOP)
  data.table::setkey(exons_copy, CHR, START, STOP)
  overlap <- data.table::foverlaps(sumstats_copy, exons_copy, type = "any", nomatch = 0L)
  overlap <- overlap[!is.na(START) & !BP == START]
  filtered <- sumstats_copy[!overlap[, .(CHR, i.START, i.STOP)], ]
  data.table::setkey(filtered, CHR, START, STOP)
  filtered
}

#' Load and combine epigenomic BED files, filter outliers
#'
#' @param dir_epi Directory containing BED files
#' @param pattern File pattern to match BED files
#' @param outlier_length Length threshold for excluding outlier peaks
#'
#' @return A cleaned data.table of combined tracks
#' @export
load_and_filter_tracks <- function(dir_epi, pattern = ".bed", outlier_length = 500000) {
  bed_files <- list.files(dir_epi, pattern = pattern, full.names = TRUE)
  bed_list <- lapply(bed_files, function(bed) {
    bedname <- gsub(pattern, "", basename(bed))
    bed_dt <- data.table::fread(bed, col.names = c("CHR", "START", "STOP"))
    bed_dt[, `:=`(CHR = as.numeric(gsub("chr", "", CHR)), START = as.numeric(START), STOP = as.numeric(STOP))]
    bed_dt <- bed_dt[CHR %in% 1:22]
    bed_dt[, track := bedname]
    bed_dt
  })
  combined <- data.table::rbindlist(bed_list)
  combined[, LENGTH := STOP - START]
  outlier_tracks <- combined[LENGTH > outlier_length, unique(track)]
  cleaned <- combined[!track %in% outlier_tracks]
  data.table::setkey(cleaned, CHR, START, STOP)
  cleaned
}

#' Load metadata file linking tracks to tissues
#'
#' @param dir_epi Directory containing metadata file
#' @param pattern Character pattern to match filename (default: "metadata")
#'
#' @return A data.table of metadata
#' @export
load_custom_metadata <- function(dir_epi, pattern = "metadata") {
  metadata_files <- list.files(dir_epi, pattern = pattern, full.names = TRUE)
  if (length(metadata_files) == 0) stop("No metadata file found matching pattern.")
  data.table::fread(metadata_files[1])
}

#' Load sample IDs associated with epigenomic tracks
#'
#' @param dir_epi Directory containing sample ID file
#'
#' @return A data.table of sample IDs
#' @export
load_custom_sample_ids <- function(dir_epi) {
  sample_path <- file.path(dir_epi, "Epimap_sampleID.txt")
  data.table::fread(sample_path, header = FALSE)
}

#' Compute SNP-track overlaps from BED data
#'
#' @param sumstats A data.table of summary statistics
#' @param tracks A data.table of epigenomic tracks
#'
#' @return A data.table with overlap annotations
#' @export
compute_track_overlaps <- function(sumstats, tracks) {
  sumstats_copy <- data.table::as.data.table(copy(sumstats))
  sumstats_copy[, `:=`(CHR = as.numeric(CHR), START = BP, STOP = BP)]
  data.table::setkey(sumstats_copy, CHR, START, STOP)
  data.table::setkey(tracks, CHR, START, STOP)
  overlaps <- data.table::foverlaps(sumstats_copy, tracks, type = "any", nomatch = NA)
  overlaps <- overlaps[!is.na(START) & !BP == START]
  overlap_pairs <- unique(overlaps[, .(RSID, track)])
  wide <- data.table::dcast(overlap_pairs, RSID ~ track, value.var = "track", fun.aggregate = length)
  #wide[, paste0("NA", "") := NULL]
  data.table::setkey(sumstats_copy, RSID)
  data.table::setkey(wide, RSID)
  merged <- wide[sumstats_copy]
  merged
}

#' Compute background expectation from 1KG chr1 and tracks
#'
#' @param dir_epi Directory containing 1KG .bim file
#' @param tracks Combined epigenomic tracks
#' @param exons A data.table of exon coordinates
#'
#' @return A named numeric vector of background expectations
#' @export
compute_background_matrix <- function(dir_epi, tracks, exons) {
  bckg_path <- file.path(dir_epi, "1000G.EUR.QC.1.bim")
  bckg <- data.table::fread(bckg_path)
  colnames(bckg) <- c("CHR", "RSID", "CM", "BP", "REF", "ALT")
  bckg[, `:=`(CHR = as.numeric(CHR), START = BP, STOP = BP)]
  overlaps <- compute_track_overlaps(bckg, tracks)
  data.table::setkey(overlaps, CHR, START, STOP)
  data.table::setkey(exons, CHR, START, STOP)
  no_exon <- data.table::foverlaps(overlaps, exons, type = "any", nomatch = 0L)
  bckg_noexon <- overlaps[!no_exon[, .(CHR, i.START, i.STOP)], ]
  compute_bckg(bckg_noexon)
}

#' Build track-to-tissue metadata from sample IDs
#'
#' @param sumstats Summary stats with overlap columns
#' @param metadata Metadata table linking tracks to tissues
#' @param sampleID Sample ID table
#'
#' @return A data.table linking tracks to tissue groupings
#' @export
build_track_metadata <- function(sumstats, metadata, sampleID) {
  ID <- data.table::as.data.table(copy(jpepR::sampleID))
  meta <- data.table::as.data.table(copy(jpepR::metadata))
  ID[, id := 1:.N]
  data.table::setnames(ID, "V1", "SampleID")
  tracks <- setdiff(colnames(sumstats), c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA", "START", "STOP"))
  track_id <- as.numeric(sapply(strsplit(tracks, "[.]"), `[[`, 3))
  track_info <- data.table::data.table(track = tracks, id = track_id)
  track_info <- merge(track_info, ID, by = "id", all.x = TRUE)
  meta <- meta[id %in% track_info$SampleID, ]
  meta[, SampleID := as.character(id)]
  meta[is.na(SampleID), SampleID := Annot]
  track_info <- merge(track_info, meta[, .(SampleID, GROUP, infoline)], by = "SampleID", all.x = TRUE)
  track_info[GROUP != "Cancer"]
}

#' Prepare EM input matrix from filtered overlaps
#'
#' @param sumstats_epi Annotated SNP-track matrix
#' @param track_info Track metadata table
#'
#' @return A binary matrix ready for EM
#' @export
prepare_em_input <- function(sumstats_epi, track_info) {
  sum <- data.table::as.data.table(copy(sumstats_epi))
  valid_tracks <- intersect(track_info$track, colnames(sum))
  mat <- as.matrix(sum[, ..valid_tracks])
  mat[is.na(mat)] <- 0
  mat
}

#' Normalize conditional probability table (CPT)
#'
#' @param mat Matrix of SNP-track binary overlaps
#'
#' @return A normalized CPT matrix (row-normalized)
#' @export
normalize_cpt <- function(mat) {
  cpt <- mat / rowSums(mat)
  cpt[is.na(cpt)] <- 0
  cpt
}

#' Run EM algorithm to estimate SNP-to-track posteriors
#'
#' @param cpt Normalized conditional probability table
#' @param pip Vector of PIPs for each SNP
#' @param theta0 Initial prior values
#' @param bckg_expectation Vector of background expectations
#' @param tol Convergence tolerance
#' @param maxiter Maximum number of iterations
#'
#' @return A list with final posterior matrix and iteration count
#' @export
run_em_algorithm <- function(cpt, pip, theta0, bckg_expectation, tol = 1e-8, maxiter = 200) {
  apt <- cpt
  theta <- theta0
  diff <- 2
  iter <- 0
  while (diff > tol && iter < maxiter) {
    apt_prev <- apt
    denominator <- 0
    for (t in seq_len(ncol(apt))) {
      weight_pip <- pip
      weight_size <- bckg_expectation[[colnames(apt)[t]]]
      numerator <- sum(apt[, t] * weight_pip / weight_size)
      denominator <- denominator + numerator
      theta[t] <- numerator
    }
    theta <- theta / denominator
    for (t in seq_len(ncol(apt))) {
      apt[, t] <- cpt[, t] * theta[t]
    }
    apt <- apt / rowSums(apt)
    apt[is.na(apt)] <- 0
    diff <- mean((apt - apt_prev)^2)
    iter <- iter + 1
  }
  list(apt = apt, iter = iter)
}

#' Aggregate SNP-level posteriors into tissue groups
#'
#' @param em_matrix Matrix of SNP-by-track posteriors
#' @param snp_info Summary stats with SNP identifiers
#' @param track_info Table linking tracks to tissues
#'
#' @return A data.table with SNP-by-tissue enrichment scores
#' @export
aggregate_by_tissue <- function(em_matrix, snp_info, track_info) {
  tissues <- unique(track_info$GROUP)
  snp <- copy(snp_info)
  v_tissue <- snp[, .(CHR, BP, RSID, A1, A2, PIP, BETA, START, STOP)]
  for (t in tissues) {
    tissue_tracks <- track_info[GROUP == t, track]
    indices <- which(colnames(em_matrix) %in% tissue_tracks)
    if (length(indices) > 1) {
      snp_count <- rowSums(em_matrix[, indices, drop = FALSE])
    } else {
      snp_count <- em_matrix[, indices]
    }
    v_tissue[[t]] <- snp_count
  }
  v_tissue[, `:=`(START = NULL, STOP = NULL)]
  v_tissue
}
