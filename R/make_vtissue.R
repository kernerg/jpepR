#' Construct V_tissue matrix from summary statistics and epigenomic tracks
#'
#' @param focal_trait Name of the focal trait
#' @param output_dir Output directory
#' @param epi_dir Path to custom epigenomic BED folder or "default" for built-in EpiMap
#' @param pip_threshold Minimum PIP threshold to retain SNPs (default = 0.01)
#'
#' @import data.table
#' @return Path to saved V_tissues file
#' @export
make_vtissue <- function(focal_trait, output_dir, epi_dir = "default", pip_threshold = 0.01) {
  use_default <- epi_dir == "default"
  dir_output <- file.path(output_dir, "V_tissues")
  dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

  # ─────────────────────────────────────────────────────
  # Step 1: Load and Filter Summary Statistics From Fine-Mapping file
  sumstats <- jpepR::load_summary_stats(focal_trait, output_dir)
  sumstats <- sumstats[PIP > pip_threshold]
  
  # ─────────────────────────────────────────────────────
  # Step 2: Exclude Exonic Regions
  exon_dat <- jpepR::exons
  colnames(exon_dat) <- c("CHR", "START", "STOP", "gene", "b", "strand")
  suppressWarnings(exon_dat[, CHR := as.numeric(gsub("chr", "", CHR))])
  exon_dat <- exon_dat[CHR %in% 1:22]
  sumstats <- jpepR::exclude_exonic_regions(sumstats, exons = exon_dat)

  # ─────────────────────────────────────────────────────
  # Step 3: Load Tracks (default or custom)
  if (use_default) {
    message("Using default EpiMap tracks from jpepR package")
    #combined_bed     <- jpepR::tracks                      # BED-style data
    combined_bed     <- fread("data-raw/Epimap_tracks.tsv.gz")
    bckg_expectation <- jpepR::bckg_expectation            # background expectations
    metadata         <- jpepR::metadata 
    sampleID         <- jpepR::sampleID
  } else {
    message("Using custom epigenomic tracks from: ", epi_dir)
    combined_bed     <- jpepR::load_and_filter_tracks(epi_dir)
    metadata         <- jpepR::load_custom_metadata(epi_dir)
    sampleID         <- jpepR::load_custom_sample_ids(epi_dir)
    bckg_expectation <- jpepR::compute_background_matrix(epi_dir, combined_bed, exon_dat)
  }

  # Validate
  jpepR::validate_epigenomic_tracks(combined_bed)

  # ─────────────────────────────────────────────────────
  # Step 4: Compute Overlap with Summary Stats
  sumstats_epi <- jpepR::compute_track_overlaps(sumstats, combined_bed)

  # ─────────────────────────────────────────────────────
  # Step 5: Prepare Metadata
  track_info <- jpepR::build_track_metadata(sumstats_epi, metadata, sampleID)
  track_info <- track_info[track %in% names(bckg_expectation)]

  # ─────────────────────────────────────────────────────
  # Step 6: Run EM Algorithm
  mat_EM <- jpepR::prepare_em_input(sumstats_epi, track_info)
  pip    <- sumstats_epi$PIP
  theta0 <- rep(1 / ncol(mat_EM), ncol(mat_EM))
  cpt    <- jpepR::normalize_cpt(mat_EM)
  em_out <- jpepR::run_em_algorithm(cpt, pip, theta0, bckg_expectation)

  # ─────────────────────────────────────────────────────
  # Step 7: Aggregate by Tissue
  v_tissue <- jpepR::aggregate_by_tissue(
    em_matrix  = em_out$apt,
    snp_info   = sumstats_epi,
    track_info = track_info
  )

  # ─────────────────────────────────────────────────────
  # Step 9: Write Output
  output_file <- file.path(dir_output, sprintf("%s_V_tissues_thres%s.tsv.gz", focal_trait, pip_threshold))
  data.table::fwrite(v_tissue, output_file, sep = "\t")

  message("✅ V_tissue matrix saved to: ", output_file)
  return(output_file)
}
