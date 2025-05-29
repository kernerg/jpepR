#' Construct V_trait matrix from PIP effects
#'
#' For a given focal trait and a list of auxiliary traits, this function computes the
#' directional PIP product for all SNPs shared across the traits, and saves the result.
#' If custom fine-mapped files are provided, they are converted to SQLite first.
#'
#' @param focal_trait Name of the focal trait
#' @param aux_traits Character vector of auxiliary trait names
#' @param output_dir Path to output directory
#' @param pip_threshold Minimum PIP for the focal trait
#' @param aux_pip_thres Minimum PIP for auxiliary traits
#' @param aux_nbr_thres Minimum number of auxiliary traits with strong signal per SNP
#' @param fine_mapped_files Optional named list with paths to custom fine-mapped files. Names must match traits.
#' @param full_matrix_path Optional path to a precomputed full pleiotropy matrix (with columns RSID, *_pos, *_neg, A1)
#' @param subset Optional vector with a subset of traits from the precomputed full pleiotropy matrix (character vector with names)
#'
#' @import data.table
#' @import DBI
#' @import RSQLite
#' @return Path to the saved V_trait matrix
#' @export
make_vtrait <- function(focal_trait,
                        aux_traits,
                        output_dir,
                        pip_threshold = 0.01,
                        aux_pip_thres = 0.5,
                        aux_nbr_thres = 1,
                        fine_mapped_files = NULL,
                        full_matrix_path = NULL,
                        subset = NULL) {

  # Ensure output directories exist
  dir_sqlite <- file.path(output_dir, "SQLite")
  dir_output <- file.path(output_dir, "V_traits")
  dir.create(dir_sqlite, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

  # ─────────────────────────────────────────────────────
  # Case 1: Use precomputed full matrix
  if (!is.null(full_matrix_path)) {
    message("Using full precomputed pleiotropy matrix")
    full_mat <- fread(full_matrix_path)

    # Select a subset of provided auxiliary traits
    if (!is.null(subset)) {
        # Expand subset to include _pos and _neg columns
        trait_cols <- unlist(lapply(subset, function(tr) paste0(tr, c("_pos", "_neg"))))

        # Display message with first few traits
        message("Including auxiliary traits: ", paste(head(subset, 5), collapse = ", "),
                if (length(subset) > 5) ", ...")

        # Ensure we retain RSID and A1 for orientation
        keep <- c("RSID", "A1", trait_cols)

        # Subset the matrix
        full_mat <- full_mat[, ..keep]
    }

    # if custom focal trait fine-mapping file is provided
    db_path <- file.path(output_dir, "SQLite", sprintf("%s_SQLite.db", focal_trait))
    if (!file.exists(db_path)) {
      if (!is.null(fine_mapped_files) && focal_trait %in% names(fine_mapped_files)) {
        message("Converting custom fine-mapped file to SQLite for: ", focal_trait)

        # Validate
        jpepR::validate_finemap(fine_mapped_files[[focal_trait]])
                
        # Save
        finemap_file <- data.table::fread(fine_mapped_files[[focal_trait]])
        finemap_file <- finemap_file[PIP > pip_threshold]
        con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
        dbWriteTable(con, "trait_data", finemap_file, overwrite = TRUE)
        dbDisconnect(con)

      } else {
          stop("No SQLite found for `", focal_trait, "` and no custom file provided.")
      }
    }

    # Load focal trait separately
    focal_data <- jpepR::load_summary_stats(focal_trait, output_dir)
    focal_data <- focal_data[PIP > pip_threshold]
    focal_data <- unique(focal_data, by = "RSID")
    setkey(focal_data, RSID)
    setkey(full_mat, RSID)
    
     # Repair errors
    snps <- "rs11624006"
    snps2 <- "rs11711477"
    thetraits <- c("PASS_MSCV_Vuckovic2020", "PASS_RBC_Vuckovic2020", "PASS_MCV_Vuckovic2020", "PASS_MRV_Vuckovic2020", "PASS_MCH_Vuckovic2020")
    thetraits2 <- c("bmd_HEEL_TSCOREz")

    # Define trait column names
    traits_pos <- paste0(thetraits, "_pos")
    traits_neg <- paste0(thetraits, "_neg")
    traits_pos2 <- paste0(thetraits2, "_pos")
    traits_neg2 <- paste0(thetraits2, "_neg")

    # Swap for first SNP
    tmp_pos <- as.vector(full_mat[RSID %in% snps, ..traits_pos])
    tmp_neg <- as.vector(full_mat[RSID %in% snps, ..traits_neg])

    full_mat[RSID %in% snps, (traits_pos) := tmp_neg]
    full_mat[RSID %in% snps, (traits_neg) := tmp_pos]

    # Swap for second SNP
    tmp_pos2 <- as.vector(full_mat[RSID %in% snps2, ..traits_pos2])
    tmp_neg2 <- as.vector(full_mat[RSID %in% snps2, ..traits_neg2])

    full_mat[RSID %in% snps2, (traits_pos2) := tmp_neg2]
    full_mat[RSID %in% snps2, (traits_neg2) := tmp_pos2]    

    # Keep only SNPs in the focal trait
    merged <- merge(full_mat, focal_data[, .(CHR, BP, RSID, A1, A2, PIP, BETA)], by = "RSID")

    # Reweight all *_pos and *_neg columns using focal trait sign
    trait_cols <- grep("_(pos|neg)$", colnames(merged), value = TRUE)
    flip <- merged$A1.y != merged$A1.x
    merged[flip, BETA := -BETA]
    merged[flip, A1.y := A1.x]

    # Identify auxiliary traits before splitting into _pos/_neg
    aux_base <- unique(gsub("_(pos|neg)$", "", trait_cols))

    # Apply the per-row update rules
    for (trait in aux_base) {
      pos_col <- paste0(trait, "_pos")
      neg_col <- paste0(trait, "_neg")

      # If both columns exist in the merged data
      if (all(c(pos_col, neg_col) %in% colnames(merged))) {

        # Case 1: _pos > 0 & BETA > 0 → scale _pos
        idx0 <- merged[[pos_col]] > 0 & sign(merged$BETA) > 0
        
        # Case 2: _pos > 0 & BETA < 0 → flip to _neg
        idx1 <- merged[[pos_col]] > 0 & sign(merged$BETA) < 0
        
        # Case 3: _neg > 0 & BETA < 0 → scale _neg
        idx2 <- merged[[neg_col]] > 0 & sign(merged$BETA) < 0
        
        # Case 4: _neg > 0 & BETA > 0 → flip to _pos
        idx3 <- merged[[neg_col]] > 0 & sign(merged$BETA) > 0
        
        # Case 5: PIP below threshold → zero out both
        idx4 <- merged$PIP <= pip_threshold

        # Updates for all cases
        # Case 1:
        merged[idx0, paste0(pos_col, "") := get(pos_col) * PIP]
        
        # Case 2:
        merged[idx1, paste0(neg_col, "") := get(pos_col) * PIP]
        merged[idx1, paste0(pos_col, "") := 0]
        
        # Case 3:
        merged[idx2, paste0(pos_col, "") := get(neg_col) * PIP]
        merged[idx2, paste0(neg_col, "") := 0]

        # Case 4:
        merged[idx3, paste0(neg_col, "") := get(neg_col) * PIP]

        # Case 5: 
        merged[idx4, c(pos_col, neg_col) := 0]
      }
    }

    # Order columns
    trait_cols <- c(paste0(subset, "_pos"), paste0(subset, "_neg"))

    # Format output
    V_traits <- merged[, c("CHR", "BP", "RSID", "A1.y", "A2", "PIP", "BETA", trait_cols), with = FALSE]
    setnames(V_traits, "A1.y", "A1")

    # Save and return
    output_file <- file.path(dir_output, sprintf("%s_V_traits_thres%s.tsv.gz", focal_trait, pip_threshold))
    fwrite(V_traits, file = output_file, sep = "\t")
    message("✅ V_trait matrix saved to: ", output_file)
    return(output_file)
  }

  # ─────────────────────────────────────────────────────
  # Case 2: Build matrix from focal + auxiliary traits

  # ─────────────────────────────────────────────────────
  # Step 1: Convert fine-mapped files to SQLite format if needed
  all_traits <- unique(c(focal_trait, aux_traits))

  for (trait in all_traits) {
    db_path <- file.path(output_dir, "SQLite", sprintf("%s_SQLite.db", trait))
    if (!file.exists(db_path)) {
      if (!is.null(fine_mapped_files) && trait %in% names(fine_mapped_files)) {
        message("Converting custom fine-mapped file to SQLite for trait: ", trait)

        # Validate
        jpepR::validate_finemap(fine_mapped_files[[trait]])
        
        # Save
        finemap_file <- data.table::fread(fine_mapped_files[[trait]])
        finemap_file <- finemap_file[PIP > pip_threshold]
        con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
        dbWriteTable(con, "trait_data", finemap_file, overwrite = TRUE)
        dbDisconnect(con)

      } else {
        stop("No SQLite found for trait `", trait, "` and no custom file provided.")
      }
    }
  }

  # ─────────────────────────────────────────────────────
  # Step 2: Load focal trait
  focal_data <- jpepR::load_summary_stats(focal_trait, output_dir)
  focal_data <- focal_data[PIP > pip_threshold]
  focal_data <- unique(focal_data, by = "RSID")
  data.table::setorder(focal_data, RSID)

  V_traits <- focal_data[, .(CHR, BP, RSID, A1, A2, PIP, BETA)]

  # ─────────────────────────────────────────────────────
  # Step 3: Process auxiliary traits
  aux_traits <- setdiff(aux_traits, focal_trait)
  for (aux in aux_traits) {
    message("Processing auxiliary trait: ", aux)
    aux_data <- jpepR::load_summary_stats(aux, output_dir)[, .(RSID, A1, A2, PIP, BETA)]

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
  # Step 4: Filter SNPs and traits
  aux_cols <- setdiff(colnames(V_traits), c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA"))
  keep_rows <- rowSums(V_traits[, ..aux_cols]) > 0
  V_traits <- V_traits[keep_rows]

  filt_aux <- V_traits[, lapply(.SD, function(x) sum(x > aux_pip_thres)), .SDcols = aux_cols]
  keep_cols <- names(filt_aux)[filt_aux > aux_nbr_thres]
  keep_all <- c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA", keep_cols)
  V_traits <- V_traits[, ..keep_all]

  # ─────────────────────────────────────────────────────
  # Step 5: Save
  output_file <- file.path(dir_output, sprintf("%s_V_traits_thres%s.tsv.gz", focal_trait, pip_threshold))
  data.table::fwrite(V_traits, file = output_file, sep = "\t")
  message("✅ V_trait matrix saved to: ", output_file)

  return(output_file)
}
