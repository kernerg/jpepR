#' Merge V_trait and V_tissue into matched SNPs
#'
#' @name intersect_vmatrices
#' @param focal_trait Trait name
#' @param output_dir Output dir
#' @param pip_threshold Minimum PIP
#'
#' @return A list with V_trait and V_tissue matrices (SNPs × traits or tissues)
#' @export 
intersect_vmatrices <- function(focal_trait, output_dir, pip_threshold = 0.01) {
  info_cols <- c("CHR", "BP", "RSID", "A1", "A2", "PIP", "BETA")
  
  vtrait_file <- file.path(output_dir, "V_traits", glue::glue("{focal_trait}_V_traits_thres{pip_threshold}.tsv.gz"))
  vtissue_file <- file.path(output_dir, "V_tissues", glue::glue("{focal_trait}_V_tissues_thres{pip_threshold}.tsv.gz"))

  # Load matrices
  vt <- data.table::fread(vtrait_file)
  vs <- data.table::fread(vtissue_file)

  # Merge on CHR and BP
  merge_cols <- c("CHR", "BP")
  merged <- merge(vt, vs, by = merge_cols)

  # Remove rows with all-zero V_trait values
  row_zero <- which(rowSums(merged[, setdiff(colnames(vt), info_cols), with = FALSE]) == 0)
  if (length(row_zero) > 0) {
    merged <- merged[-row_zero]
  }

  # Construct CHR labels like "chr1", "chr10", etc.
  chr_str <- paste0("chr", merged$CHR)

  # Final rownames: chrCHR;BP;BP
  row_ids <- sprintf("%s;%s;%s", chr_str, merged$BP - 1, merged$BP)

  # Apply lexicographic ordering on chrCHR and numeric ordering on BP
  order_idx <- order(chr_str, as.integer(merged$BP))

  # Order rownames and matrices
  row_ids <- row_ids[order_idx]
  V_trait <- as.matrix(merged[, setdiff(colnames(vt), info_cols), with = FALSE])[order_idx, , drop = FALSE]
  V_tissue <- as.matrix(merged[, setdiff(colnames(vs), info_cols), with = FALSE])[order_idx, , drop = FALSE]

  rownames(V_trait) <- row_ids
  rownames(V_tissue) <- row_ids

  list(
    V_trait = V_trait,
    V_tissue = V_tissue
  )
}
