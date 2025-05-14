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
  
  vt <- data.table::fread(vtrait_file)
  vs <- data.table::fread(vtissue_file)
  merged <- merge(vt, vs, by = info_cols)
  
  list(
    V_trait = as.matrix(merged[, setdiff(colnames(vt), info_cols), with = FALSE]),
    V_tissue = as.matrix(merged[, setdiff(colnames(vs), info_cols), with = FALSE])
  )
}
