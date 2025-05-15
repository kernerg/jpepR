#' Construct V_tissue matrix
#'
#' @param focal_trait Trait name
#' @param output_dir Output directory
#' @param epi_dir Path to epigenomic BED files or "default"
#' @param pip_threshold Minimum PIP
#' @export
make_vtissue <- function(focal_trait, output_dir, epi_dir = "default", pip_threshold = 0.01) {

  script <- system.file("scripts", "make_vtissue_refactored.R", package = "jpepR")

  cmd <- glue::glue(
    "Rscript {shQuote(script)} {shQuote(focal_trait)} {shQuote(epi_dir)} ",
    "{shQuote(pip_threshold)} {shQuote(output_dir)}"
  )

  system(cmd)

}
