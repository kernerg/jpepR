#' Construct V_trait matrix
#'
#' @param focal_trait Name of focal trait
#' @param aux_traits Vector of auxiliary trait names
#' @param output_dir Output directory
#' @param pip_threshold Minimum PIP threshold
#' @param aux_pip_thres Min PIP for aux trait
#' @param aux_nbr_thres Min number of traits with strong PIP per SNP
#' @export
make_vtrait <- function(focal_trait, aux_traits, output_dir,
                        pip_threshold = 0.01, aux_pip_thres = 0.5, aux_nbr_thres = 1) {

  script <- system.file("scripts", "make_vtrait_refactored.R", package = "jpepR")

  aux_traits_string <- paste(aux_traits, collapse = ",")
  cmd <- glue::glue(
    "Rscript {shQuote(script)} {shQuote(focal_trait)} {shQuote(aux_traits_string)} ",
    "{shQuote(pip_threshold)} {shQuote(aux_pip_thres)} {shQuote(aux_nbr_thres)} {shQuote(output_dir)}"
  )
  
  system(cmd)

}
