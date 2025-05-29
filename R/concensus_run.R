#' Retrieve concensus run
#'
#' @param model Either JPEP, Pleiotropic or Epigenomic
#' @param runs bNMF runs
#' @param V_trait Trait matrix
#' @param V_tissue Tissue matrix
#' @export
concensus_run <- function(model, runs, V_trait, V_tissue) {
  source(system.file("scripts", "bNMF_functions.R", package = "jpepR"))
  process_results(model_name = model, runs, V_trait, V_tissue)
}
