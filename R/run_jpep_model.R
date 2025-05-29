#' Run J-PEP bNMF joint decomposition
#'
#' @param model Either JPEP, Pleiotropic or Epigenomic
#' @param V_trait Trait matrix
#' @param V_tissue Tissue matrix
#' @param Tolerance Convergence tolerance
#' @param k Number of clusters
#' @param NbrReps Number of repetitions
#' @param Alpha_proj Projection penalty
#' @param Alpha_cross Cross-correlation penalty
#' @param w_tis Weight on tissue component
#' @export
run_jpep_model <- function(model, V_trait, V_tissue, Tolerance = 1e-7, k = 15, NbrReps = 10, Alpha_proj = 10, Alpha_cross = 0, w_tis = 0.5, attempts = 30) {
  suppressPackageStartupMessages({
    library(igraph)
  })
  source(system.file("scripts", "bNMF_functions.R", package = "jpepR"))
  run_single_bNMF_model(model_name = model, V_trait, V_tissue, Tolerance, k, NbrReps, Alpha_proj, Alpha_cross, w_tis, attempts)
}
