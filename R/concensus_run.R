#' Retrieve concensus run
#'
#' @param runs bNMF runs
#' @param V_trait Trait matrix
#' @param V_tissue Tissue matrix
#' @export
concensus_run <- function(runs, V_trait, V_tissue) {
  source(system.file("scripts", "bNMF_functions.R", package = "jpepR"))
  process_results("JPEP", runs, V_trait, V_tissue)
}
