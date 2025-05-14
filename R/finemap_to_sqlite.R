#' Fine-map summary statistics and save to SQLite
#'
#' @param trait Trait name
#' @param sumstat_file Path to .sumstats.gz file
#' @param ld_blocks Path to LD block file
#' @param output_dir Output directory
#' @param pip_threshold Minimum PIP threshold
#' @param force Whether to overwrite existing SQLite file
#' @export
finemap_to_sqlite <- function(trait, sumstat_file, ld_blocks, output_dir, pip_threshold = 0.01, force = FALSE) {
  # internally call R script or embed logic here if standalone is desired
  script <- system.file("scripts", "gwas_finemap_to_sqlite.R", package = "jpepR")

  cmd <- glue::glue(
    "Rscript {shQuote(script)} {shQuote(trait)} {shQuote(ld_blocks)} ",
    "{shQuote(sumstat_file)} {shQuote(output_dir)} {shQuote(pip_threshold)} {shQuote(force)}"
  )

  system(cmd)

}
