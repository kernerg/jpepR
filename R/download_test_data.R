#' Download example test data for jpepR
#'
#' Downloads test files needed to run the vignette and example workflows.
#'
#' @param dest_dir Local directory to save downloaded files (default: "data-raw")
#' @param overwrite Logical. If TRUE, overwrite existing files. Default is FALSE.
#'
#' @return Invisibly returns the path to the downloaded folder
#' @export
download_test_data <- function(dest_dir = "data-raw", overwrite = FALSE) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  base_url <- "https://alkesgroup.broadinstitute.org/J-PEP/test_data"
  files <- c(
    "T2D_postfinemap.tsv",
    "T2D_auxtraits.rda",
    "big_pleio_matrix.tsv.gz",
    "Epimap_tracks.tsv.gz",
    "Epimap_bckg_expectation.rds",
    "Epimap_metadata.txt",
    "Epimap_sampleID.txt"
  )

  for (file in files) {
    dest_file <- file.path(dest_dir, file)
    if (!file.exists(dest_file) || overwrite) {
      message("Downloading: ", file)
      utils::download.file(
        url = paste0(base_url, "/", file),
        destfile = dest_file,
        mode = "wb", quiet = TRUE
      )
    } else {
      message("File already exists, skipping: ", file)
    }
  }

  invisible(normalizePath(dest_dir))
}
