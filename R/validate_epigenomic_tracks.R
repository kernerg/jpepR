#' Validate the format of epigenomic BED tracks
#'
#' Ensures that a custom `combined_bed` object resembles the EpiMap-like format expected
#' by the `make_vtissue()` pipeline.
#'
#' @param combined_bed A data.table or data.frame containing epigenomic regions
#'
#' @return Invisibly returns TRUE if valid, otherwise throws an informative error
#' @export
validate_epigenomic_tracks <- function(combined_bed) {
  required_cols <- c("CHR", "START", "STOP", "track")
  
  # Check for required columns
  if (!all(required_cols %in% names(combined_bed))) {
    stop("`combined_bed` must contain the following columns: ", paste(required_cols, collapse = ", "))
  }

  # CHR: should be numeric and in 1:22
  if (!is.numeric(combined_bed$CHR) || any(!combined_bed$CHR %in% 1:22)) {
    stop("`CHR` must be numeric and contain only values from 1 to 22.")
  }

  # START, STOP: positive integers, STOP >= START
  if (!is.numeric(combined_bed$START) || any(combined_bed$START < 0)) {
    stop("`START` must be a positive numeric value.")
  }
  if (!is.numeric(combined_bed$STOP) || any(combined_bed$STOP < combined_bed$START)) {
    stop("`STOP` must be numeric and greater than or equal to `START` for all rows.")
  }

  # track: character
  if (!is.character(combined_bed$track)) {
    stop("`track` must be a character vector indicating the track label.")
  }

  invisible(TRUE)
}
