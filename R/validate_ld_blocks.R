#' Validate custom LD blocks file or object
#'
#' Ensures that the LD block object has the correct format and values.
#'
#' @param ld_blocks A data.frame or data.table with LD blocks
#'
#' @return Invisibly returns TRUE if valid, otherwise throws an error
#' @export
validate_ld_blocks <- function(ld_blocks) {
  required_cols <- c("LOCUS", "CHR", "START", "STOP")

  # Check column names
  if (!all(required_cols %in% names(ld_blocks))) {
    stop("`ld_blocks` must contain the following columns: ", paste(required_cols, collapse = ", "))
  }

  # Check CHR
  if (!is.numeric(ld_blocks$CHR) || any(!ld_blocks$CHR %in% 1:22)) {
    stop("Column `CHR` must be numeric and contain only values 1 through 22.")
  }

  # Check START and STOP
  if (!is.numeric(ld_blocks$START) || !is.numeric(ld_blocks$STOP)) {
    stop("Columns `START` and `STOP` must be numeric hg19 genomic positions.")
  }
  if (any(ld_blocks$STOP <= ld_blocks$START)) {
    stop("Each LD block must have STOP > START.")
  }

  # Check LOCUS
  if (any(is.na(ld_blocks$LOCUS)) || anyDuplicated(ld_blocks$LOCUS)) {
    stop("Column `LOCUS` must contain unique, non-missing identifiers for each block.")
  }

  invisible(TRUE)
}
