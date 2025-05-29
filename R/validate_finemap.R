#' Validate a fine-mapped file format
#'
#' Ensures that the custom fine-mapped file contains required columns and that
#' `PIP` values are properly defined between 0 and 1.
#'
#' @param file Path to a fine-mapped file (e.g. .tsv or .tsv.gz)
#'
#' @return Invisibly returns TRUE if valid, otherwise throws an error
#' @export
validate_finemap <- function(file) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }

  dat <- data.table::fread(file)

  required_cols <- c("CHR", "BP", "RSID", "A1", "A2", "BETA", "PIP")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  # Check CHR: numeric 1–22
  if (!is.numeric(dat$CHR) || any(!dat$CHR %in% 1:22)) {
    stop("Column `CHR` must be numeric and contain only values from 1 to 22.")
  }

  # Check BP: positive integer
  if (!is.numeric(dat$BP) || any(dat$BP <= 0)) {
    stop("Column `BP` must be a positive numeric value (hg19 position).")
  }

  # Check RSID: character or factor, no NA
  if (!is.character(dat$RSID) && !is.factor(dat$RSID)) {
    stop("Column `RSID` must be character or factor.")
  }
  if (any(is.na(dat$RSID))) {
    stop("Column `RSID` contains missing values.")
  }

  # Check BETA: numeric
  if (!is.numeric(dat$BETA)) {
    stop("Column `BETA` must be numeric.")
  }

  # Check PIP: numeric in (0, 1]
  if (!is.numeric(dat$PIP)) {
    stop("Column `PIP` must be numeric.")
  }
  if (any(dat$PIP < 0 | dat$PIP > 1, na.rm = TRUE)) {
    stop("Column `PIP` must contain values >= 0 and ≤ 1.")
  }

  invisible(TRUE)
}
