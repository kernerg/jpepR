#' Validate summary statistics file for fine-mapping
#'
#' Ensures that summary statistics have the required columns with proper types and values.
#'
#' @param sumstats A data.frame or data.table containing GWAS summary statistics
#'
#' @return Invisibly returns TRUE if valid, otherwise throws an informative error
#' @export
validate_sumstats <- function(sumstats) {
  required_cols <- c("CHR", "BP", "RSID", "A1", "A2", "BETA", "P")

  # Check required columns
  if (!all(required_cols %in% names(sumstats))) {
    stop("`sumstats` must contain the following columns: ", paste(required_cols, collapse = ", "))
  }

  # CHR: numeric 1 to 22
  if (!is.numeric(sumstats$CHR) || any(!sumstats$CHR %in% 1:22)) {
    stop("Column `CHR` must be numeric and contain only values 1 through 22.")
  }

  # BP: positive integer
  if (!is.numeric(sumstats$BP) || any(sumstats$BP < 1)) {
    stop("Column `BP` must be a positive numeric hg19 base-pair position.")
  }

  # RSID: non-missing character (can be rsID or custom)
  if (!is.character(sumstats$RSID) && !is.factor(sumstats$RSID)) {
    stop("Column `RSID` must be character or factor.")
  }
  if (any(is.na(sumstats$RSID))) {
    stop("Column `RSID` contains missing values.")
  }

  # BETA: numeric
  if (!is.numeric(sumstats$BETA)) {
    stop("Column `BETA` must be numeric (GWAS effect size for A1).")
  }

  # P: numeric between 0 and 1 (exclusive)
  if (!is.numeric(sumstats$P) || any(sumstats$P <= 0 | sumstats$P >= 1, na.rm = TRUE)) {
    stop("Column `P` must be numeric with values in the open interval (0, 1).")
  }

  invisible(TRUE)
}
