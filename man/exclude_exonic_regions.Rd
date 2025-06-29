% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/v_tissue_functions.R
\name{exclude_exonic_regions}
\alias{exclude_exonic_regions}
\title{Exclude SNPs overlapping exonic regions}
\usage{
exclude_exonic_regions(sumstats, exons)
}
\arguments{
\item{sumstats}{A data.table of summary statistics}

\item{exons}{A data.table of exon coordinates}
}
\value{
Filtered summary statistics data.table
}
\description{
Exclude SNPs overlapping exonic regions
}
