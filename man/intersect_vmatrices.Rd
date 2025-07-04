% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intersect_vmatrices.R
\name{intersect_vmatrices}
\alias{intersect_vmatrices}
\title{Merge V_trait and V_tissue into matched SNPs}
\usage{
intersect_vmatrices(focal_trait, output_dir, pip_threshold = 0.01)
}
\arguments{
\item{focal_trait}{Trait name}

\item{output_dir}{Output dir}

\item{pip_threshold}{Minimum PIP}
}
\value{
A list with V_trait and V_tissue matrices (SNPs × traits or tissues)
}
\description{
Merge V_trait and V_tissue into matched SNPs
}
