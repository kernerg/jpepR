% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_jpep_model.R
\name{run_jpep_model}
\alias{run_jpep_model}
\title{Run J-PEP bNMF joint decomposition}
\usage{
run_jpep_model(
  model,
  V_trait,
  V_tissue,
  Tolerance = 1e-07,
  k = 15,
  NbrReps = 10,
  Alpha_proj = 10,
  Alpha_cross = 0,
  w_tis = 0.5,
  attempts = 30
)
}
\arguments{
\item{model}{Either JPEP, Pleiotropic or Epigenomic}

\item{V_trait}{Trait matrix}

\item{V_tissue}{Tissue matrix}

\item{Tolerance}{Convergence tolerance}

\item{k}{Number of clusters}

\item{NbrReps}{Number of repetitions}

\item{Alpha_proj}{Projection penalty}

\item{Alpha_cross}{Cross-correlation penalty}

\item{w_tis}{Weight on tissue component}
}
\description{
Run J-PEP bNMF joint decomposition
}
