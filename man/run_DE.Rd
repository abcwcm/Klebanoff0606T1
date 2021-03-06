% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tgrams.R
\name{run_DE}
\alias{run_DE}
\title{findMarkers wrapper}
\usage{
run_DE(
  sce.object,
  group_identifier = "Sample",
  direction = "up",
  logfc = 0.2,
  FDR = 0.01,
  rank = Inf,
  comp_name = NULL
)
}
\arguments{
\item{sce.object}{SingleCellExperiment object}

\item{group_identifier}{define the variable that will be used to determine which
samples belong to the mutant or the wildtype condition, default: "Sample"}

\item{direction}{parameter for \code{scran::findMarkers}; default: "up"}

\item{logfc}{parameter for code{scran::findMarkers}}

\item{FDR}{parameter for \code{extract_markers}, default: 0.01}

\item{rank}{parameter for \code{extract_markers}, default: Inf}

\item{comp_name}{parameter for \code{extract_markers}; this should be
the clonotype for which the comparison is being done

#@return list with both the findMarkers results (which will consist of lists of 2)
and the extract_markers results, which will be a data.table with gene names,
but without the statistical metrics}
}
\description{
Wrapper function for running \code{scran::findMarkers}.
Will only return the results for the MUT condition, i.e. logFC values
represent MUT/WT.
}
