% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{fx.extract_exprs}
\alias{fx.extract_exprs}
\title{Extract expression values from SCE object}
\usage{
fx.extract_exprs(object, exprs_values, genes)
}
\arguments{
\item{object}{SCE object}

\item{exprs_values}{which type of \code{assay} values to extract}

\item{genes}{vector of gene names (subset of rownames)}
}
\value{
molten data.table with one row per gene and cell and the corresponding
\code{exprs_values} and their type in a third and fourth column.
}
\description{
Wrapper function for extracting expression values into a data.table
from an SCE object including checking the output.
}
\seealso{
\code{\link{make_long_dt}}
}
