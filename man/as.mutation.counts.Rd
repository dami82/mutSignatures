% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s02_all_methods.R
\name{as.mutation.counts}
\alias{as.mutation.counts}
\alias{as.mutation.counts,data.frame-method}
\alias{as.mutation.counts,data.frame,ANY,ANY-method}
\title{Method as.mutation.counts.}
\usage{
as.mutation.counts(x, rownames = NULL, colnames = NULL)

\S4method{as.mutation.counts}{data.frame}(x, rownames = NULL, colnames = NULL)
}
\arguments{
\item{x}{an object to extract Signature Identifiers from, i.e. a mutSignExposures-class}

\item{rownames}{character vector to overwrite data row names. Can be NULL if rownames(x) is not NULL.}

\item{colnames}{character vector to overwrite data column names. Can be NULL if colnames(x) is not NULL.}
}
\description{
Cast a data.frame into a mutationCounts-class object.
}
