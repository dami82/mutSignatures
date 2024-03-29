% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s03_all_functions.R
\name{matchSignatures}
\alias{matchSignatures}
\title{Match Mutational Signatures.}
\usage{
matchSignatures(
  mutSign,
  reference = NULL,
  method = "cosine",
  threshold = 0.5,
  plot = TRUE
)
}
\arguments{
\item{mutSign}{a mutationSignatures object}

\item{reference}{a mutationSignatures object. If NULL, COSMIC signatures will be retrieved}

\item{method}{distance method used to compute similarity (1 - distance)}

\item{threshold}{signal (similarity) upper threshold for maxing the signal}

\item{plot}{logical, shall a heatmap be plotted}
}
\value{
list, including distance matrix and a heatmap plot
}
\description{
Analyze the similarity between mutational signatures from different analyses/runs.
This function can be helpful to match de novo extracted signatures with previously
described signatures (such as COSMIC), or to reveal signatures that can be
identified with alternative NMF algorithms, or that may be due to an algorithm bias.
}
\details{
This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
}
\references{
More information and examples about mutational signature analysis can be found here:
\enumerate{
  \item \bold{Official website}: \url{http://www.mutsignatures.org}
  \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
  \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
 }
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
