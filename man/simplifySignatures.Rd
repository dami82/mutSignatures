% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s03_all_functions.R
\name{simplifySignatures}
\alias{simplifySignatures}
\title{Simplify Mutational Signatures.}
\usage{
simplifySignatures(x, asMutationSignatures = TRUE)
}
\arguments{
\item{x}{a mutationSignatures-class object}

\item{asMutationSignatures}{logical, shall the results be returned as a mutationSignatures-class object}
}
\value{
object including simplified mutational signatures data
}
\description{
This function is useufl when working with non-standard muation types, such as
tetra-nnucleotide mutation types or mutation types with long/complex context.
THe goal of this function is to aggregated together mutations that can be simplified
because of a common mutation core.
For example, mutation types AA[A>T]A, TA[A>T]A, CA[A>T]A, and GA[A>T]A can be simplified to the
core tri-nucleotide mutation A[A>T]A. THis function identifies mergeable mutation types,
and aggregates the corresponding counts/freqs.
}
\details{
This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
}
\examples{
A <- data.frame(Sig1=1:5, Sig2=5:1, Sig3=1:5)
A <- A/apply(A, 2, sum)
rownames(A) <- c("AA[C>A]A", "CA[C>A]A", "TA[C>A]A", "TA[C>G]A", "A[C>G]AT")
A <- mutSignatures::as.mutation.signatures(A)
mutSignatures::simplifySignatures(x = A)


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
