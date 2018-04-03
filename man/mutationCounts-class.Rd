\name{mutationCounts-class}
\Rdversion{1.1}
\docType{class}
\alias{mutationCounts-class}
\alias{as.data.frame,mutationCounts-method}
\alias{as.matrix,mutationCounts-method}
\alias{[,mutationCounts,numeric,ANY,ANY-method}
\alias{coerceObj,mutationCounts,character-method}
\alias{getCounts,mutationCounts-method}
\alias{getMutationTypes,mutationCounts-method}
\alias{getSampleIdentifiers,mutationCounts-method}
\alias{initialize,mutationCounts-method}
\alias{plot,mutationCounts,ANY-method}
\alias{print,mutationCounts-method}
\alias{show,mutationCounts-method}

\title{Class \code{"mutationCounts"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mutationCounts", x, muts, samples)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{counts}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mutTypes}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{sampleId}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "mutationCounts")}: ... }
    \item{as.matrix}{\code{signature(x = "mutationCounts")}: ... }
    \item{[}{\code{signature(x = "mutationCounts", i = "numeric", j = "ANY", drop = "ANY")}: ... }
    \item{coerceObj}{\code{signature(x = "mutationCounts", to = "character")}: ... }
    \item{getCounts}{\code{signature(x = "mutationCounts")}: ... }
    \item{getMutationTypes}{\code{signature(x = "mutationCounts")}: ... }
    \item{getSampleIdentifiers}{\code{signature(x = "mutationCounts")}: ... }
    \item{initialize}{\code{signature(.Object = "mutationCounts")}: ... }
    \item{plot}{\code{signature(x = "mutationCounts", y = "ANY")}: ... }
    \item{print}{\code{signature(x = "mutationCounts")}: ... }
    \item{show}{\code{signature(object = "mutationCounts")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("mutationCounts")
}
\keyword{classes}
