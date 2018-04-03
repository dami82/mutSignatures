\name{mutSignExposures-class}
\Rdversion{1.1}
\docType{class}
\alias{mutSignExposures-class}
\alias{as.data.frame,mutSignExposures-method}
\alias{[,mutSignExposures,numeric,ANY,ANY-method}
\alias{coerceObj,mutSignExposures,character-method}
\alias{getSampleIdentifiers,mutSignExposures-method}
\alias{getSignatureIdentifiers,mutSignExposures-method}
\alias{initialize,mutSignExposures-method}
\alias{plot,mutSignExposures,ANY-method}
\alias{print,mutSignExposures-method}
\alias{setSignatureIdentifiers,mutSignExposures,character-method}
\alias{show,mutSignExposures-method}

\title{Class \code{"mutSignExposures"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mutSignExposures", x, samples, signNames)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{exposures}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{sampleId}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{signatureId}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "mutSignExposures")}: ... }
    \item{[}{\code{signature(x = "mutSignExposures", i = "numeric", j = "ANY", drop = "ANY")}: ... }
    \item{coerceObj}{\code{signature(x = "mutSignExposures", to = "character")}: ... }
    \item{getSampleIdentifiers}{\code{signature(x = "mutSignExposures")}: ... }
    \item{getSignatureIdentifiers}{\code{signature(x = "mutSignExposures")}: ... }
    \item{initialize}{\code{signature(.Object = "mutSignExposures")}: ... }
    \item{plot}{\code{signature(x = "mutSignExposures", y = "ANY")}: ... }
    \item{print}{\code{signature(x = "mutSignExposures")}: ... }
    \item{setSignatureIdentifiers}{\code{signature(x = "mutSignExposures", names = "character")}: ... }
    \item{show}{\code{signature(object = "mutSignExposures")}: ... }
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
showClass("mutSignExposures")
}
\keyword{classes}
