\name{mutationSignatures-class}
\Rdversion{1.1}
\docType{class}
\alias{mutationSignatures-class}
\alias{as.data.frame,mutationSignatures-method}
\alias{as.list,mutationSignatures-method}
\alias{[,mutationSignatures,numeric,ANY,ANY-method}
\alias{cbind2,mutationSignatures,mutationSignatures-method}
\alias{coerceObj,mutationSignatures,character-method}
\alias{getMutationTypes,mutationSignatures-method}
\alias{getSignatureIdentifiers,mutationSignatures-method}
\alias{initialize,mutationSignatures-method}
\alias{plot,mutationSignatures,ANY-method}
\alias{print,mutationSignatures-method}
\alias{setSignatureIdentifiers,mutationSignatures,character-method}
\alias{show,mutationSignatures-method}

\title{Class \code{"mutationSignatures"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mutationSignatures", x, muts, signNames)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{mutationFreq}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mutTypes}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{signatureId}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "mutationSignatures")}: ... }
    \item{as.list}{\code{signature(x = "mutationSignatures")}: ... }
    \item{[}{\code{signature(x = "mutationSignatures", i = "numeric", j = "ANY", drop = "ANY")}: ... }
    \item{cbind2}{\code{signature(x = "mutationSignatures", y = "mutationSignatures")}: ... }
    \item{coerceObj}{\code{signature(x = "mutationSignatures", to = "character")}: ... }
    \item{getMutationTypes}{\code{signature(x = "mutationSignatures")}: ... }
    \item{getSignatureIdentifiers}{\code{signature(x = "mutationSignatures")}: ... }
    \item{initialize}{\code{signature(.Object = "mutationSignatures")}: ... }
    \item{plot}{\code{signature(x = "mutationSignatures", y = "ANY")}: ... }
    \item{print}{\code{signature(x = "mutationSignatures")}: ... }
    \item{setSignatureIdentifiers}{\code{signature(x = "mutationSignatures", names = "character")}: ... }
    \item{show}{\code{signature(object = "mutationSignatures")}: ... }
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
showClass("mutationSignatures")
}
\keyword{classes}
