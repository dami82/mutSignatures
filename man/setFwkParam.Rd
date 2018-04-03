\name{setFwkParam}
\alias{setFwkParam}
\title{
Set or Adjust Single Framework Parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setFwkParam(x, label, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (x, label, value) 
{
    standardGeneric("setFwkParam")
  }, generic = structure("setFwkParam", package = "mutSignatures"), package = "mutSignatures", group = list(), valueClass = character(0), signature = c("x", 
"label", "value"), default = `\001NULL\001`, skeleton = (function (x, 
    label, value) 
stop("invalid call in method dispatch to 'setFwkParam' (no default method)", 
    domain = NA))(x, label, value), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
