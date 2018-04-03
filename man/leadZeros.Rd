\name{leadZeros}
\alias{leadZeros}
\title{
Add Leading Zeros to Numbers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
leadZeros(n, m, char = "0", na.value = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{char}{
%%     ~~Describe \code{char} here~~
}
  \item{na.value}{
%%     ~~Describe \code{na.value} here~~
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
function (n, m, char = "0", na.value = NA) 
{
    max.zeros <- nchar(as.character(round(m)))
    tmp.digits <- nchar(as.character(round(n)))
    zeros.toAdd <- max.zeros - tmp.digits
    returnVect <- sapply(1:length(n), function(i) {
        if (zeros.toAdd[i] >= 0) {
            paste(c(rep(char, zeros.toAdd[i]), as.character(round(n[i]))), 
                sep = "", collapse = "")
        }
        else {
            na.value
        }
    })
    return(returnVect)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
