\name{frequencize}
\alias{frequencize}
\title{
Convert Counts to Frequencies 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
frequencize(countMatrix, permille = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countMatrix}{
%%     ~~Describe \code{countMatrix} here~~
}
  \item{permille}{
%%     ~~Describe \code{permille} here~~
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
function (countMatrix, permille = TRUE) 
{
    out <- list()
    cf <- ifelse(permille, 1000, 1)
    out[["colSums"]] <- apply(countMatrix, 2, sum)
    out[["freqs"]] <- cf * apply(countMatrix, 2, (function(clmn) {
        clmn/sum(clmn)
    }))
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
