\name{filterOutIterations}
\alias{filterOutIterations}
\title{
Filter Out Iterations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filterOutIterations(wall, hall, cnt_errors, cnt_reconstructed, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wall}{
%%     ~~Describe \code{wall} here~~
}
  \item{hall}{
%%     ~~Describe \code{hall} here~~
}
  \item{cnt_errors}{
%%     ~~Describe \code{cnt_errors} here~~
}
  \item{cnt_reconstructed}{
%%     ~~Describe \code{cnt_reconstructed} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
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
function (wall, hall, cnt_errors, cnt_reconstructed, params) 
{
    num_processesToExtract <- params$num_processesToExtract
    thresh_removeLastPercent <- params$thresh_removeLastPercent
    num_totIterations <- ncol(wall)/num_processesToExtract
    tot.rm.iterations <- round(thresh_removeLastPercent * num_totIterations)
    if (tot.rm.iterations > 0) {
        closeness.mutCounts <- matrix(0, nrow = num_totIterations, 
            ncol = 1)
        for (i in 1:num_totIterations) {
            closeness.mutCounts[i, ] <- base::norm(cnt_errors[[i]], 
                "F")
        }
        indexClosenessGenomes <- order(closeness.mutCounts, decreasing = TRUE)
        removeIterations <- indexClosenessGenomes[1:tot.rm.iterations]
        removeIterationSets <- matrix(0, nrow = (num_processesToExtract * 
            tot.rm.iterations), ncol = 1)
        for (i in 1:tot.rm.iterations) {
            iStart <- num_processesToExtract * (removeIterations[i] - 
                1) + 1
            iEnd <- num_processesToExtract * removeIterations[i]
            tmpRowRange <- (num_processesToExtract * (i - 1) + 
                1):(num_processesToExtract * i)
            removeIterationSets[tmpRowRange, ] <- iStart:iEnd
        }
        wall <- wall[, -removeIterationSets]
        hall <- hall[-removeIterationSets, ]
        cnt_errors <- cnt_errors[-removeIterations]
        cnt_reconstructed <- cnt_reconstructed[-removeIterations]
    }
    res.list <- list()
    res.list$Wall <- wall
    res.list$Hall <- hall
    res.list$mutCounts.errors <- cnt_errors
    res.list$mutCounts.reconstructed <- cnt_reconstructed
    return(res.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
