\name{removeWeak}
\alias{removeWeak}
\title{
Remove Uncommon Mutation Types 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
removeWeak(input_mutCounts, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_mutCounts}{
%%     ~~Describe \code{input_mutCounts} here~~
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
function (input_mutCounts, params) 
{
    thresh_removeWeakMutTypes <- params$thresh_removeWeakMutTypes
    sum.counts <- apply(input_mutCounts, 1, sum)
    sum.counts.idx <- order(sum.counts, decreasing = FALSE)
    sorted.sum.counts <- sum.counts[sum.counts.idx]
    tot.mut.counts <- sum(input_mutCounts)
    tot.muttypes.toremove <- sum((sapply(1:length(sorted.sum.counts), 
        (function(i) {
            sum(sorted.sum.counts[1:i])
        }))/tot.mut.counts) < thresh_removeWeakMutTypes)
    return.list <- list()
    if (tot.muttypes.toremove > 0) {
        removed.mutset <- sum.counts.idx[c(1:tot.muttypes.toremove)]
        input_mutCounts <- input_mutCounts[-removed.mutset, ]
        return.list$removed.mutset <- removed.mutset
    }
    else {
        return.list$removed.mutset <- (-1)
    }
    return.list$output.mutCounts <- input_mutCounts
    return(return.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
