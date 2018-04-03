\name{addWeak}
\alias{addWeak}
\title{
Add Weak Mutation Types
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addWeak(mutationTypesToAddSet, processes_I, processesStd_I, Wall_I, genomeErrors_I, genomesReconstructed_I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationTypesToAddSet}{
%%     ~~Describe \code{mutationTypesToAddSet} here~~
}
  \item{processes_I}{
%%     ~~Describe \code{processes_I} here~~
}
  \item{processesStd_I}{
%%     ~~Describe \code{processesStd_I} here~~
}
  \item{Wall_I}{
%%     ~~Describe \code{Wall_I} here~~
}
  \item{genomeErrors_I}{
%%     ~~Describe \code{genomeErrors_I} here~~
}
  \item{genomesReconstructed_I}{
%%     ~~Describe \code{genomesReconstructed_I} here~~
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
function (mutationTypesToAddSet, processes_I, processesStd_I, 
    Wall_I, genomeErrors_I, genomesReconstructed_I) 
{
    if (length(mutationTypesToAddSet) > 0 & mutationTypesToAddSet[1] > 
        0) {
        totalMutTypes <- nrow(Wall_I) + length(mutationTypesToAddSet)
        processes <- matrix(0, nrow = totalMutTypes, ncol = ncol(processes_I))
        processesStd <- matrix(0, nrow = totalMutTypes, ncol = ncol(processesStd_I))
        Wall <- matrix(0, nrow = totalMutTypes, ncol = ncol(Wall_I))
        genomeErrors <- lapply(1:length(genomeErrors_I), (function(i) {
            matrix(0, nrow = totalMutTypes, ncol = ncol(genomeErrors_I[[1]]))
        }))
        genomesReconstructed <- lapply(1:length(genomesReconstructed_I), 
            (function(i) {
                matrix(0, nrow = totalMutTypes, ncol = ncol(genomesReconstructed_I[[1]]))
            }))
        origArrayIndex <- 1
        for (i in 1:totalMutTypes) {
            if (!(i \%in\% mutationTypesToAddSet)) {
                processes[i, ] <- processes_I[origArrayIndex, 
                  ]
                processesStd[i, ] <- processesStd_I[origArrayIndex, 
                  ]
                Wall[i, ] <- Wall_I[origArrayIndex, ]
                for (j in 1:length(genomeErrors_I)) {
                  genomeErrors[[j]][i, ] <- genomeErrors_I[[j]][origArrayIndex, 
                    ]
                }
                for (j in 1:length(genomesReconstructed_I)) {
                  genomesReconstructed[[j]][i, ] <- genomesReconstructed_I[[j]][origArrayIndex, 
                    ]
                }
                origArrayIndex <- origArrayIndex + 1
            }
        }
    }
    else {
        processes <- processes_I
        processesStd <- processesStd_I
        Wall <- Wall_I
        genomeErrors <- genomeErrors_I
        genomesReconstructed <- genomesReconstructed_I
    }
    weakAdded.list <- list()
    weakAdded.list$processes <- processes
    weakAdded.list$processesStd <- processesStd
    weakAdded.list$Wall <- Wall
    weakAdded.list$mutCountErrors <- genomeErrors
    weakAdded.list$mutCountReconstructed <- genomesReconstructed
    return(weakAdded.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
