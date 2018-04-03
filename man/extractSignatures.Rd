\name{extractSignatures}
\alias{extractSignatures}
\title{
Extract Signatures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractSignatures(mutCountMatrix, params, bootStrap = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutCountMatrix}{
%%     ~~Describe \code{mutCountMatrix} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{bootStrap}{
%%     ~~Describe \code{bootStrap} here~~
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
function (mutCountMatrix, params, bootStrap = TRUE) 
{
    num_processesToExtract <- params$num_processesToExtract
    approach <- params$approach
    algorithm <- params$algorithm
    eps <- params$eps
    if (bootStrap) {
        bstrpd.result <- bootstrapCancerGenomes(mutCountMatrix)
    }
    else {
        bstrpd.result <- mutCountMatrix
    }
    bstrpd.result[bstrpd.result < eps] <- eps
    if (algorithm \%in\% c("brunet", "alexa")) {
        nmf.results <- alexaNMF(v = bstrpd.result, r = num_processesToExtract, 
            params = params)
    }
    else {
        nmf.results <- chihJenNMF(v = bstrpd.result, r = num_processesToExtract, 
            params = params)
    }
    tmp.w <- nmf.results$w
    tmp.h <- nmf.results$h
    for (jj in 1:num_processesToExtract) {
        tmp.tot <- sum(tmp.w[, jj])
        tmp.w[, jj] <- tmp.w[, jj]/tmp.tot
        tmp.h[jj, ] <- tmp.h[jj, ] * tmp.tot
    }
    mutCountMatrix.reconstructed <- tmp.w \%*\% tmp.h
    result.list <- list()
    result.list$Wk <- tmp.w
    result.list$Hk <- tmp.h
    result.list$mutCounts.reconstructed <- mutCountMatrix.reconstructed
    result.list$mutCounts.errors <- bstrpd.result - mutCountMatrix.reconstructed
    if (params$logIterations != "lite") {
        result.list$inputMatrix <- bstrpd.result
        result.list$cosDist <- proxy::dist(rbind(as.vector(bstrpd.result), 
            as.vector(mutCountMatrix.reconstructed)), method = "cosine")[1]
    }
    return(result.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
