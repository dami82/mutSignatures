\name{resolveMutSignatures}
\alias{resolveMutSignatures}
\title{
Resolve Mutation Signatures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
resolveMutSignatures(mutCountData, signFreqData, byFreq = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutCountData}{
%%     ~~Describe \code{mutCountData} here~~
}
  \item{signFreqData}{
%%     ~~Describe \code{signFreqData} here~~
}
  \item{byFreq}{
%%     ~~Describe \code{byFreq} here~~
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
function (mutCountData, signFreqData, byFreq = TRUE) 
{
    if (class(mutCountData) == "mutationCounts") 
        mutCountData <- as.data.frame(mutCountData)
    if (class(signFreqData) == "mutationSignatures") 
        signFreqData <- as.data.frame(signFreqData)
    if (!(sum(!rownames(mutCountData) \%in\% rownames(signFreqData)) == 
        0 & sum(!rownames(signFreqData) \%in\% rownames(mutCountData)) == 
        0)) {
        stop(paste("There is an issue with the mutTypes.", "MutTypes in the mutType Count Matrix", 
            "have to match those in the signature", "Matrix... check rownames()"))
    }
    mutCountData <- mutCountData[rownames(signFreqData), ]
    if (byFreq) {
        full.Y <- apply(mutCountData, 2, (function(clmn) {
            1000 * clmn/sum(clmn)
        }))
    }
    else {
        full.Y <- as.matrix(mutCountData)
    }
    mutSums <- apply(mutCountData, 2, sum)
    my.signatures <- apply(signFreqData, 2, (function(clmn) {
        clmn/sum(clmn)
    }))
    X <- as.matrix(my.signatures)
    out <- list()
    res <- NMF::fcnnls(x = X, y = full.Y)
    beta.hat <- data.frame(t(res$x/ifelse(byFreq, 1000, 1)), 
        stringsAsFactors = FALSE)
    colnames(beta.hat) <- colnames(signFreqData)
    rownames(beta.hat) <- colnames(mutCountData)
    out$results <- list()
    out$coeffs <- list()
    out$coeffs$beta.hat <- beta.hat
    out$coeffs$unexplained.mutNum <- round((1 - apply(beta.hat, 
        1, sum)) * mutSums, digits = 0)
    out$coeffs$unexplained.mutFrac <- out$coeffs$unexplained.mutNum/mutSums
    if (byFreq) {
        for (i in 1:nrow(beta.hat)) {
            beta.hat[i, ] <- beta.hat[i, ] * mutSums[i]
        }
    }
    out$results$count.result <- as.mutation.exposure(beta.hat, 
        samplesAsCols = FALSE)
    out$results$freq.result <- as.mutation.exposure(do.call(rbind, 
        lapply(1:nrow(beta.hat), (function(jjj) {
            beta.hat[jjj, ]/sum(beta.hat[jjj, ])
        }))), samplesAsCols = FALSE)
    out$results$fitted <- res$fitted
    out$results$residuals <- res$residuals
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
