\name{prelimProcessAssess}
\alias{prelimProcessAssess}
\title{
Assess Preliminary Processes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prelimProcessAssess(input, maxProcess = 6, approach = "counts", plot = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{maxProcess}{
%%     ~~Describe \code{maxProcess} here~~
}
  \item{approach}{
%%     ~~Describe \code{approach} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (input, maxProcess = 6, approach = "counts", plot = TRUE, 
    verbose = TRUE) 
{
    tmpParams <- setMutClusterParams(num_processesToExtract = 1, 
        approach = approach, stopconv = 10000, niter = 1e+05, 
        thresh_removeWeakMutTypes = 0.01, thresh_removeLastPercent = 0.025, 
        num_totIterations = 10, num_parallelCores = 1, debug = FALSE, 
        logIterations = "full")
    tmpParams <- as.list(tmpParams)
    input_mutCounts <- as.matrix(input)
    if (tmpParams$approach != "counts") {
        tmpFRQ.input <- frequencize(input_mutCounts)
        tmpParams$approach <- "counts"
        input_mutCounts <- tmpFRQ.input$freqs
    }
    bckgrnd.removed.mutCounts <- removeWeak(input_mutCounts, 
        tmpParams)
    bckgrnd.removed.mutset <- bckgrnd.removed.mutCounts$removed.mutset
    bckgrnd.removed.mutCounts <- bckgrnd.removed.mutCounts$output.mutCounts
    total.mutationTypes <- nrow(bckgrnd.removed.mutCounts)
    total.samples <- ncol(bckgrnd.removed.mutCounts)
    tmpParams$num_totIterations <- 1
    medianMaxErr <- sum(t(sapply(1:nrow(bckgrnd.removed.mutCounts), 
        (function(i) {
            ((median(bckgrnd.removed.mutCounts[i, ]) - bckgrnd.removed.mutCounts[i, 
                ]))^2
        }))))
    if (verbose) 
        message("Preliminary Mutational Process Assessment: ", 
            appendLF = FALSE)
    outRes <- lapply(1:maxProcess, (function(i) {
        tmpParams$num_processesToExtract <- i
        tmpRes <- suppressWarnings(suppressMessages(extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts, 
            params = tmpParams, bootStrap = FALSE)))
        tmpErr <- sum((bckgrnd.removed.mutCounts - tmpRes$mutCounts.reconstructed)^2)
        if (verbose) 
            message(".", appendLF = FALSE)
        tmpErr
    }))
    if (verbose) 
        message("", appendLF = TRUE)
    err.points <- c(medianMaxErr, do.call(c, outRes))
    err.points <- (-1) * (err.points - max(err.points))
    err.points <- err.points/max(err.points, na.rm = TRUE)
    if (plot) {
        plot(err.points ~ c(0:(length(err.points) - 1)), type = "n", 
            las = 1, axes = FALSE, ylab = "", xlab = "", main = "Preliminary Mutational Process Assessment")
        axis(side = 1, at = NULL, cex = 0.65, font = 3)
        axis(side = 2, at = seq(0, 1, by = 0.2), cex = 0.65, 
            font = 3, labels = format(seq(1, 0, by = -0.2), digits = 2, 
                nsmall = 2), las = 1)
        lines(c(0:(length(err.points) - 1)), err.points, lwd = 1.5, 
            col = "red2")
        points(c(0:(length(err.points) - 1)), err.points, pch = 19, 
            col = "gray30")
        title(xlab = "Num of Signatures", line = 2.2, cex.lab = 1.2, 
            font = 2)
        title(ylab = "Error (\% vs. Median)", line = 3.1, cex.lab = 1.2, 
            font = 2)
        box()
    }
    return(data.frame(numProcess = c(0:(length(err.points) - 
        1)), percentErr = (1 - err.points), stringsAsFactors = TRUE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
