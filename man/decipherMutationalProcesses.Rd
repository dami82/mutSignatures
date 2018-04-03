\name{decipherMutationalProcesses}
\alias{decipherMutationalProcesses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decipher Mutational Processes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
decipherMutationalProcesses(input, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
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
function (input, params) 
{
    if (class(params) == "mutFrameworkParams" & class(input) == 
        "mutationCounts") {
        paramsList <- as.list(params)
        inputMAT <- as.matrix(input)
    }
    else {
        stop("Malformed Input")
    }
    if (paramsList$approach != "counts") {
        freq.input <- frequencize(inputMAT)
        inputMAT <- freq.input$freqs
        inputColsums <- freq.input$colSums
    }
    currentWarnings <- options()$warn
    options(warn = -1)
    if (is.numeric(paramsList$num_processesToExtract)) {
        paramsList$analyticApproach <- "denovo"
    }
    else {
        stop("An error occurred!")
    }
    if (paramsList$analyticApproach == "denovo") {
        deconvData <- deconvoluteMutCounts(input_mutCounts = inputMAT, 
            params = paramsList)
    }
    else {
        stop("An error occurred!")
    }
    mutProcesses <- list()
    final.proc <- data.frame(deconvData$processes, stringsAsFactors = FALSE)
    colnames(final.proc) <- paste("Sign.", sapply(1:ncol(final.proc), 
        (function(n) {
            leadZeros(n, (10 * ncol(final.proc)))
        })), sep = "")
    rownames(final.proc) <- getMutationTypes(input)
    signResult <- as.mutation.signatures(final.proc)
    final.expo <- data.frame(deconvData$exposure, stringsAsFactors = FALSE)
    if (getFwkParam(params, "approach") != "counts") {
        final.expo <- data.frame(sapply(1:ncol(final.expo), function(cjj) {
            inputColsums[cjj] * final.expo[, cjj]/sum(final.expo[, 
                cjj])
        }), row.names = NULL, stringsAsFactors = FALSE)
    }
    tryCatch({
        colnames(final.expo) <- getSampleIdentifiers(input)
    }, error = function(e) {
        NULL
    })
    rownames(final.expo) <- colnames(final.proc)
    expoResult <- as.mutation.exposure(x = final.expo)
    mutProcesses$Results <- list()
    mutProcesses$Results$signatures <- signResult
    mutProcesses$Results$exposures <- expoResult
    mutProcesses$RunSpecs <- list()
    mutProcesses$RunSpecs$input <- input
    mutProcesses$RunSpecs$params <- params
    if (getFwkParam(params, "logIterations") != "lite") {
        mutProcesses$Supplementary <- list()
        mutProcesses$Supplementary$allProcesses <- deconvData$Wall
        mutProcesses$Supplementary$allExposures <- deconvData$Hall
        mutProcesses$Supplementary$idx <- deconvData$idx
        mutProcesses$Supplementary$mutCountErrors <- deconvData$mutCountErrors
        mutProcesses$Supplementary$mutCountReconstructed <- deconvData$mutCountReconstructed
        mutProcesses$Supplementary$processStab <- deconvData$processStab
        mutProcesses$Supplementary$processStabAvg <- deconvData$processStabAvg
    }
    options(warn = currentWarnings)
    return(mutProcesses)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
