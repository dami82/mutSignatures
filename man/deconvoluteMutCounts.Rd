\name{deconvoluteMutCounts}
\alias{deconvoluteMutCounts}
\title{
Deconvolute Mutation Counts
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deconvoluteMutCounts(input_mutCounts, params)
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
    num_totIterations <- params$num_totIterations
    num_processesToExtract <- params$num_processesToExtract
    distanceFunction <- params$distanceFunction
    thresh_removeWeakMutTypes <- params$thresh_removeWeakMutTypes
    num_parallelCores <- params$num_parallelCores
    guided <- params$guided
    num_totReplicates <- params$num_totReplicates
    thresh_removeLastPercent <- params$thresh_removeLastPercent
    colnames(input_mutCounts) <- NULL
    rownames(input_mutCounts) <- NULL
    input_mutCounts <- as.matrix(input_mutCounts)
    bckgrnd.removed.mutCounts <- removeWeak(input_mutCounts, 
        params)
    bckgrnd.removed.mutset <- bckgrnd.removed.mutCounts$removed.mutset
    bckgrnd.removed.mutCounts <- bckgrnd.removed.mutCounts$output.mutCounts
    total.mutationTypes <- nrow(bckgrnd.removed.mutCounts)
    total.samples <- ncol(bckgrnd.removed.mutCounts)
    if (guided) {
        guide.W <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts, 
            params = params, bootStrap = FALSE)
        guide.W <- guide.W$Wk
    }
    else {
        guide.W <- 0
    }
    if (num_parallelCores < 2) {
        muCounts.checkDF <- tryCatch(lapply(1:num_totIterations, 
            (function(j) {
                if (j \%in\% as.integer(seq(1, num_totIterations, 
                  length.out = 100))) {
                  message(paste("(", j, ")", sep = ""), appendLF = FALSE)
                }
                tmp.out <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts, 
                  bootStrap = TRUE, params = params)
                if (guided) {
                  re.ORD <- rep(0, num_processesToExtract)
                  for (ki in 1:num_processesToExtract) {
                    my.i <- order(apply(abs(tmp.out$Wk - guide.W[, 
                      ki]), 2, sum))
                    if (ki > 1) {
                      my.i[re.ORD[1:(ki - 1)]] <- max(my.i) + 
                        1
                    }
                    re.ORD[ki] <- which.min(my.i)
                  }
                }
                else {
                  re.ORD <- 1:num_processesToExtract
                }
                if (num_processesToExtract > 1) {
                  tmp.out$Wk <- tmp.out$Wk[, re.ORD]
                  tmp.out$Hk <- tmp.out$Hk[re.ORD, ]
                }
                else {
                  tmp.out$Wk <- tmp.out$Wk
                  tmp.out$Hk <- rbind(tmp.out$Hk)
                }
                tmp.out
            })), error = (function(e) {
            print(e)
        }))
        message("Done!", appendLF = TRUE)
    }
    else {
        max.cores <- parallel::detectCores()
        max.cores <- max.cores - 1
        max.cores <- ifelse(max.cores < 1, 1, max.cores)
        use.cores <- ifelse(1 <= num_parallelCores & num_parallelCores <= 
            max.cores, num_parallelCores, max.cores)
        cl <- suppressMessages(parallel::makeCluster(use.cores, 
            outfile = ""))
        print(paste("Extracting", num_processesToExtract, "mutational signatures X", 
            num_totIterations, "iterations using", use.cores, 
            "cores"))
        suppressMessages(doParallel::registerDoParallel(cl))
        stuffToExp <- c("alexaNMF", "leadZeros", "extractSignatures", 
            "frequencize", "bootstrapCancerGenomes", "chihJenNMF", 
            "params")
        suppressMessages(parallel::clusterExport(cl, stuffToExp))
        muCounts.checkDF <- tryCatch(foreach::foreach(j = (1:num_totIterations), 
            .verbose = TRUE, .packages = "stats") \%dopar\% {
            if (j \%in\% as.integer(seq(1, num_totIterations, length.out = 100))) {
                message(paste("(", j, ")", sep = ""), appendLF = FALSE)
            }
            tmp.out <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts, 
                params = params)
            if (guided) {
                re.ORD <- rep(0, num_processesToExtract)
                for (ki in 1:num_processesToExtract) {
                  my.i <- order(apply(abs(tmp.out$Wk - guide.W[, 
                    ki]), 2, sum))
                  if (ki > 1) {
                    my.i[re.ORD[1:(ki - 1)]] <- max(my.i) + 1
                  }
                  re.ORD[ki] <- which.min(my.i)
                }
            }
            else {
                re.ORD <- 1:num_processesToExtract
            }
            if (num_processesToExtract > 1) {
                tmp.out$Wk <- tmp.out$Wk[, re.ORD]
                tmp.out$Hk <- tmp.out$Hk[re.ORD, ]
            }
            else {
                tmp.out$Wk <- tmp.out$Wk
                tmp.out$Hk <- rbind(tmp.out$Hk)
            }
            tmp.out
        }, error = (function(e) {
            print(e)
        }), finally = (function(f) {
            parallel::stopCluster(cl)
        }))
        message("Done!", appendLF = TRUE)
    }
    W.all <- do.call(cbind, lapply(muCounts.checkDF, (function(tmp) {
        tmp$Wk
    })))
    H.all <- do.call(rbind, lapply(muCounts.checkDF, (function(tmp) {
        tmp$Hk
    })))
    errors.all <- lapply(muCounts.checkDF, (function(tmp) {
        tmp$mutCounts.errors
    }))
    reconstruct.all <- lapply(muCounts.checkDF, (function(tmp) {
        tmp$mutCounts.reconstructed
    }))
    fltr.mutCounts.data <- filterOutIterations(wall = W.all, 
        hall = H.all, cnt_errors = errors.all, cnt_reconstructed = reconstruct.all, 
        params)
    stability.check <- evaluateStability(wall = fltr.mutCounts.data$Wall, 
        hall = fltr.mutCounts.data$Hall, params)
    final.mutCounts.data <- addWeak(mutationTypesToAddSet = bckgrnd.removed.mutset, 
        processes_I = stability.check$centroids, processesStd_I = stability.check$centroidStd, 
        Wall_I = fltr.mutCounts.data$Wall, genomeErrors_I = fltr.mutCounts.data$mutCounts.errors, 
        genomesReconstructed_I = fltr.mutCounts.data$mutCounts.reconstructed)
    deconvoluted.results <- list()
    deconvoluted.results$Wall <- final.mutCounts.data$Wall
    deconvoluted.results$Hall <- fltr.mutCounts.data$Hall
    deconvoluted.results$mutCountErrors <- final.mutCounts.data$mutCountErrors
    deconvoluted.results$mutCountReconstructed <- final.mutCounts.data$mutCountReconstructed
    deconvoluted.results$idx <- stability.check$idx
    deconvoluted.results$idxS <- stability.check$idxS
    deconvoluted.results$processes <- final.mutCounts.data$processes
    deconvoluted.results$processesStd <- final.mutCounts.data$processesStd
    deconvoluted.results$exposure <- stability.check$exposure
    deconvoluted.results$exposureStd <- stability.check$exposureStd
    deconvoluted.results$processStab <- stability.check$processStab
    deconvoluted.results$processStabAvg <- stability.check$processStabAvg
    deconvoluted.results$clusterCompactness <- stability.check$clusterCompactness
    return(deconvoluted.results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
