\name{evaluateStability}
\alias{evaluateStability}
\title{
Evaluate Stability
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evaluateStability(wall, hall, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wall}{
%%     ~~Describe \code{wall} here~~
}
  \item{hall}{
%%     ~~Describe \code{hall} here~~
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
function (wall, hall, params) 
{
    BIG_NUMBER <- 100
    CONVERG_ITER <- 10
    CONVERG_CUTOFF <- 0.005
    TOTAL_INIT_CONDITIONS <- 5
    num_processesToExtract <- params$num_processesToExtract
    num_totReplicates <- params$num_totReplicates
    distanceFunction <- params$distanceFunction
    minClusterDist <- BIG_NUMBER
    totalIter <- ncol(wall)/num_processesToExtract
    idx = matrix(0, nrow = nrow(hall), ncol = 1)
    clusterCompactness <- matrix(0, nrow = num_processesToExtract, 
        ncol = totalIter)
    iStartDataSet = seq(1, ncol(wall), by = num_processesToExtract)
    iStartingDataSet = iStartDataSet[sample(1:totalIter)]
    for (iInitData in 1:min(c(TOTAL_INIT_CONDITIONS, totalIter))) {
        iStartingData <- iStartingDataSet[iInitData]
        iEnd <- iStartingData + num_processesToExtract - 1
        centroids <- cbind(wall[, iStartingData:iEnd])
        centroidsTest <- sapply(1:ncol(centroids), (function(kk) {
            runif(nrow(centroids))
        }))
        countIRep <- 0
        for (iRep in 1:num_totReplicates) {
            tmp.tab <- t(cbind(centroids, wall))
            tmp.pdist <- as.vector(proxy::dist(tmp.tab, distanceFunction))
            if (num_processesToExtract > 1) 
                tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
            allDist <- pracma::squareform(tmp.pdist)
            cd.colRange <- (ncol(centroids) + 1):ncol(allDist)
            centroidDist = t(allDist[1:ncol(centroids), cd.colRange])
            jRange <- sort(1:num_processesToExtract)
            for (jIndex in 1:num_processesToExtract) {
                j <- jRange[jIndex]
                for (i in seq(1, ncol(wall), by = num_processesToExtract)) {
                  iRange = i:(i + num_processesToExtract - 1)
                  tmp.min <- min(centroidDist[iRange, j], na.rm = TRUE)
                  Ind <- which(centroidDist[iRange, j] == tmp.min)[1]
                  centroidDist[iRange[Ind], ] <- BIG_NUMBER
                  idx[iRange[Ind], 1] <- j
                }
            }
            maxDistToNewCentroids <- 0
            for (i in 1:ncol(centroids)) {
                tmp.dset <- wall[, as.vector(idx == i)]
                centroids[, i] <- apply(tmp.dset, 1, mean)
                tmp.dset <- t(cbind(centroids[, i], centroidsTest[, 
                  i]))
                tmp.pdist <- as.vector(proxy::dist(tmp.dset, 
                  distanceFunction))
                tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
                maxDistToNewCentroids <- max(maxDistToNewCentroids, 
                  tmp.pdist, na.rm = TRUE)
            }
            if (maxDistToNewCentroids < CONVERG_CUTOFF) {
                countIRep <- countIRep + 1
            }
            else {
                countIRep <- 0
                centroidsTest <- centroids
            }
            if (countIRep == CONVERG_ITER) {
                break
            }
        }
        for (i in 1:ncol(centroids)) {
            tmp.tab <- t(cbind(centroids[, i], wall[, as.vector(idx == 
                i)]))
            tmp.pdist <- as.vector(proxy::dist(tmp.tab, distanceFunction))
            tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
            clusterDist <- pracma::squareform(tmp.pdist)
            clusterCompactness[i, ] = clusterDist[1, 2:ncol(clusterDist)]
        }
        dist.test <- apply(clusterCompactness, 2, (function(clm) {
            mean(clm, na.rm = TRUE)
        }))
        if (sum(minClusterDist > dist.test) == length(dist.test)) {
            centroidsFinal <- centroids
            idxFinal <- idx
            clusterCompactnessFinal <- clusterCompactness
        }
    }
    centroids <- t(centroidsFinal)
    idx <- idxFinal
    clusterCompactness <- clusterCompactnessFinal
    centDist <- apply(clusterCompactness, 1, (function(tmprw) {
        mean(tmprw, na.rm = TRUE)
    }))
    centDistInd <- order(centDist, decreasing = FALSE)
    clusterCompactness <- clusterCompactness[centDistInd, ]
    centroids <- centroids[centDistInd, ]
    idxNew <- idx
    for (i in 1:num_processesToExtract) {
        idxNew[as.vector(idx == centDistInd[i]), 1] <- i
    }
    idx <- idxNew
    if (num_processesToExtract > 1) {
        processStab <- silhouetteMLB(data = t(wall), fac = idx, 
            distanceFunction)
        processStabAvg <- matrix(0, nrow = 1, ncol = num_processesToExtract)
        for (i in 1:num_processesToExtract) {
            processStabAvg[1, i] = mean(processStab[idx == i])
        }
    }
    else {
        tmp.tab <- t(cbind(centroids, wall))
        tmp.pdist <- as.vector(proxy::dist(tmp.tab, distanceFunction))
        tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
        allDist <- pracma::squareform(tmp.pdist)
        processStab <- 1 - t(allDist[1, 2:ncol(allDist)])
        xrange <- c(min(processStab), max(processStab))
        xrange[1] <- ifelse(xrange[1] > 0, 0, (-1.15) * abs(xrange[1]))
        xrange[2] <- 1.15
        graphics::barplot(sort(as.numeric(processStab), decreasing = TRUE), 
            col = "gray20", xlim = xrange, horiz = TRUE, xlab = "Silhouette Value", 
            ylab = "", main = "Silhouette Plot", border = "gray20")
        graphics::abline(v = 0)
        graphics::title(ylab = "Iter. Results (by Group)", line = 1, 
            cex.lab = 1, font = 2)
        processStabAvg <- apply(processStab, 2, (function(clmn) {
            base::mean(clmn, na.rm = TRUE)
        }))
    }
    if (num_processesToExtract > 1) {
        centroidStd <- matrix(0, nrow = nrow(centroids), ncol = ncol(centroids))
    }
    else {
        centroidStd <- matrix(0, ncol = length(centroids), nrow = 1)
    }
    for (i in 1:num_processesToExtract) {
        centroidStd[i, ] <- apply(wall[, idx == i], 1, (function(rw) {
            stats::sd(rw, na.rm = TRUE)
        }))
    }
    centroids <- t(cbind(centroids))
    centroidStd <- t(centroidStd)
    idxS <- matrix(0, nrow = length(idx), ncol = 1)
    for (i in seq(1, ncol(wall), by = num_processesToExtract)) {
        iEnd <- i + num_processesToExtract - 1
        idxG <- idx[i:iEnd]
        for (j in 1:num_processesToExtract) {
            idxS[(i + j - 1), ] = which(idxG == j)
        }
    }
    exposure <- matrix(0, nrow = max(idxS), ncol(hall))
    exposureStd <- matrix(0, nrow = max(idxS), ncol(hall))
    for (i in 1:max(idxS)) {
        exposure[i, ] <- apply(hall[idx == i, ], 2, (function(cl) {
            mean(cl, na.rm = TRUE)
        }))
        exposureStd[i, ] <- apply(hall[idx == i, ], 2, (function(cl) {
            sd(cl, na.rm = TRUE)
        }))
    }
    if (num_processesToExtract < 2) {
        centroids <- t(centroids)
    }
    result.list <- list()
    result.list$centroids <- centroids
    result.list$centroidStd <- centroidStd
    result.list$exposure <- exposure
    result.list$exposureStd <- exposureStd
    result.list$idx <- idx
    result.list$idxS <- idxS
    result.list$processStab <- processStab
    result.list$processStabAvg <- processStabAvg
    result.list$clusterCompactness <- clusterCompactness
    return(result.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
