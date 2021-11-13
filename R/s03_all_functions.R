#
##  ~~~~~~~~~~~~~~~~~~~~~~~~~
### ~~~~~ MutSignatures ~~~~~
##  ...all functions...
#   ~~~~~~~~~~~~~~~~~~~~~~~~~

# Damiano Fantini, Ph.D.
# 2021-Nov-13

###
##### Custom Functions - core package
###

#' Add Weak Mutation TYpes
#'
#' Restore Mutation Types that were initially excluded because a low number of total counts.
#'
#'
#' @param mutationTypesToAddSet Set of mutations to restore
#' @param processes_I Set of Mutational Processes
#' @param processesStd_I Set of standard deviations of all Mutational Processes
#' @param Wall_I Set of all W matrices previously extracted
#' @param genomeErrors_I Set of all residuals
#' @param genomesReconstructed_I Fitted Values according to the most likely Model
#'
#' @return Output is the final result of the deconvolution process
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @keywords internal
addWeak <- function (mutationTypesToAddSet,
                     processes_I,
                     processesStd_I,
                     Wall_I,
                     genomeErrors_I,
                     genomesReconstructed_I)
{
  if (length(mutationTypesToAddSet) > 0 & mutationTypesToAddSet[1] > 0) {
    totalMutTypes <- nrow(Wall_I) + length(mutationTypesToAddSet)
    processes <- matrix(0, nrow = totalMutTypes, ncol = ncol(processes_I))
    processesStd <- matrix(0, nrow = totalMutTypes, ncol = ncol(processesStd_I))
    Wall <- matrix(0, nrow = totalMutTypes, ncol = ncol(Wall_I))
    genomeErrors <- lapply(1:length(genomeErrors_I), (function(i) {
      matrix(0, nrow = totalMutTypes, ncol = ncol(genomeErrors_I[[1]]))
    }))
    genomesReconstructed <- lapply(1:length(genomesReconstructed_I),
                                   (function(i) {
                                     matrix(0,
                                            nrow = totalMutTypes,
                                            ncol = ncol(genomesReconstructed_I[[1]]))
                                   }))
    origArrayIndex <- 1
    for (i in 1:totalMutTypes) {
      #message(i)
      if (!(i %in% mutationTypesToAddSet)) {
        processes[i, ] <- processes_I[origArrayIndex, ]
        processesStd[i, ] <- processesStd_I[origArrayIndex,]
        Wall[i, ] <- Wall_I[origArrayIndex, ]
        for (j in 1:length(genomeErrors_I)) {
          genomeErrors[[j]][i, ] <- genomeErrors_I[[j]][origArrayIndex,  ]
        }
        for (j in 1:length(genomesReconstructed_I)) {
          genomesReconstructed[[j]][i, ] <- genomesReconstructed_I[[j]][origArrayIndex,   ]
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

#' Bootstrap a Mutation Count Matrix.
#'
#' Rearrange a Mutation count Matrix using the multivariate normal distribution.
#' The function returns a bootstrapped Mutation Count matrix whose dimensions
#' are identical to the input matrix.
#'
#' @param genomes a numeric matrix of Mutation Counts.
#' Rows correspond to Mutation Types, columns to different samples.
#' @param seed integer, set a seed to obtain reproducible results. Defaulted to NULL
#'
#' @return a numeric matrix of bootstrapped Mutation Counts.
#' Rows correspond to Mutation Types, columns to different samples.
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' x <- cbind(c(10, 100, 20, 200, 30, 5),
#'            c(100, 90, 80, 100, 11, 9))
#' mutSignatures:::bootstrapCancerGenomes(x)
#'
#'
#' @importFrom stats rmultinom
#'
#' @keywords internal
bootstrapCancerGenomes <- function (genomes, seed = NULL)
{
  if (!is.null(seed))
    try(set.seed(seed), silent = TRUE)

  genome.col.sums <- suppressWarnings(apply(genomes, 2, sum))
  my.denom <- matrix(genome.col.sums,
                     ncol = ncol(genomes),
                     nrow = nrow(genomes),
                     byrow = TRUE)
  norm.genomes <- suppressWarnings(genomes/my.denom)

  bootstrapGenomes <- suppressWarnings(
    lapply(1:ncol(genomes), (function(i) {
      tmp <- norm.genomes[,i]
      tmp <- stats::rmultinom(1, size = genome.col.sums[i], prob = tmp)
      as.numeric(tmp[,1])
    })))
  bootstrapGenomes <- suppressWarnings(do.call(cbind, bootstrapGenomes))
  return(bootstrapGenomes)
}

#' Evaluate Results Stability.
#'
#' Perform a final Stability check comparing the results from all iterations of the analysis.
#'
#' @param wall numeric matrix including the w results from all the iterations of the analysis
#' @param hall numeric matrix including the h results from all the iterations of the analysis
#' @param params list including all the parameters required for running tha analysis
#'
#' @return list including all results from the stability checks.
#' This includes the most likely signatures (cen-troids) and exposures.
#' All information for plotting the silhoueette plot will also be returned.
#'
#' @details The function evaluates the results from all iterations by performing
#' a silhouette check. A silhouette plot will also be plotted.
#' This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' # Obtain sample data
#' TMP <- mutSignatures:::getTestRunArgs("evaluateStability")
#' Y <- mutSignatures:::evaluateStability(wall = TMP$W,
#'                                        hall = TMP$H,
#'                                        params = TMP$params)
#'
#'
#' @importFrom stats runif sd
#' @importFrom proxy dist
#' @importFrom pracma squareform
#' @importFrom graphics barplot abline title
#'
#' @keywords internal
evaluateStability <- function (wall,
                               hall,
                               params)
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
  clusterCompactness <- matrix(0,
                               nrow = num_processesToExtract,
                               ncol = totalIter)
  iStartDataSet = seq(1, ncol(wall), by = num_processesToExtract)
  iStartingDataSet = iStartDataSet[sample(1:totalIter)]
  for (iInitData in 1:min(c(TOTAL_INIT_CONDITIONS, totalIter))) {
    #message(iInitData)
    iStartingData <- iStartingDataSet[iInitData]
    iEnd <- iStartingData + num_processesToExtract - 1
    centroids <- cbind(wall[, iStartingData:iEnd])
    centroidsTest <- sapply(1:ncol(centroids), (function(kk) {
      stats::runif(nrow(centroids))
    }))
    countIRep <- 0
    for (iRep in 1:num_totReplicates) {
      #message(iRep)
      tmp.tab <- base::t(cbind(centroids, wall))
      tmp.pdist <- base::as.vector(proxy::dist(tmp.tab, distanceFunction))
      if (num_processesToExtract > 1)
        tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
      allDist <- pracma::squareform(tmp.pdist)

      cd.colRange <- (ncol(centroids) + 1):ncol(allDist)
      centroidDist = t(allDist[1:ncol(centroids), cd.colRange])

      jRange <- sort(1:num_processesToExtract)
      for (jIndex in 1:num_processesToExtract) {
        j <- jRange[jIndex]
        for (i in seq(1, ncol(wall), by = num_processesToExtract)) {
          #message(i)
          iRange = i:(i + num_processesToExtract - 1)
          tmp.min <- min(centroidDist[iRange, j], na.rm = TRUE)
          Ind <- which(centroidDist[iRange, j] == tmp.min)[1]
          centroidDist[iRange[Ind], ] <- BIG_NUMBER
          idx[iRange[Ind], 1] <- j
        }
      }
      maxDistToNewCentroids <- 0
      for (i in 1:ncol(centroids)) {
        tmp.dset <- wall[, base::as.vector(idx == i)]
        centroids[, i] <- apply(tmp.dset, 1, mean)
        tmp.dset <- base::t(cbind(centroids[, i],
                                  centroidsTest[,i]))
        tmp.pdist <- base::as.vector(proxy::dist(tmp.dset, distanceFunction))
        tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
        maxDistToNewCentroids <- max(maxDistToNewCentroids,
                                     tmp.pdist,
                                     na.rm = TRUE)
      }
      if (maxDistToNewCentroids < CONVERG_CUTOFF) {
        countIRep <- countIRep + 1

      } else {
        countIRep <- 0
        centroidsTest <- centroids
      }

      if (countIRep == CONVERG_ITER) {
        break
      }
    }
    for (i in 1:ncol(centroids)) {
      tmp.tab <- base::t(cbind(centroids[, i],
                               wall[, base::as.vector(idx == i)]))
      tmp.pdist <- base::as.vector(proxy::dist(tmp.tab, distanceFunction))
      tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
      clusterDist <- pracma::squareform(tmp.pdist)
      clusterCompactness[i, ] = clusterDist[1, 2:ncol(clusterDist)]
    }
    dist.test <- apply(clusterCompactness, 2, (function(clm) {
      base::mean(clm, na.rm = TRUE)
    }))
    if (sum(minClusterDist > dist.test) == length(dist.test)) {
      centroidsFinal <- centroids
      idxFinal <- idx
      clusterCompactnessFinal <- clusterCompactness
    }
  }
  centroids <- base::t(centroidsFinal)
  idx <- idxFinal
  clusterCompactness <- clusterCompactnessFinal
  centDist <- apply(clusterCompactness, 1, (function(tmprw) {
    base::mean(tmprw, na.rm = TRUE)
  }))
  centDistInd <- base::order(centDist, decreasing = FALSE)
  clusterCompactness <- clusterCompactness[centDistInd, ]
  centroids <- centroids[centDistInd, ]
  idxNew <- idx
  for (i in 1:num_processesToExtract) {
    idxNew[base::as.vector(idx == centDistInd[i]), 1] <- i
  }
  idx <- idxNew
  if (num_processesToExtract > 1) {
    processStab <- silhouetteMLB(data = t(wall), fac = idx,
                                 distanceFunction)
    processStabAvg <- matrix(0, nrow = 1, ncol = num_processesToExtract)
    for (i in 1:num_processesToExtract) {
      processStabAvg[1, i] = base::mean(processStab[idx == i])
    }
  } else {
    # Adjusted params for a 1-class silhouette!
    tmp.tab <- base::t(cbind(centroids, wall))
    tmp.pdist <- base::as.vector(proxy::dist(tmp.tab, distanceFunction))
    tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
    allDist <- pracma::squareform(tmp.pdist)
    processStab <- 1 - base::t(allDist[1, 2:ncol(allDist)])

    # Silhouette plot
    xrange <- c(min(processStab), max(processStab))
    xrange[1] <- ifelse(xrange[1] > 0, 0, (-1.15) * abs(xrange[1]))
    xrange[2] <- 1.15
    graphics::barplot(base::sort(base::as.numeric(processStab), decreasing = TRUE),
                      col = "gray20",
                      xlim = xrange,
                      horiz = TRUE, xlab = "Silhouette Value",
                      ylab = "",
                      main = "Silhouette Plot",
                      border = "gray20")
    graphics::abline(v=0)
    graphics::title(ylab="Iter. Results (by Group)", line=1, cex.lab=1, font = 2)

    processStabAvg <- apply(processStab, 2, (function(clmn) {
      base::mean(clmn, na.rm = TRUE)
    }))
  }

  if (num_processesToExtract > 1) {
    centroidStd <- matrix(0, nrow = nrow(centroids), ncol = ncol(centroids))
  } else {
    centroidStd <- matrix(0, ncol = length(centroids), nrow = 1)
  }

  for (i in 1:num_processesToExtract) {
    centroidStd[i, ] <- apply(wall[, idx == i], 1, (function(rw) {
      stats::sd(rw, na.rm = TRUE)
    }))
  }
  centroids <- base::t(cbind(centroids))
  centroidStd <- base::t(centroidStd)
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
      base::mean(cl, na.rm = TRUE)
    }))
    exposureStd[i, ] <- apply(hall[idx == i, ], 2, (function(cl) {
      sd(cl, na.rm = TRUE)
    }))
  }

  # Fix to sign.to.extract.num = 1
  if (num_processesToExtract < 2){
    centroids <- base::t(centroids)
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

#' Remove Iterations that Generated Outlier Results.
#'
#' Internal function from the WTSI framework, ported to R. This is a core
#' function called from within a deconvoluteMutCounts() call. This
#' function removes iterations that generated outlier results.
#'
#' @param wall numeric matrix combining w results from all iterations
#' @param hall numeric matrix combining h results from all iterations
#' @param cnt_errors numeric matrix combining all residuals from all iterations
#' @param cnt_reconstructed numeric matrix combining fitted values from all iterations
#' @param params list including alll parameters for running the analysis
#'
#' @return list including all data required for running the subsequent stability check
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item Mutational signatures operative in bladder cancer (\bold{Oncogene paper}): \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#'
#' @keywords internal
filterOutIterations <- function (wall,
                                 hall,
                                 cnt_errors,
                                 cnt_reconstructed,
                                 params)
{
  num_processesToExtract <- params$num_processesToExtract
  thresh_removeLastPercent <- params$thresh_removeLastPercent
  num_totIterations <- ncol(wall)/num_processesToExtract
  tot.rm.iterations <- round(thresh_removeLastPercent * num_totIterations)
  if (tot.rm.iterations > 0) {
    closeness.mutCounts <- matrix(0,
                                  nrow = num_totIterations,
                                  ncol = 1)
    for (i in 1:num_totIterations) {
      closeness.mutCounts[i, ] <- base::norm(cnt_errors[[i]], "F")
    }
    indexClosenessGenomes <- order(closeness.mutCounts,
                                   decreasing = TRUE)
    removeIterations <- indexClosenessGenomes[1:tot.rm.iterations]
    removeIterationSets <- matrix(0,
                                  nrow = (num_processesToExtract * tot.rm.iterations),
                                  ncol = 1)
    for (i in 1:tot.rm.iterations) {
      iStart <- num_processesToExtract * (removeIterations[i] - 1) + 1
      iEnd <- num_processesToExtract * removeIterations[i]
      tmpRowRange <- (num_processesToExtract * (i - 1) + 1):(num_processesToExtract * i)
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

#' Generate Arguments for Running Examples and Mock Runs.
#'
#' This function generates objects that can be used for
#' running the examples included in the package documentation files,
#' as well as some simple mutSignature analyses. Note that his function is not exported.
#'
#' @param testN string, name of the function that we want to test
#'
#' @return an object (typically, a list) including sample data to run
#' analyses or to test mutSignatures functions.
#'
#' @details This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @keywords internal
getTestRunArgs <- function (testN = "evaluateStability")
{
  out <- list()

  if (testN == "evaluateStability") {
    out$W <- do.call(cbind, lapply(1:20, (function(i) {
      cbind(c(runif(4, 0.05, 0.15),
              c(1e-15 * runif(1,1, 9)),
              runif(3, 0.003, 0.007),
              runif(9, 0.04, 0.09)),
            c(runif(3, 0.08, 0.18),
              c(1e-15 * runif(1, 1, 9)),
              runif(1, 0.03, 0.07),
              runif(3, 0.004, 0.009),
              runif(9, 0.04, 0.09)))
    })))
    out$H <- do.call(rbind, lapply(1:40, (function(i) {
      c(runif(1, 0.005, 0.099),
        runif(2, 50, 800),
        runif(2, 0.005, 0.099),
        runif(1, 50, 750),
        runif(2, 0.005, 0.099),
        runif(1, 20, 800))
    })))
    out$params$num_processesToExtract <- 2
    out$params$num_totReplicates <- 20
    out$params$distanceFunction <- "cosine"
    out$params$analyticApproach <- "denovo"

  } else if (testN == "removeWeak") {

    out <- list()
    out$data <- sapply(1:12, function(i) {sample(50:99, size = 50, replace = TRUE) })
    out$data[6, ] <- sample(0:5, size = 12, replace = TRUE)
    out$params <- as.list(setMutClusterParams())

  } else if (testN == "silhouetteMLB") {

    out <- list()
    out$data <- rbind(sapply(1:10, function(i) {sample(x = 0:15, size = 10)}),
                      sapply(1:10, function(i) {sample(x = 11:29, size = 10)}))
    out$fac <- factor(c(rep(1, 10), rep(2, 10)))


  } else if (testN %in% c("alexaNMF", "chihJenNMF")) {

    out <- list()

    zz1 <- sapply(1:4, function(i) {base::sample(30:58, replace = TRUE, size = 20)})
    zz2 <- sapply(1:2, function(i) {base::sample(52:80, replace = TRUE, size = 20)})
    zz3 <- sapply(1:4, function(i) {base::sample(60:120, replace = TRUE, size = 15)})
    zz4 <- sapply(1:2, function(i) {base::sample(15:40, replace = TRUE, size = 15)})

    out$v <- rbind(
      cbind(zz1, zz2),
      cbind(zz3, zz4))

    out$r <- 2
    out$params <- as.list(setMutClusterParams())

  } else if (testN %in% c("decipherMutationalProcesses", "deconvoluteMutCounts",
                          "extractSignatures")) {

    out <- list()

    muts <- cbind(
      sapply(1:12, function(i) {c(sample(5:15, size = 20, replace = TRUE),
                                  sample(10:25, size = 10))}),
      sapply(1:8, function(i) {c(sample(10:25, size = 20, replace = TRUE),
                                 sample(5:15, size = 10))}))
    rownames(muts) <- c("A[A>C]A", "C[A>C]A", "G[A>C]A", "T[A>C]A", "T[A>C]C", "T[A>C]G",
                        "A[A>G]A", "C[A>G]A", "G[A>G]A", "T[A>G]A", "T[A>G]C", "T[A>G]G",
                        "A[A>T]A", "C[A>T]A", "G[A>T]A", "T[A>T]A", "T[A>T]C", "T[A>T]G",
                        "A[G>C]A", "C[G>C]A", "G[G>C]A", "T[G>C]A", "T[G>C]C", "T[G>C]G",
                        "A[G>T]A", "C[G>T]A", "G[G>T]A", "T[G>T]A", "T[G>T]C", "T[G>T]G")
    colnames(muts) <- paste0("SMPL.", 1001:1020)
    out$muts <- mutSignatures::as.mutation.counts(as.data.frame(muts))
    out$params <- setMutClusterParams(2)
    out$params@params$num_totIterations <- 10
    out$params@params$niter <- 1000

  } else if (testN == "custom_fcnnls") {

    out <- list()
    out$muts <- cbind(A = c(10, 15, 23, 15, 23,  5, 2, 6, 8),
                      B = c(12, 13, 22, 14, 18,  1, 5, 2, 7),
                      C = c(12,  5, 10,  5, 13, 22, 1, 1, 9))
    out$signs <- cbind(S1 = c(0.16, 0.06, 0.13, 0.06, 0.18, 0.28, 0.01, 0.01, 0.11),
                       S2 = c(0.09, 0.14, 0.22, 0.15, 0.21, 0.04, 0.02, 0.06, 0.07))

  } else if (testN == "custom_cssls") {

    CtC <- cbind(c(0.1728, 0.1179), c(0.1179, 0.1532))
    CtA <- cbind(c(10.76, 14.37), c(13.33, 9.05))
    Pset <- cbind(c(FALSE, TRUE), c(TRUE, FALSE))

    out <- list(CtC = CtC,
                CtA = CtA,
                Pset = Pset)

  } else if (testN == "extractXvarlinkData") {

    out <- c("getma.org/?cm=var&var=hg19,9,107576738,C,T&fts=all",
             "",
             "getma.org/?cm=var&var=hg19,9,107594970,G,A&fts=all",
             "getma.org/?cm=var&var=hg19,9,107599275,T,A&fts=all",
             "getma.org/?cm=var&var=hg19,9,107591257,C,T&fts=all",
             "",
             "getma.org/?cm=var&var=hg19,16,2336768,G,A&fts=all",
             "getma.org/?cm=var&var=hg19,16,2347484,G,C&fts=all",
             "getma.org/?cm=var&var=hg19,16,2339472,A,G&fts=all",
             "getma.org/?cm=var&var=hg19,16,2358451,C,T&fts=all",
             "getma.org/?cm=var&var=hg19,1,94522197,G,T&fts=all",
             "getma.org/?cm=var&var=hg19,1,94548924,G,A&fts=all")

  } else if (testN == "removeMismatchMut") {

    out <- data.frame(CHROM = c("chr1", "chr1", "chr2", "chr3", "chr3", "chr5"),
                      POS = c(12144, 155464, 4232, 35222, 35663, 244425),
                      REF = c("A", "T", "G", "G", "T", "G"),
                      ALT = c("G", "G", "T", "A", "C", "A"),
                      INFO = c(".", ".", ".", "ref133121", ".", "."),
                      context = c("TAA", "ATA", "CGG", "CCG", "ATA", "AAT"),
                      stringsAsFactors = FALSE)

  } else if (testN == "resolveMutSignatures") {
    out <- list()

    muts <- cbind(
      sapply(1:12, function(i) {c(sample(5:15, size = 20, replace = TRUE),
                                  sample(10:25, size = 10))}),
      sapply(1:8, function(i) {c(sample(10:25, size = 20, replace = TRUE),
                                 sample(5:15, size = 10))}))
    rownames(muts) <- c("A[A>C]A", "C[A>C]A", "G[A>C]A", "T[A>C]A", "T[A>C]C", "T[A>C]G",
                        "A[A>G]A", "C[A>G]A", "G[A>G]A", "T[A>G]A", "T[A>G]C", "T[A>G]G",
                        "A[A>T]A", "C[A>T]A", "G[A>T]A", "T[A>T]A", "T[A>T]C", "T[A>T]G",
                        "A[G>C]A", "C[G>C]A", "G[G>C]A", "T[G>C]A", "T[G>C]C", "T[G>C]G",
                        "A[G>T]A", "C[G>T]A", "G[G>T]A", "T[G>T]A", "T[G>T]C", "T[G>T]G")
    colnames(muts) <- paste0("SMPL.", 1001:1020)

    sigs <- cbind(
      SIG.1 = c(sample(5:15, size = 20, replace = TRUE), sample(10:25, size = 10)),
      SIG.2 = c(sample(10:25, size = 20, replace = TRUE),sample(5:15, size = 10)))
    for ( i in 1:ncol(sigs)) {
      sigs[, i] <- sigs[, i] / sum(sigs[, i])
    }

    rownames(sigs) <- c("A[A>C]A", "C[A>C]A", "G[A>C]A", "T[A>C]A", "T[A>C]C", "T[A>C]G",
                        "A[A>G]A", "C[A>G]A", "G[A>G]A", "T[A>G]A", "T[A>G]C", "T[A>G]G",
                        "A[A>T]A", "C[A>T]A", "G[A>T]A", "T[A>T]A", "T[A>T]C", "T[A>T]G",
                        "A[G>C]A", "C[G>C]A", "G[G>C]A", "T[G>C]A", "T[G>C]C", "T[G>C]G",
                        "A[G>T]A", "C[G>T]A", "G[G>T]A", "T[G>T]A", "T[G>T]C", "T[G>T]G")

    out$muts <- mutSignatures::as.mutation.counts(as.data.frame(muts))
    out$sigs <- mutSignatures::as.mutation.signatures(as.data.frame(sigs))


  } else if (testN == "countMutTypes") {

    out <- data.frame(CHROM = "chr1",
                      mutation = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
                                   "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
                                   "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                                   "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T",
                                   "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T"),
                      sample = paste0("S.", c(rep(1, 10), rep(2, 10))))


  } else if (testN == "filterSNV") {

    out <- data.frame(CHR = "chr1",
                      POS = c(972116, 1647600, 11529902, 16624394, 21617326,
                              26923997, 30718807, 36309634, 44805336, 99288016,
                              108181979, 146846852, 150067867, 152843497, 154054711,
                              158292407, 160895402, 161042822, 165667287, 166620791),
                      REF = c("G","C", "T", "G", "G", "A", "C", "G", "C", "G",
                              "C", "G", "T", "C", "A", "A", "G", "C", "C", "T"),
                      ALT = c("A","G", "C", "GC", "A", "T", "G", "T", "A", "T",
                              "A", "C", "G", "G", "G", "TT", "T", "G", "G", "A"),
                      SAMPLEID = paste0("smpl.", c(rep(1, 10), rep(2, 10))),
                      stringsAsFactors = FALSE)

  } else  {
    out <- NULL
  }
  return(out)
}

#' Add Leading Zeros to Numbers.
#'
#' Internal function to convert a numeric vector into a character vector, where all elements
#' have the same number of characters (nchar). This is obtained by pasting a series of leading zeros
#' (or other character) to each number in the input vector.
#'
#' @param n numeric vector whose numbers are to be transformed
#' @param m maximum number that will be used to define how many leading zeros to attach
#' @param char string (typically, a single character). This character is used to fill the leading space.
#' Defaults to 0.
#' @param na.value value used to fill mising values. Defaults to NA
#'
#' @return numeric vector of length equal to length(n), where all numbers are converted
#' to character and modified by attaching the required number of leading zeros (characters).
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' n = c(5:12, NA, 9:11)
#' m = 111
#' mutSignatures:::leadZeros(n=n, m=m)
#'
#' @keywords internal
leadZeros <- function (n, m, char = "0", na.value = NA)
{
  max.zeros <- nchar(base::as.character(round(m)))
  tmp.digits <- nchar(base::as.character(round(n)))
  zeros.toAdd <- max.zeros - tmp.digits

  returnVect <- sapply(1:length(n), function(i){
    if (is.na(zeros.toAdd[i])) {
      na.value
    } else if (zeros.toAdd[i] >= 0) {
      paste(c(rep(char, zeros.toAdd[i]), base::as.character(round(n[i]))), sep = "", collapse = "")
    } else {
      na.value
    }
  })
  return(returnVect)
}


#' Remove Mutation Types Not Meeting the Threshold.
#'
#' Remove mutation types that account for a total number of mutations below a defined threshold.
#'
#' @param input_mutCounts numeric matrix of Mutation Counts
#' @param params object (list) including all parameters required for running the analysis
#'
#' @return List including two elements:
#' \enumerate{
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("removeWeak")
#' nrow(x$data)
#' y <- mutSignatures:::removeWeak(input_mutCounts = x$data, params = x$params)
#' nrow(y)
#'
#' @keywords internal
removeWeak <- function (input_mutCounts,
                        params)
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

#' Silhouette Analysis.
#'
#' Analyze the clustering quality and generate a Silhouette Plot.
#'
#' @param data numeric matrix
#' @param fac clustering factor
#' @param method method to be used as distance function. Defaults to c("cosine")
#' @param plot logical, shall a barplot showing the cluster silhouettes be printed
#'
#' @return numeric vector including the silhouette values of the data poointts in the input matrix
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item \bold{Silhouette analysis in R}: \url{http://www.biotechworld.it/bioinf/2017/01/20/translating-matlabs-silhouette-function-to-r/}
#'  }
#'
#' @importFrom cluster silhouette
#' @importFrom proxy dist
#' @importFrom graphics barplot abline title
#'
#' @examples
#' library(mutSignatures)
#' x <- mutSignatures:::getTestRunArgs("silhouetteMLB")
#' y <- silhouetteMLB(data = x$data, fac = x$fac)
#' y
#'
#' @export
silhouetteMLB <- function (data,
                           fac,
                           method = "cosine",
                           plot = TRUE)
{
  if (nrow(data) != length(fac))
    stop("Bad input!")
  dist.matrix <- base::as.matrix(proxy::dist(x = data, method = method))
  sil.check <- cluster::silhouette(x = base::as.numeric(base::as.factor(fac)),
                                   dist = dist.matrix)
  if (plot == TRUE) {
    tmp <- lapply(unique(sil.check[, 1]), (function(clid) {
      part.out <- sil.check[sil.check[, 1] == clid, ]
      part.out[order(part.out[, 3], decreasing = TRUE),
      ]
    }))
    tmp <- do.call(rbind, tmp)
    xrange <- c(min(tmp[,3]), max(tmp[,3]))
    xrange[1] <- ifelse(xrange[1] > 0, 0, (-1.15) * abs(xrange[1]))
    xrange[2] <- 1.15
    graphics::barplot(tmp[nrow(tmp):1, 3],
                      col = base::as.factor(tmp[nrow(tmp):1,1]),
                      xlim = xrange,
                      horiz = TRUE, xlab = "Silhouette Value",
                      ylab = "",
                      main = "Silhouette Plot",
                      border = base::as.factor(tmp[nrow(tmp):1,1]))
    graphics::abline(v=0)
    graphics::title(ylab="Iter. Results (by Group)", line=1, cex.lab=1, font = 2)
  }
  return(base::as.vector(sil.check[, 3]))
}

#' Perform Non-negative Matrix Factorization using Brunet's Algotithm.
#'
#' Perform Non-negative Matrix Factorization.
#'
#' @param v numeric matrix of Mutation Type Counts
#' @param r numeric, number of signatures to extract
#' @param params list including all paramaters for running the analysis
#'
#' @return list including all paramaters for running the analysis:
#' \enumerate{
#'   \item \bold{W} extracted signatures
#'   \item \bold{H} contribution of each signature in all the samples of the input mut count matrix
#' }
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#'
#' @importFrom graphics plot points
#' @importFrom proxy dist
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("alexaNMF")
#' y <- mutSignatures:::alexaNMF(v = x$v, r = x$r, params = x$params)
#' y$w[1:5, ]
#'
#' @keywords internal
alexaNMF <- function(v, r, params)
{
  # Brunet algorithm, from Alexandrov NMF WTSI code
  # define params (hard-set)
  debug <- params$debug
  chk.step <- 50
  dot.eachSteps <- 2000

  # retireve user-defined params
  eps <- params$eps
  num.processes <- r
  stopconv <- params$stopconv
  niter <- params$niter
  err.threshold <- 1e-10
  stopRule <- ifelse(params$stopRule == "LA", "LA", "DF")

  # set.seed(231082)
  v <- base::as.matrix(v)
  rownames(v) <- NULL
  colnames(v) <- NULL

  # Double check input matrix
  if (min(v) < 0)
    stop("Matrix entries cannot be negative")
  if (min(apply(v, 1, sum)) == 0)
    stop("Entries cannot all be equal to 0")

  # Initialize W0 and H0
  W.k <- do.call(cbind, lapply(1:num.processes, (function(i){
    out.o <- runif(n = nrow(v), min = eps, max = 100)
    out.o/sum(out.o)
  })))
  H.k <- matrix((1/num.processes), nrow = num.processes, ncol = ncol(v))

  # Initialize looping vars
  itr <- 1
  chk.j <- 1
  stationary.chk <- 0
  force.out <- 1

  # Debugging plot
  if (debug)
    graphics::plot(-10, xlim = c(1000, niter),
                   ylim = c( ifelse(stopRule == "DF", (0.1*err.threshold), 0.001) ,
                             ifelse(stopRule == "DF", 10, ncol(H.k))), log = "xy",
                   xlab = "Iteration", ylab = "Variation", main = "Convergence")

  # Initialize the objects for comparing dissimilarity
  # DF approach
  cons.old <- base::as.vector(W.k)
  # Alexandrov approach
  consold <- matrix(0, nrow = ncol(H.k), ncol = ncol(H.k))

  while (itr < niter){

    if (itr %% dot.eachSteps == 0) {
      if (stationary.chk > chk.step) {
        message(":", appendLF = FALSE)
      } else {
        message(".", appendLF = FALSE)
      }
    }

    delta.01 <- apply(W.k, 2, sum)
    H.k <- H.k * (base::t(W.k) %*% (v/(W.k %*% H.k))) / delta.01
    H.k[H.k < eps] <- eps

    W.tmp <- W.k * ((v/(W.k %*% H.k)) %*% base::t(H.k))
    W.k <- do.call(cbind, lapply(1:ncol(W.tmp), (function(ci){
      W.tmp[,ci] / sum(H.k[ci,])
    })))
    W.k[W.k<eps] <- eps

    # check convergence every 'chk.step' iterations
    if (itr > stopconv & itr %% chk.step == 0 & stopRule == "DF") {
      chk.j <- chk.j + 1
      H.k[H.k < eps] <- eps
      W.k[W.k < eps] <- eps
      cons <- base::as.vector(W.k)

      # compare to consold and reorder
      dist.measure <- proxy::dist(rbind(cons, cons.old), method = "cosine") [1]
      cons.old <- cons

      if (debug)
        points(itr, (dist.measure + (err.threshold*0.1)), pch = 19, col = "red2")

      # evaluate distance
      if (dist.measure < err.threshold) {
        stationary.chk <- stationary.chk + 1
      } else {
        stationary.chk <- 0
      }
      if (stationary.chk > (stopconv / chk.step)) {
        force.out <- 0
        message("$", appendLF = FALSE)
        break()
      }
    } else if (itr > stopconv & itr %% chk.step == 0 & stopRule == "LA") {
      chk.j <- chk.j + 1
      H.k[H.k < eps] <- eps
      W.k[W.k < eps] <- eps
      y <- apply(H.k, 2, max)
      index <- apply(H.k, 2, (function(dt) {
        which.max(dt)[1]
      }))
      mat1 = base::t(sapply(1:ncol(H.k), (function(ii) {
        index
      })))
      mat2 = sapply(1:ncol(H.k), (function(ii) {
        index
      }))
      cons <- mat1 == mat2
      if (sum(cons != consold) == 0) {
        stationary.chk <- stationary.chk + 1
      } else {
        stationary.chk <- 0
      }
      #
      consold <- cons

      if (debug)
        points(itr, (sum(cons != consold) / 100), pch = 19, col = "red2")

      if (stationary.chk > (stopconv/chk.step)) {
        force.out <- 0
        message("$", appendLF = FALSE)
        break()
      }
    }
    itr <- itr + 1
  }
  if (force.out == 1) {
    message("!", appendLF = FALSE)
  }
  output <- list()
  output$w <- W.k
  output$h <- H.k
  return(output)
}


#' Perform Non-negative Matrix Factorization using Chih-Jen Lin's Algotithm.
#'
#' Perform Non-negative Matrix Factorization (alternative approach).
#'
#' @param v numeric matrix of Mutation Type Counts
#' @param r numeric, number of signatures to extract
#' @param params list including all paramaters for running the analysis
#'
#' @return list including all paramaters for running the analysis:
#' \enumerate{
#'   \item \bold{W} extracted signatures
#'   \item \bold{H} contribution of each signature in all the samples of the input mut count matrix
#' }
#'
#' @details This is a core internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item \bold{Chih-Jen Lin original paper}: \url{http://ieeexplore.ieee.org/document/4359171/}
#'  }
#'
#'
#' @importFrom graphics plot points
#' @importFrom proxy dist
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("chihJenNMF")
#' y <- mutSignatures:::chihJenNMF(v = x$v, r = x$r, params = x$params)
#' y$w[1:10, ]
#'
#' @keywords internal
chihJenNMF <- function(v, r, params) {

  # http://ieeexplore.ieee.org/document/4359171/
  # alternative approach for NMF
  # define params (hard-set)
  debug <- params$debug
  chk.step <- 50
  dot.eachSteps <- 2000

  # retireve user-defined params
  eps <- params$eps
  num.processes <- r
  stopconv <- params$stopconv
  niter <- params$niter
  err.threshold <- 1e-10

  # set.seed(231082)
  v <- base::as.matrix(v)
  rownames(v) <- NULL
  colnames(v) <- NULL

  # Double check input matrix
  if (min(v) < 0)
    stop("Matrix entries cannot be negative")
  if (min(apply(v, 1, sum)) == 0)
    stop("Entries cannot all be equal to 0")

  # Debugging plot
  if (debug)
    graphics::plot(-10, xlim = c(1000, niter), ylim = c( (0.1*err.threshold), 10), log = "xy",
                   xlab = "Iteration", ylab = "Variation", main = "Convergence")

  # Initialize W0 and H0
  W.k <- do.call(cbind, lapply(1:num.processes, (function(i){
    out.o <- runif(n = nrow(v), min = eps, max = 100)
    out.o/sum(out.o)
  })))
  H.k <- matrix((1/num.processes), nrow = num.processes, ncol = ncol(v))

  # Initialize looping vars
  itr <- 1
  chk.j <- 1
  cons.old <- base::as.vector(W.k)
  stationary.chk <- 0
  force.out <- 1
  while (itr < niter){

    if (itr %% dot.eachSteps == 0) {
      if (stationary.chk > chk.step) {
        message(":", appendLF = FALSE)
      } else {
        message(".", appendLF = FALSE)
      }
    }

    WtW <- base::t(W.k) %*% W.k
    gradH <- ((WtW %*% H.k) - (base::t(W.k) %*% v))
    H.b <- H.k; H.b[H.b<eps] <- eps;
    H.k <- H.k - (H.b / ((WtW %*% H.b) + eps)) * gradH
    HHt <- H.k %*% base::t(H.k)
    gradW <- (W.k %*% HHt) - (v %*% base::t(H.k))

    W.b <- W.k; W.b[W.b < eps] <- eps
    W.k <- W.k - (W.b / ((W.b %*% HHt) + eps)) * gradW
    S <- apply(W.k, 2, sum)

    W.k <- do.call(cbind, lapply(1:ncol(W.k), (function(ci){
      W.k[,ci] / S[ci]
    })))
    H.k <- do.call(rbind, lapply(1:nrow(H.k), (function(ri){
      H.k[ri,] * S[ri]
    })))

    # optional ? Keep as is for now
    H.k <- do.call(cbind, lapply(1:ncol(H.k), (function(ci){
      H.k[,ci] / sum(H.k[,ci])
    })))

    # Final non-negative check
    H.k[H.k < eps] <- eps
    W.k[W.k < eps] <- eps

    # check convergence every 'chk.step' iterations
    if (itr > stopconv & itr %% chk.step == 0) {
      chk.j <- chk.j + 1
      W.k[W.k < eps] <- eps
      cons <- base::as.vector(W.k)
      # compare to consold and reorder
      dist.measure <- proxy::dist(rbind(cons, cons.old), method = "cosine") [1]
      cons.old <- cons

      if (debug)
        points(itr, (dist.measure + (err.threshold*0.1)), pch = 19, col = "red2")

      # evaluate distance
      if (dist.measure < err.threshold) {
        stationary.chk <- stationary.chk + 1
      } else {
        stationary.chk <- 0
      }
      if (stationary.chk > (stopconv / chk.step)) {
        force.out <- 0
        message("$", appendLF = FALSE)
        break()
      }
    }
    itr <- itr + 1
  }
  if (force.out == 1) {
    message("!", appendLF = FALSE)
  }
  output <- list()
  output$w <- W.k
  output$h <- H.k
  return(output)
}

#' Decipher Mutational Processes Contributing to a Collection of Genomic Mutations.
#'
#' Decipher Mutational ProCancer cells accumulate DNA mutations as result of DNA
#' damage and DNA repair processes. Thiscomputational framework allows to decipher
#' mutational processes from cancer-derived somatic mutational catalogs.
#'
#' @param input a mutationCounts-class object, including a mutation counts data.
#' @param params a mutFrameworkParams-class object including all the parameters required for
#' running the mutational signature analysis.
#'
#' @return list including all results of the analysis.
#' The extracted signatures (processes) are included in the "processes" element of the list.
#' The relative contribution of each signature in each sample is summarized in the
#' "exposures" element.
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is the main user interface for the mutSignatures analysis.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' library(mutSignatures)
#' x <- mutSignatures:::getTestRunArgs("decipherMutationalProcesses")
#' x$muts
#' y <- mutSignatures::decipherMutationalProcesses(input = x$muts,
#'                                                 params = x$params)
#' y$Results$signatures
#'
#' @export
decipherMutationalProcesses <- function (input,
                                         params)
{
  if (class(params) == "mutFrameworkParams" & class(input) == "mutationCounts") {
    paramsList <- coerceObj(params, to = "list")
    inputMAT <- coerceObj(input, to = "matrix", keepNames = FALSE)
  } else {
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
  } else {
    stop("An error occurred!")
  }

  # Make sure that we have at least 2 samples
  if (length(mutSignatures::getSampleIdentifiers(input)) < 2) {

    stop("mutSignatures does not allow to extract mutational signatures from a single sample.
         Please, add more samples and run again!")
  }

  # Make sure that we are extracting at least 2 signatures
  paramsList$num_processesToExtract <- as.integer(paramsList$num_processesToExtract)
  if (paramsList$num_processesToExtract < 1) {

    stop("mutSignatures does not allow to extract mutational signatures from a single sample.
         Please, add more samples and run again!")
  } else if (paramsList$num_processesToExtract == 1 &&
             paramsList$analyticApproach == "denovo") {

    zz <- apply(inputMAT, 1, sum, na.rm = TRUE)
    tot.zz <- sum(zz, na.rm = TRUE)

    zz2 <- data.frame(Sign.1 = (zz / tot.zz))
    zz2 <- data.frame(Sign.1 = (zz / tot.zz))
    zzr <- getMutationTypes(input)
    rownames(zz2) <- zzr

    out.sig <- as.mutation.signatures(zz2)

    message("Since the user requested to extract a single mutational signature,
            bootstrapping was not performed. A normalized average signature
            is returned instead!")

    return(out.sig)
  }

  if (paramsList$analyticApproach == "denovo") {
    deconvData <- deconvoluteMutCounts(input_mutCounts = inputMAT,
                                       params = paramsList)
  } else {
    stop("An error occurred!")
  }

  mutProcesses <- list()
  final.proc <- data.frame(deconvData$processes, stringsAsFactors = FALSE)
  colnames(final.proc) <- paste("Sign.", sapply(1:ncol(final.proc),
                                                (function(n) {
                                                  leadZeros(n, (10 * ncol(final.proc)))
                                                })), sep = "")
  rownames(final.proc) <- getMutationTypes(input)
  signResult <- mutSignatures::as.mutation.signatures(final.proc)
  final.expo <- data.frame(deconvData$exposure, stringsAsFactors = FALSE)
  if (getFwkParam(params, "approach") != "counts") {
    final.expo <- data.frame(sapply(1:ncol(final.expo), function(cjj) {
      inputColsums[cjj] * final.expo[, cjj]/sum(final.expo[, cjj])
    }), row.names = NULL, stringsAsFactors = FALSE)
  }
  tryCatch({
    colnames(final.expo) <- getSampleIdentifiers(input)[deconvData$includedSampleId]
  }, error = function(e) {
    NULL
  })
  rownames(final.expo) <- colnames(final.proc)
  expoResult <- mutSignatures::as.mutsign.exposures(x = final.expo, samplesAsCols = TRUE) #added TRUE
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


#' Deconvolute Mutation Counts.
#'
#' Characterize mutational signatures from cancer-derived somatic mutational catalogs.
#'
#' @param input_mutCounts numeric matrix of Mutation Type Counts
#' @param params list including all parameters required for running the analysis
#'
#' @return list including all the results from the deconvolution analysis.
#' This function is called within thedecipherMutationalProcesses() function
#' after parameters and input data have been validated
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("deconvoluteMutCounts")
#' y <- mutSignatures:::deconvoluteMutCounts(input_mutCounts = as.matrix(x$muts),
#'                                           params = as.list(x$params))
#' y$processes[1:10,]
#'
#'
#' @import foreach
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @keywords internal
deconvoluteMutCounts <- function(input_mutCounts, params)
{
  j <- NULL
  num_totIterations <- params$num_totIterations
  num_processesToExtract <- params$num_processesToExtract
  distanceFunction <- params$distanceFunction
  thresh_removeWeakMutTypes <- params$thresh_removeWeakMutTypes
  num_parallelCores <- params$num_parallelCores
  guided <- params$guided
  debugStatus <- params$debug
  num_totReplicates <- params$num_totReplicates
  thresh_removeLastPercent <- params$thresh_removeLastPercent
  colnames(input_mutCounts) <- NULL
  rownames(input_mutCounts) <- NULL
  input_mutCounts <- base::as.matrix(input_mutCounts)
  bckgrnd.removed.mutCounts <- removeWeak(input_mutCounts, params)

  # Make sure we have enough counts
  # else, remove samples
  # Carry sample identifiers forard as well
  includedSampleId <- 1:ncol(input_mutCounts)
  if (sum(apply(bckgrnd.removed.mutCounts$output.mutCounts, 2, sum) < 1) > 0) {

    # Make a working copy of the matrix,
    # do this iteratively, just cause
    tmp_input_counts <- input_mutCounts
    tmp_includedSampleId <- includedSampleId
    RMV <- which(apply(bckgrnd.removed.mutCounts$output.mutCounts, 2, sum) < 1)
    while(length(RMV) > 0) {
      tmp_input_counts <- tmp_input_counts[, -RMV]
      tmp_includedSampleId <- tmp_includedSampleId[-RMV]
      bckgrnd.removed.mutCounts <- removeWeak(tmp_input_counts, params)
      RMV <- which(apply(bckgrnd.removed.mutCounts$output.mutCounts, 2, sum) < 1)
    }

    # Inform about removed samples
    message(paste("Some of the samples (n=",
                  (ncol(input_mutCounts) - ncol(tmp_input_counts)),
                  ") were ineligible for analysis and were removed.", sep = ""))
    input_mutCounts <- tmp_input_counts
    includedSampleId <- tmp_includedSampleId
  }

  # Back to analysis
  bckgrnd.removed.mutset <- bckgrnd.removed.mutCounts$removed.mutset
  bckgrnd.removed.mutCounts <- bckgrnd.removed.mutCounts$output.mutCounts
  total.mutationTypes <- nrow(bckgrnd.removed.mutCounts)
  total.samples <- ncol(bckgrnd.removed.mutCounts)
  if (guided) {
    if (!debugStatus) {
      guide.W <- suppressMessages(
        extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                          params = params, bootStrap = FALSE))
    } else {
      guide.W <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                   params = params, bootStrap = FALSE)
    }
    guide.W <- guide.W$Wk
  } else {
    guide.W <- 0
  }
  if (num_parallelCores < 2) {
    muCounts.checkDF <- tryCatch(
      lapply(1:num_totIterations, (function(j) {
        if (debugStatus) {
          if (j %in% base::as.integer(seq(1, num_totIterations, length.out = 100))) {
            message(paste("(", j, ")", sep = ""), appendLF = FALSE)
          }
        }
        if (!debugStatus) {
          tmp.out <- suppressMessages(
            extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                              bootStrap = TRUE, params = params))

        } else {

          tmp.out <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                       bootStrap = TRUE, params = params)
        }

        if (guided) {
          re.ORD <- rep(0, num_processesToExtract)

          for (ki in 1:num_processesToExtract) {

            my.i <- order(apply(abs(tmp.out$Wk - guide.W[, ki]), 2, sum))

            if (ki > 1) {
              my.i[re.ORD[1:(ki - 1)]] <- max(my.i) + 1
            }
            re.ORD[ki] <- which.min(my.i)
          }

        } else {
          re.ORD <- 1:num_processesToExtract
        }
        if (num_processesToExtract > 1) {
          tmp.out$Wk <- tmp.out$Wk[, re.ORD]
          tmp.out$Hk <- tmp.out$Hk[re.ORD, ]
        } else {
          tmp.out$Wk <- tmp.out$Wk
          tmp.out$Hk <- rbind(tmp.out$Hk)
        }
        tmp.out
      })), error = (function(e) {
        print(e)
      }))
    if (debugStatus) {
      message("Done!", appendLF = TRUE)
    }
  } else {
    max.cores <- parallel::detectCores()
    max.cores <- max.cores - 1
    max.cores <- ifelse(max.cores < 1, 1, max.cores)
    use.cores <- ifelse(1 <= num_parallelCores & num_parallelCores <= max.cores,
                        num_parallelCores, max.cores)
    if (debugStatus) {
      cl <- suppressMessages(
        parallel::makeCluster(use.cores, outfile = ""))
    } else {
      cl <- suppressMessages(parallel::makeCluster(use.cores))
    }
    print(paste("Extracting", num_processesToExtract, "mutational signatures X",
                num_totIterations, "iterations using", use.cores,
                "cores"))
    suppressMessages(doParallel::registerDoParallel(cl))

    # Prep before exporting
    alexaNMF <- alexaNMF
    leadZeros <- leadZeros
    extractSignatures <- extractSignatures
    frequencize <- frequencize
    bootstrapCancerGenomes <- bootstrapCancerGenomes
    chihJenNMF <- chihJenNMF

    stuffToExp <- c()
    #stuffToExp <- c("alexaNMF", "leadZeros", "extractSignatures",
    #                "frequencize", "bootstrapCancerGenomes", "chihJenNMF")
    suppressMessages(parallel::clusterExport(cl, stuffToExp))
    muCounts.checkDF <- tryCatch(
      foreach::foreach(j = (1:num_totIterations),
                       .verbose = TRUE,
                       .packages = c("stats", "mutSignatures")) %dopar% {

                         if (j %in% base::as.integer(seq(1, num_totIterations, length.out = 100))) {
                           message(paste("(", j, ")", sep = ""), appendLF = FALSE)
                         }
                         tmp.out <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                                      params = params)
                         if (guided) {
                           re.ORD <- rep(0, num_processesToExtract)
                           for (ki in 1:num_processesToExtract) {
                             my.i <- order(apply(abs(tmp.out$Wk - guide.W[, ki]), 2, sum))
                             if (ki > 1) {
                               my.i[re.ORD[1:(ki - 1)]] <- max(my.i) + 1
                             }
                             re.ORD[ki] <- which.min(my.i)
                           }
                         } else {
                           re.ORD <- 1:num_processesToExtract
                         }
                         if (num_processesToExtract > 1) {
                           tmp.out$Wk <- tmp.out$Wk[, re.ORD]
                           tmp.out$Hk <- tmp.out$Hk[re.ORD, ]
                         } else {
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
                                             hall = H.all,
                                             cnt_errors = errors.all,
                                             cnt_reconstructed = reconstruct.all,
                                             params)
  stability.check <- evaluateStability(wall = fltr.mutCounts.data$Wall,
                                       hall = fltr.mutCounts.data$Hall, params)
  final.mutCounts.data <- addWeak(mutationTypesToAddSet = bckgrnd.removed.mutset,
                                  processes_I = stability.check$centroids,
                                  processesStd_I = stability.check$centroidStd,
                                  Wall_I = fltr.mutCounts.data$Wall,
                                  genomeErrors_I = fltr.mutCounts.data$mutCounts.errors,
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
  deconvoluted.results$includedSampleId <- includedSampleId
  return(deconvoluted.results)
}

#' Extract Signatures from Genomic Mutational Catalogs.
#'
#' Extract mutational signatures after the input Data and the
#' input parameters have been checked andvalidated.
#'
#' @param mutCountMatrix numeric matrix of mutation counts
#' @param params list including all parameters for performing the analysis
#' @param bootStrap logical, shall bootstrapping be performed
#'
#' @return list including the following elements
#' \enumerate{
#'    \item \bold{Wall}:  all extracted signatures
#'    \item \bold{Hall}: all extracted exposures
#'    \item \bold{mutCounts.reconstructed}: fitted values
#'    \item \bold{mutCounts.errors}: residuals
#'  }
#'
#' @details This is one of the core functions included in the original mutSignatures R library,
#' and in the WTSI MATLAB framework. This is an internal function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("extractSignatures")
#' y <- mutSignatures:::extractSignatures(mutCountMatrix = as.matrix(x$muts),
#'                                        params = as.list(x$params), bootStrap = TRUE)
#' y$Wk[1:10,]
#'
#'
#' @importFrom proxy dist
#'
#' @keywords internal
extractSignatures <- function (mutCountMatrix,
                               params,
                               bootStrap = TRUE)
{
  num_processesToExtract <- params$num_processesToExtract
  approach <- params$approach
  algorithm <- params$algorithm
  eps <- params$eps

  if (bootStrap) {
    bstrpd.result <- bootstrapCancerGenomes(mutCountMatrix)
  } else {
    bstrpd.result <- mutCountMatrix
  }

  #if (approach != "counts") {
  #  frq.bstrpd <- frequencize(bstrpd.result)
  #  bstrpd.result <-  frq.bstrpd$freqs
  #}
  bstrpd.result[bstrpd.result < eps] <- eps

  if (algorithm %in% c("brunet", "alexa")) {
    nmf.results <- alexaNMF(v = bstrpd.result,
                            r = num_processesToExtract,
                            params = params)
  } else {
    nmf.results <- chihJenNMF(v =  bstrpd.result,
                              r = num_processesToExtract,
                              params = params)
  }
  # nmf.results
  tmp.w <- nmf.results$w
  tmp.h <- nmf.results$h

  # This step seems useless for the new approach, as it divides or multiplies by 1
  # Keep for consistency and cause we let the user select what approach to use
  for (jj in 1:num_processesToExtract) {
    tmp.tot <- sum(tmp.w[, jj])
    tmp.w[, jj] <- tmp.w[, jj]/tmp.tot
    tmp.h[jj, ] <- tmp.h[jj, ] * tmp.tot
  }

  # modify for frequentized approach
  #if (approach != "counts"){
  #  if (length(frq.bstrpd$colSums) == ncol(tmp.h)) {
  #    tmp.h <- sapply(1:length(frq.bstrpd$colSums), function(ai) {
  #      frq.bstrpd$colSums[ai] * tmp.h[,ai] / sum(tmp.h[,ai], na.rm = TRUE)
  #    })
  #  }
  #}

  mutCountMatrix.reconstructed <- tmp.w %*% tmp.h

  result.list <- list()
  result.list$Wk <- tmp.w
  result.list$Hk <- tmp.h
  result.list$mutCounts.reconstructed <- mutCountMatrix.reconstructed
  result.list$mutCounts.errors <- bstrpd.result - mutCountMatrix.reconstructed

  if (params$logIterations != "lite") {
    result.list$inputMatrix <- bstrpd.result
    result.list$cosDist <- proxy::dist(rbind(base::as.vector(bstrpd.result),
                                             base::as.vector(mutCountMatrix.reconstructed)),
                                       method = "cosine")[1]
  }

  return(result.list)
}


#' Set Parameters for Extracting Mutational Signatures.
#'
#' Create an object including all parameters required for running the mutSignatures framework.
#'
#' @param num_processesToExtract integer, number of mutational signatures to extract
#' @param num_totIterations integer, total number of iterations (bootstrapping)
#' @param num_parallelCores integer, number of cores to use for the analysis
#' @param thresh_removeWeakMutTypes numeric, threshold for filtering out under-represented mutation types
#' @param thresh_removeLastPercent numeric, threshold for removing outlier iteration results
#' @param distanceFunction string, method for calculating distances. Default method is "cosine"
#' @param num_totReplicates integer, number of replicates while checking stability
#' @param eps numeric, close-to-zero positive numeric value for replacing zeros and
#' preventing negative values to appear in the matrix during NMF
#' @param stopconv integer, max number of stable iterations before termination. Defaults to 20000.
#' @param niter integer, max number of iterations to run. Defaults to 1000000
#' @param guided logical, shall clustering be guided to improve aggregation upon bootstrapping
#' @param debug logical, shall the analysis be run in DEBUG mode
#' @param approach string, indicating whether to model absolute counts ("counts") or
#' per_mille frequency ("freq"). Defaults to "freq".
#' @param stopRule = string, use the sub-optimal termination rule ("AL") from the WTSI package
#' (actually, iterations won't terminate, so niter will most certainly reached) or our
#' efficient termination rule ("DF"). Defaults to "DF". The "AL" option is implemented for
#' compatibility reasons, but not recommended.
#' @param algorithm string, algorithm to be used for NMF. Set to "brunet", or "alexa" for using the standard algorithm (Brunet's),
#' otherwise the alternative "chihjen" algorithm will be used.
#' @param logIterations string indicating if storing and returining all intermediates,
#' or only final results. Defaults to "lite", i.e. returns full output and limited intermediates.
#' Alternatively, set to "full".
#' @param seed integer, seed to set for reproducibility
#'
#' @return Object including all parameters for running the analysis
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'   \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
#'  }
#'
#' @examples
#' library(mutSignatures)
#' # defaults params
#' A <- setMutClusterParams()
#' A
#' # A second example, set num_processes
#' B <- setMutClusterParams(num_processesToExtract = 5)
#' B
#'
#'
#' @export
setMutClusterParams <- function(num_processesToExtract = 2,
                                num_totIterations = 10,
                                num_parallelCores = 1,
                                thresh_removeWeakMutTypes = 0.01,
                                thresh_removeLastPercent = 0.07,
                                distanceFunction = "cosine",
                                num_totReplicates = 100,
                                eps = 2.2204e-16,
                                stopconv = 20000,
                                niter = 1000000,
                                guided = TRUE,
                                debug = FALSE,
                                approach = "freq",
                                stopRule = "DF",
                                algorithm = "brunet",
                                logIterations = "lite",
                                seed = 12345)
{
  #
  # Step-by-step Parameter Validation and preparation
  paramList <- list()
  #
  if (!((is.numeric(num_processesToExtract[1]) & num_processesToExtract[1] > 0) ))
    stop("Provide a reasonable number of signatures/processes to extract")
  paramList$num_processesToExtract <- round(num_processesToExtract[1])
  #
  if (!(is.numeric(num_totIterations[1]) & num_totIterations[1] > 0))
    stop("Provide a reasonable number of iterations to run (Bootstrapping)")
  paramList$num_totIterations <- round(num_totIterations[1])
  #
  if (!(is.numeric(num_parallelCores[1]) & num_parallelCores[1] > 0))
    stop("Provide a reasonable number of CPU cores to use for the analysis")
  paramList$num_parallelCores <- round(num_parallelCores[1])
  #
  #paramList$perCore.iterations <- as.integer(paramList$num_totIterations / paramList$num_parallelCores)
  #paramList$perCore.iterations <- ifelse(paramList$perCore.iterations < 1, 1, paramList$perCore.iterations)
  #
  if (!(is.numeric(thresh_removeWeakMutTypes[1]) & thresh_removeWeakMutTypes[1] >= 0 & thresh_removeWeakMutTypes[1] < 1))
    stop("Provide a reasonable (0.00-0.99) number of low-occurring mutation types to remove from the input before starting the analysis")
  paramList$thresh_removeWeakMutTypes <- thresh_removeWeakMutTypes[1]
  #
  if (!(is.numeric(thresh_removeLastPercent[1]) & thresh_removeLastPercent[1] >= 0 & thresh_removeLastPercent[1] < 1))
    stop("Provide a reasonable (0.00-0.99) number for filtering out poor iterations")
  paramList$thresh_removeLastPercent <- thresh_removeLastPercent[1]
  #
  allowed.dist.methods <- c("Braun-Blanquet", "Chi-squared", "correlation", "cosine", "Cramer", "Dice",
                            "eDice", "eJaccard", "Fager", "Faith", "Gower", "Hamman", "Jaccard",
                            "Kulczynski1", "Kulczynski2", "Michael", "Mountford", "Mozley", "Ochiai",
                            "Pearson", "Phi", "Phi-squared", "Russel", "simple matching", "Simpson",
                            "Stiles", "Tanimoto", "Tschuprow", "Yule", "Yule2", "Bhjattacharyya",
                            "Bray", "Canberra", "Chord", "divergence", "Euclidean", "fJaccard",
                            "Geodesic", "Hellinger", "Kullback", "Levenshtein", "Mahalanobis",
                            "Manhattan", "Minkowski", "Podani", "Soergel","supremum", "Wave", "Whittaker")
  if (!(is.character(distanceFunction[1]) & distanceFunction[1] %in% allowed.dist.methods))
    stop("Unknown method for calculating distances. For options, run: <<summary(proxy::pr_DB)>>")
  paramList$distanceFunction <- distanceFunction[1]
  #
  if (!(is.numeric(num_totReplicates[1]) & num_totReplicates[1] > 99))
    stop("Provide a reasonable number of replicates for stability evaluation of the results")
  paramList$num_totReplicates <- round(num_totReplicates[1])
  #
  if (!(is.numeric(eps[1]) & eps[1] > 0 & eps[1] < 0.0001))
    stop("Provide a reasonably small number (0 < n < 0.0001) for data overflow prevention")
  paramList$eps <- eps[1]
  #
  if (!(is.numeric(stopconv[1]) & stopconv[1] >= 500))
    stop("Provide a reasonable large number: number of 'conn-matrix-stable' iterations before stopping NMF")
  paramList$stopconv <- round(stopconv[1])
  #
  if (!(is.numeric(niter[1]) & niter[1] >= 20000))
    stop("Provide a reasonable large number: total NMF iterations")
  paramList$niter <- round(niter[1])
  #
  paramList$guided <- ifelse(guided, TRUE, FALSE)
  paramList$debug <- ifelse(debug, TRUE, FALSE)
  paramList$approach <- ifelse (approach == "counts", "counts", "freq")
  paramList$stopRule <- ifelse (stopRule == "LA", "LA", "DF")
  paramList$algorithm <- ifelse (tolower(algorithm) %in% c("brunet", "alexa"), "brunet", "chihjen")
  paramList$logIterations <- ifelse(tolower(logIterations) %in% c("lite", "light", "li"), "lite", "full")
  #
  return(new(Class = "mutFrameworkParams", params = paramList))
}

## ~~~~~~~~~~~~~~~~~~~~~

#' Custom Fast Combinatorial Nonnegative Least-Square.
#'
#' This function contributes to solve a least square linear
#' problem using the fast combinatorial strategy from Van
#' Benthem et al. (2004). This implementation is similar to that included in the
#' NMF R package by Renaud Gaujoux and Cathal Seoighe,
#' and it is tailored to the data used in the mutational signature analysis.
#' For more info, see: \url{https://CRAN.R-project.org/package=NMF}
#'
#' @param mutCounts numeric matrix including mutation counts
#' @param signatures numeric matrix including mutation signatures
#'
#' @return list, including: (K) a numeric matrix of estimated exposures; and (Pset) a Pset numeric matrix
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs(testN = "custom_fcnnls")
#' y <- mutSignatures:::custom_fcnnls(mutCounts = x$muts, signatures = x$signs)
#' y$coef
#'
#' @keywords internal
custom_fcnnls <- function (mutCounts, signatures)
{

  # Hardcode eps
  eps <-  0

  # Initial checks
  if (sum(c(dim(signatures), dim(mutCounts)) < 1) > 0) {
    stop("Bad input!")
  }
  if (nrow(signatures) != nrow(mutCounts)) {
    stop("Bad input: matrices have imcompatible size")
  }

  # Define W, and comute matrix cross-prods
  W <- matrix(0, ncol(signatures), ncol(mutCounts))
  maxiter = 5 * ncol(signatures)
  xxp <- base::crossprod(signatures)
  xyp <- base::crossprod(signatures, mutCounts)

  # Solve, aka compute K
  K <- base::solve(xxp, xyp)

  # Prep for looping
  Ps <- K > 0
  K[!Ps] <- 0
  D <- K
  Fset <- which(apply(Ps, 2, sum) != ncol(signatures))
  oitr <- 0
  iter <- 0

  while (length(Fset) > 0) {
    oitr <- oitr + 1
    K[, Fset] <- custom_cssls(xxp,
                              xyp[, Fset, drop = FALSE],
                              Ps[, Fset, drop = FALSE])

    keep <- apply(K[, Fset, drop = FALSE] < eps, 2, sum) > 0
    Hset <- Fset[keep]
    if (length(Hset) > 0) {
      nHset <- length(Hset)
      alpha <- matrix(0, ncol(signatures), nHset)
      while (nHset > 0 && (iter < maxiter)) {
        iter <- iter + 1
        alpha[, 1:nHset] <- Inf
        ij <- which(Ps[, Hset, drop = FALSE] & (K[, Hset, drop = FALSE] < eps), arr.ind = TRUE)
        i <- ij[, 1]
        j <- ij[, 2]
        if (length(i) == 0)
          break
        hIdx <- (j - 1) * ncol(signatures) + i
        negIdx <- (Hset[j] - 1) * ncol(signatures) + i
        alpha[hIdx] <- D[negIdx]/(D[negIdx] - K[negIdx])
        alpha.inf <- alpha[, 1:nHset, drop = FALSE]
        minIdx <- max.col(-t(alpha.inf))
        alphaMin <- alpha.inf[minIdx + (0:(nHset - 1) * ncol(signatures))]
        alpha[, 1:nHset] <- matrix(alphaMin, ncol(signatures), nHset, byrow = TRUE)
        D[, Hset] <- D[, Hset, drop = FALSE] - alpha[, 1:nHset, drop = FALSE] *
          (D[, Hset, drop = FALSE] - K[, Hset, drop = FALSE])
        idx2zero <- (Hset - 1) * ncol(signatures) + minIdx
        D[idx2zero] <- 0
        Ps[idx2zero] <- FALSE

        K[, Hset] <- custom_cssls(xxp,
                                  xyp[, Hset, drop = FALSE],
                                  Ps[, Hset, drop = FALSE])

        Hset <- which(apply(K < eps, 2, sum) > 0)
        nHset <- length(Hset)
      }
    }
    W[, Fset] <- xyp[, Fset, drop = FALSE] - xxp %*% K[, Fset, drop = FALSE]

    tmp_coeff <- ifelse(!(Ps[, Fset, drop = FALSE]), 1, 0)
    Jset <- which(apply((tmp_coeff *  W[, Fset, drop = FALSE]) > eps, 2, sum) == 0)

    Fset <- base::setdiff(Fset, Fset[Jset])
    if (length(Fset) > 0) {
      tmp_coeff <- ifelse(!Ps[, Fset, drop = FALSE], 1, 0)
      mxidx <- max.col(t(tmp_coeff * W[, Fset, drop = FALSE]))
      Ps[(Fset - 1) * ncol(signatures) + mxidx] <- TRUE
      D[, Fset] <- K[, Fset, drop = FALSE]
    }
  }
  return(list(coef = K, Pset = Ps))
}

#' Custom CSSLS.
#'
#' This function contributes to solving a nonnegative least square linear
#' problem using normal equations and the fast combinatorial strategy from Van
#' Benthem et al. (2004). This implementation is similar to that included in the
#' NMF R package by Renaud Gaujoux and Cathal Seoighe,
#' and it is tailored to the data used in the mutational signature analysis.
#' For more info, see: \url{https://CRAN.R-project.org/package=NMF}
#'
#' @param CtC numeric matrix
#' @param CtA numeric matrix
#' @param Pset nueric matrix
#'
#' @return a numeric matrix
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs(testN = "custom_cssls")
#' y <- mutSignatures:::custom_cssls(CtC = x$CtC, CtA = x$CtA, Pset = x$Pset)
#' y
#'
#' @keywords internal
custom_cssls <- function (CtC, CtA, Pset)
{
  K = matrix(0, nrow(CtA), ncol(CtA))

  lVar <- nrow(Pset)
  pRHS <- ncol(Pset)
  codedPset <- as.numeric(2^(seq(lVar - 1, 0, -1)) %*% Pset)
  sortedPset <- sort(codedPset)
  sortedEset <- order(codedPset)
  breaks <- diff(sortedPset)
  breakIdx <- c(0, which(breaks > 0), pRHS)
  for (k in seq(1, length(breakIdx) - 1)) {
    cols2solve <- sortedEset[seq(breakIdx[k] + 1, breakIdx[k + 1])]
    vars <- Pset[, sortedEset[breakIdx[k] + 1]]
    K[vars, cols2solve] <-
      solve(CtC[vars, vars, drop = FALSE]) %*%  CtA[vars, cols2solve, drop = FALSE]
  }
  K
}


## ~~~~~~~~~~~~~~~~~~~~~

###
##### Exported functions - user interface
###


#' Attach Mutation Types.
#'
#' Modify a data.frame carrying information about DNA mutation, and add a new column that
#' stores formatted multi-nucleotide types.
#'
#' @param mutData data.frame including information about DNA mutations
#' @param ref_colName string, pointing to the column with information about the sequence of the "reference_allele"
#' @param var_colName string, pointing to the column with information about the sequence of the "variant_allele"
#' @param var2_colName string (optional), pointing to the column with information about the
#' sequence of a second "variant_allele". Can be NULL
#' @param context_colName string, pointing to the column with information about the nucleotidic "context"
#' @param format integer, indicates the desired mutation type format: (1) N[R>V]N; (2) NN.R>V; (3) R.V[NRN][NVN]
#' @param mutType_dict string, indicates the dictionary to be used for simplifying reverse-complement identical mutation types.
#' It is recommended to use the standard dictionary from COSMIC, by selecting the default value, i.e. "alexa".
#' @param mutType_colName string, column name of the new column added to the data.frame where mutTypes are stored.
#'
#' @return a data.frame including a new column with mutation Types.
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' A <- data.frame(REF = c("A", "T", "G"),
#'                 VAR = c("G", "C", "C"),
#'                 CTX = c("TAG", "GTG", "CGA"),
#'                 stringsAsFactors = FALSE)
#' mutSignatures::attachMutType(mutData = A, ref_colName = "REF",
#'                              var_colName = "VAR", context_colName = "CTX")
#'
#' @export
attachMutType <- function(mutData,
                          ref_colName = "reference_allele",
                          var_colName = "variant_allele",
                          var2_colName = NULL,
                          context_colName = "context",
                          format = 1,
                          mutType_dict = "alexa",
                          mutType_colName = "mutType")
{
  # Validate input data and fields in mutData
  if(!((is.data.frame(mutData) | is.matrix(mutData) ) &
       sum(c(ref_colName, var_colName, context_colName) %in% colnames(mutData)) == 3 ))
    stop ("Issue with the input dataset. Make sure to feed in a data.frame or
          a matrix and double check the name of the fields pointing to chromosome
          name, start and end positions")
  if (!(format %in% c(1,2)))
    stop ("Please, specify a valid format number (example: 1)")
  if (! (is.null(var2_colName))) {
    if (!(var2_colName %in% colnames(mutData)))
      stop ("Invalid var2 column")
  }

  if (!is.character(mutType_colName) |
      length(mutType_colName) > 1)
    stop("Bad mutType_colName")
  if(mutType_colName %in% colnames(mutData))
    stop ("mutType_colName already exists as column name in the current dataset")

  # convert factors to chars
  mutData <- data.frame(mutData, stringsAsFactors = FALSE, row.names = NULL)
  my.key.cols <- c(ref_colName, var_colName, var2_colName, context_colName)
  my.key.cols <- my.key.cols[!is.na(my.key.cols)]
  for (clmn in my.key.cols) {
    mutData[,clmn] <- base::as.character(base::as.vector(mutData[,clmn]))
  }

  message("Assigning mutation types ", appendLF = FALSE)

  mutData[,mutType_colName] <- sapply(1:nrow(mutData), (function(i){

    if (nrow(mutData) > 1000 & i %in% base::as.integer(seq(1, nrow(mutData),length.out = 20)))
      message(".", appendLF = FALSE)

    # first, extract elems and check for middle base to match the reference
    ctx.len <- nchar(mutData[i,context_colName])
    half.ln <- (ctx.len - 1) / 2
    mid.seq <- substr(mutData[i,context_colName], (half.ln + 1), (half.ln + 1))
    pre.seq <- substr(mutData[i,context_colName], 1, half.ln)
    post.seq <- substr(mutData[i,context_colName], (half.ln + 2), ctx.len)

    if((mid.seq !=  mutData[i,ref_colName]) ||
       (is.null(var2_colName) & mid.seq == mutData[i,var_colName]) ||
       (tryCatch({
         ChK0 <- (mid.seq == mutData[i,var_colName] & mid.seq == mutData[i,var2_colName])
         ifelse(length(ChK0) == 1, ChK0, FALSE)
       }, error = function(e) {FALSE}))) {
      # no match means --> NA
      mut.base <- NA
    } else {
      #mut.variant to use in case var2 is specified
      if (mutData[i,ref_colName] != mutData[i,var_colName]) {
        mut.base <- mutData[i,var_colName]
      } else if (!is.null(var2_colName) ){
        if ( mutData[i,ref_colName] != mutData[i,var2_colName]) {
          mut.base <- mutData[i,var2_colName]
        } else {
          mut.base <- NA
        }
      } else {
        mut.base <- NA
      }

      # match, format and return according to a standard format (for now)
      if (is.na(mut.base)) {
        NA
      } else {
        paste(mid.seq, ".", mut.base, "[", pre.seq, mid.seq, post.seq, "][", pre.seq, mut.base, post.seq, "]", sep = "", collapse = "")
      }
    }
  }))
  message(". Done!", appendLF = TRUE)

  if (sum(is.na(mutData[,mutType_colName])) > 0) {
    message(paste("Removing",sum(is.na(mutData[,mutType_colName])), "positions."))
    mutData <- mutData[!is.na(mutData[,mutType_colName]),]
  }

  # Now, apply revCompl transformation
  message("Now applying RevCompl transformation", appendLF = FALSE)
  if (mutType_dict == "alexa") {
    idx <- grep("^((G|A)\\.)", mutData[,mutType_colName])
    mutData[idx,mutType_colName] <- sapply(mutData[idx,mutType_colName], (function(seq){

      base.wt  <- revCompl(gsub("\\..+$", "", seq))
      base.mut <- revCompl(gsub("^.+\\.", "", gsub("\\[.+$", "", seq)))
      seq.wt   <- revCompl(gsub("^.+\\[", "", gsub("\\]\\[.+$", "", seq)))
      seq.mut  <- revCompl(gsub("^.+\\]\\[", "", gsub("\\]$", "", seq)))
      paste(base.wt,".",base.mut, "[", seq.wt, "][", seq.mut, "]", sep = "", collapse = "")
    }))

  } else if (mutType_dict == "custom") {
    idx <- grep("^((G|T)\\.)", mutData[,mutType_colName])
    mutData[idx,mutType_colName] <- sapply(mutData[idx,mutType_colName], (function(seq){
      base.wt  <- revCompl(gsub("\\..+$", "", seq))
      base.mut <- revCompl(gsub("^.+\\.", "", gsub("\\[.+$", "", seq)))
      seq.wt   <- revCompl(gsub("^.+\\[", "", gsub("\\]\\[.+$", "", seq)))
      seq.mut  <- revCompl(gsub("^.+\\]\\[", "", gsub("\\]$", "", seq)))
      paste(base.wt,".",base.mut, "[", seq.wt, "][", seq.mut, "]", sep = "", collapse = "")
    }))
  }

  message(". Done!", appendLF = TRUE)
  message("Final formatting", appendLF = FALSE)
  #Attach the format of interest
  mutData[,mutType_colName] <- sapply(mutData[,mutType_colName], (function(seq){
    base.wt  <- gsub("\\..+$", "", seq)
    base.mut <- gsub("^.+\\.", "", gsub("\\[.+$", "", seq))
    seq.wt   <- gsub("^.+\\[", "", gsub("\\]\\[.+$", "", seq))
    seq.mut  <- gsub("^.+\\]\\[", "", gsub("\\]$", "", seq))
    half.len <- (nchar(seq.wt) - 1 ) / 2
    pre.seq <- substr(seq.wt, 1, half.len)
    post.seq <-substr(seq.wt, half.len + 2, nchar(seq))

    if (format == 1) {
      # --> N[N>M]N
      paste(pre.seq, "[", base.wt, ">", base.mut, "]", post.seq, sep = "", collapse = "")
    } else if (format == 2) {
      # --> NN.N>M
      paste(pre.seq, post.seq, ".", base.wt, ">", base.mut, sep = "", collapse = "")
    }  else {
      # --> N.M[NNN][NMN]
      paste(base.wt,".",base.mut, "[", seq.wt, "][", seq.mut, "]", sep = "", collapse = "")
    }
  }))

  message(". Done!", appendLF = TRUE)

  return(mutData)
}


#' Extract Variants from XvarlinkData.
#'
#' Extract Variants from data stored as XvarlinkData.
#'
#' @param xvarLink_data character vector, including mutation data embedded in XvarlinkData
#'
#' @return a data.frame including mutations as well as corresponding reference nucleotides.
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("extractXvarlinkData")
#' y <- mutSignatures:::extractXvarlinkData(xvarLink_data = x)
#' y
#'
#'
#' @export
extractXvarlinkData <- function(xvarLink_data) {
  tmpVars <- sub("^.*&var=[[:alnum:]]+(,.*)&.*$", "\\1", xvarLink_data)
  tmpVars[!grepl("^,.+", tmpVars)] <- NA
  tmpVars <- strsplit(tmpVars, ",")
  tmpVars <- do.call(rbind, lapply(1:length(tmpVars), function(i){
    if (tmpVars[[i]][1] == "" & length(tmpVars[[i]]) == 5) {
      tmpVars[[i]][2:5]
    } else {
      c(NA, NA, NA, NA)
    }
  }))
  out <- data.frame(chrXvar = base::as.character(tmpVars[,1]),
                    posXvar = base::as.numeric(tmpVars[,2]),
                    refXvar = base::as.character(tmpVars[,3]),
                    mutXvar = base::as.character(tmpVars[,4]),
                    stringsAsFactors = FALSE)
  return(out)
}

#' Filter Single Nucleotide Variants.
#'
#' Remove entries corresponding to non-SNV, such as insertions and deletions.
#'
#' @param dataSet data.frame including variant information
#' @param seq_colNames character vector with the names of the columns storing variant data
#'
#' @return a filtered data.frame only including SNVs
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("filterSNV")
#' nrow(x)
#' y <- mutSignatures::filterSNV(dataSet = x,
#'                               seq_colNames = c("REF", "ALT"))
#' nrow(y)
#'
#'
#' @export
filterSNV <- function(dataSet,
                      seq_colNames)
{
  if (!(is.data.frame(dataSet) | is.matrix(dataSet) ) |
      sum(!seq_colNames %in% colnames(dataSet)) > 0 |
      length(seq_colNames) < 2) {
    stop("Bad input or seq_colNames not found")
  }
  check.tab <- sapply(1:length(seq_colNames), (function(i){
    tmp <- gsub("[[:space:]]", "", dataSet[,seq_colNames[i]])
    toupper(tmp) %in% c("A","C","G","T")
  }))
  toKeep <- apply(check.tab, 1, (function(rw){
    sum(rw) == length(rw)
  }))
  out <- dataSet[toKeep,]
  rownames(out) <- NULL
  return(out)
}

#' Convert Mutation COunts to PerMille Frequencies.
#'
#' Convert Mutation COunts to frequencies. Typically, a permille frequence is returned.
#' In other words, the resulting number indicates the expected mutation count if the genome hat a
#' total of 1000 mutations. This way, the MutSignatures analysis will be
#' less biased toward the hyper-mutator samples.
#'
#' @param countMatrix numeric matrix of mutation counts
#' @param permille ligucal, shall the permille conversion be used instead of the standard frequency
#'
#' @return list including colSums (mutation burden of each sample) and freqs (matrix of frequencies)
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' A <- cbind(c(7, 100, 90, 1000), c(1, 3, 5, 9))
#' fA <- mutSignatures::frequencize(A)
#' fA$freqs
#'
#' @export
frequencize <- function(countMatrix,
                        permille = TRUE)
{
  out <- list()
  cf <- ifelse(permille, 1000, 1)

  out[["colSums"]] <- apply(countMatrix, 2, sum)
  out[["freqs"]] <- cf * apply(countMatrix, 2, (function(clmn){clmn/sum(clmn)}))
  return(out)
}

#' Obtain COSMIC mutational Signatures.
#'
#' Obtain latest mutational Signature definitions from COSMIC. FOr more info, please visit: \url{http://cancer.sanger.ac.uk/}
#'
#' @param forceUseMirror logical, shall signatures be downloaded from a mirror. Set to TRUE if the COSMIC
#' server goes down.
#' @param asMutSign logical, shall data be returned as a mutSignatures-class object. Defaults to TRUE
#'
#' @return an object storing COSMIC mutational signature data
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @importFrom utils read.delim read.csv
#' @export
getCosmicSignatures <- function(forceUseMirror = FALSE, asMutSign = TRUE)
{
  mutType.labels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
                      "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
                      "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
                      "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                      "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
                      "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
                      "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
                      "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
                      "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  cosmic.url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"

  #initialize variable
  my_fullW <- NULL

  if (!forceUseMirror) {
    my_fullW <- tryCatch({TMP <- suppressWarnings(read.delim(cosmic.url, header = TRUE));
    rownames(TMP) <- TMP$Somatic.Mutation.Type;
    TMP <- TMP[,grep("Signature", colnames(TMP))];
    TMP},
    error = function(e) NULL)
  }

  # private mirror
  if (is.null(my_fullW)) {
    private.mirror.url <- "http://www.labwizards.com/rlib/mutSignatures/cosmic.signatures.csv"
    my_fullW <- tryCatch({suppressWarnings(read.csv(private.mirror.url,
                                                    header = TRUE, as.is = TRUE, row.names = 1))},
                         error = function(e2) { NULL })
  }
  if (is.null(my_fullW)) {
    message("An error occurred!")
    return(NULL)
  } else {
    if (sum(mutType.labels %in% rownames(my_fullW)) == length(mutType.labels)) {
      
      colnames(my_fullW) <- sub("Signature", "COSMIC", colnames(my_fullW))
      obj2rt <- my_fullW[mutType.labels,]
      if(asMutSign)
        obj2rt <- mutSignatures::as.mutation.signatures(obj2rt)

      return(obj2rt)
    } else {
      message("An error occurred!")
      return(NULL)
    }
  }
}

#' Import Mutation data from VCF files.
#'
#' Import Mutation data from VCF files. The columns are expected in the following order:
#' c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"). Optional columns
#' can be present to inform about sample ID or other info.
#'
#' @param vcfFiles character vector, includes the names of the VCF files to be analyzed
#' @param sampleNames character vector with alternative sample names (otherwise,
#' VCF file names will be ised to identify each sample).
#'
#' @return a concatenated data.frame with all variants found in the input VCF files. Sample
#' ID is stored in the "SAMPLEID" column.
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @importFrom utils read.delim
#' @export
importVCFfiles <- function(vcfFiles, sampleNames = NULL){

  my.colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                   "FILTER", "INFO", "FORMAT", "XTR1", "XTR2", "XTR3")

  if (is.null(sampleNames) | length(sampleNames) != length(vcfFiles)){

    bypassNames <- paste("sample", 1:length(vcfFiles), sep = ".")

    sampleNames <- vcfFiles
    sampleNames <- sub("\\.vcf$", "", sub("^.*(\\\\|/)", "", tolower(sampleNames)))

    sampleNames[sampleNames == ""] <- bypassNames[sampleNames == ""]

  }
  out <- sapply(1:length(vcfFiles), (function(j){
    x <- vcfFiles[j]
    if (!file.exists(x)) {
      NULL
    } else {
      tmpVCF <- read.delim(x, comment.char = '#', header = F, stringsAsFactors = F)

      for (i in 1:ncol(tmpVCF)) {
        colnames(tmpVCF)[i] <- my.colnames[i]
      }
      tmpVCF$SAMPLEID <- sampleNames[j]
      tmpVCF
    }
  }), simplify = FALSE, USE.NAMES = TRUE)
  out <- do.call(rbind, out)
  #names(out) <- sub("\\.vcf$", "", names(out))
  return(out)
}





##
### Plotting
##

#' Plot Signature Exposure Profiles.
#'
#' Build a barplot to visualize exposures to mutation signatures.
#'
#' @param mutCount a data.frame including mutation Counts
#' @param top integer, max number of samples to include in the plot
#'
#' @return a plot (ggplot2 object)
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme theme_minimal element_blank element_line unit element_text scale_y_continuous
#' @export
plotSignExposures <- function(mutCount, top = 50) {
  # avoid NOTEs
  count <- NULL
  feature <- NULL
  Signature <- NULL

  sampleLabs <- rownames(mutCount)
  rownames(mutCount) <- NULL

  nuSmplOrder <- order(apply(mutCount, 1, sum), decreasing = TRUE)
  mutCount <- mutCount[nuSmplOrder,]
  if (!is.null(sampleLabs))
    sampleLabs <- sampleLabs[nuSmplOrder]

  rownames(mutCount) <- NULL

  if (is.null(top)) {
    mutDF <- table2df(dataMatrix = mutCount)
  } else if(is.na(base::as.numeric(top[1]))){
    mutDF <- table2df(dataMatrix = mutCount)
  } else if (top[1] > nrow(mutCount) | top[1] < 2) {
    mutDF <- table2df(dataMatrix = mutCount)
  } else {
    mutDF <- table2df(dataMatrix = mutCount[1:top,])
    if (!is.null(sampleLabs))
      sampleLabs <- sampleLabs[1:top]
  }

  mutDF$sample <- 100000 + base::as.numeric(base::as.character(mutDF$sample))
  mutDF$sample <- base::as.character(mutDF$sample)

  tryCatch({mutDF$inputLabel <- do.call(c, lapply(sampleLabs,
                                                  rep, times = ncol(mutCount)))},
           error = function(e) {
             message("damn")})

  mutDF$Signature <- factor(mutDF$feature, levels = rev(colnames(mutCount)))
  bp <- ggplot2::ggplot(data=mutDF, ggplot2::aes(x=sample, y=count, fill=Signature)) +
    ggplot2::geom_bar(stat="identity")
  bp <- bp + ggplot2::theme_minimal() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(axis.line.y = ggplot2::element_line(colour = "black", size = 0.75),
                   axis.line.x = ggplot2::element_line(colour = "black", size = 0.75),
                   axis.ticks.y = ggplot2::element_line(colour = "black", size = 1),
                   axis.ticks.length = ggplot2::unit(x = 6, "points"),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_y_continuous(expand=c(0,0),
                                limits = c(0, 1.2 * max(apply(mutCount, 1, function(rx) sum(rx, na.rm = TRUE)))))

  return(bp)
}

#' Plot Mutation Signature Profiles.
#'
#' Build a barplot to visualize the relative abundance of mutation counts in a mutational
#' signature or biological sample of interest.
#'
#' @param mutCounts data.frame including mutation types counts or frequencies, such as a
#' data.frame of mutation counts from samples, or mutation type frequencies from a mutational signature.
#' @param mutLabs character vector, labels to be used for the mutation types
#' @param freq logical, shall frequency be plotted rather than counts. Defaults to TRUE
#' @param ylim values used for ylim. Defaults to "auto" (ylim automatically set)
#' @param ylab string, used as y-axis title. Defaults to "Fraction of Variants"
#' @param xlab string, used as x-axis title. Defaults to "Sequence Motifs"
#' @param xaxis_cex numeric, cex value for the xaxis
#' @param cols character vector, indicates the colors to be used for the bars. It typically requires 6 colors.
#' @param main string, tutle of the plot. Defaults to "MutType Profile"
#'
#' @return NULL. A plot is printed to the active device.
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @importFrom graphics barplot box axis mtext
#' @export
plotMutTypeProfile <- function(mutCounts,
                               mutLabs,
                               freq = TRUE,
                               ylim = "auto",
                               ylab = "Fraction of Variants",
                               xlab = "Sequence Motifs",
                               xaxis_cex = 0.475,
                               cols = c("#4eb3d3", "#040404", "#b30000", "#bdbdbd", "#41ab5d", "#dd3497"),
                               main = "MutType Profile")
{

  mutLabs <- toupper(mutLabs)
  if (!sum(regexpr("^(A|C|G|T)(\\[)(A|C|G|T)(>)(A|C|G|T)(\\])(A|C|G|T)$", toupper(mutLabs)) > 0) == length(mutLabs)) {
    message("mutTypes are in a non-standard format (ie, not Sanger 'A[G>A]T' style)")
    altPlot <- 1
  } else {
    altPlot <- 0
  }


  if (length(mutCounts) != length(mutLabs))
    stop("Mutation Type Labels and Mutation Counts do not match!")

  mutCounts[is.na(mutCounts)] <- 0

  if (freq == TRUE) {
    if (ylim[1] == "auto")
      ylim <- c(0,0.2)
    freqs <- mutCounts/sum(mutCounts, na.rm = TRUE)
  } else {
    if (ylim[1] == "auto")
      ylim = c(0,(1.05*max(mutCounts)))
    freqs <- mutCounts
  }

  tinyMut <- gsub("\\](.*)$", "",
                  gsub("^(.*)\\[", "", toupper(mutLabs)))

  # Order the frequencies
  first.out <- lapply(sort(unique(tinyMut)), (function(mt){
    idxKeep <- which(tinyMut == mt)
    tmp.order <- order(mutLabs[idxKeep])
    out <- freqs[idxKeep][tmp.order]
    names(out) <- mutLabs[idxKeep][tmp.order]
    out
  }))
  names.first.out <- sort(unique(tinyMut))

  # Add a spacer
  second.out <- lapply(first.out, (function(vct){
    c(vct, 0)
  }))

  # Wrap together
  third.out <- do.call(c, second.out)
  third.out <- third.out[-length(third.out)]

  # Define color scheme
  col.out <- lapply(1:length(first.out), (function(ii){
    rep(cols[ii], (length(first.out[[ii]]) + 1))
  }))
  col.out <- do.call(c, col.out)
  col.out <- col.out[-length(col.out)]

  third.out[third.out>max(ylim)] <- max(ylim)
  third.shortLab <- gsub("\\[([[:alpha:]]>[[:alpha:]])\\]","-", names(third.out))

  xpos <- graphics::barplot(third.out,
                            col = col.out,
                            ylim = ylim,
                            axes = FALSE,
                            names.arg = FALSE,
                            ylab = ylab,
                            border = NA,
                            main = main)

  xpos.out <- lapply(1:length(first.out), (function(ii){
    c(rep(ii, length(first.out[[ii]])),0)
  }))
  xpos.out <- do.call(c, xpos.out)
  xpos.out <- xpos.out[-length(xpos.out)]

  #(max(ylim)*0.0075*0.2)
  graphics::box()
  graphics::axis(side = 1, tick = FALSE,
                 hadj = 1, cex.axis = xaxis_cex,
                 pos = (max(ylim) * 0.035), family = "mono",
                 at = xpos,
                 labels = third.shortLab,
                 las = 2)

  my.y <- seq(min(ylim), max(ylim), length.out = 4)
  if (max(ylim) > 1.75) {
    my.y.labs <- format(round(my.y, digits = 0))
  } else {
    my.y.labs <- format(round(my.y, digits = 2))
  }

  graphics::axis(side = 2, tick = -0.005,
                 at = my.y, labels = my.y.labs,
                 las = 1, cex.axis = 0.75)

  zz <- unique(xpos.out)
  zz <- zz[zz>0]
  lab.xx <- sapply(zz, (function(i){
    base::mean(xpos[,1][xpos.out == i])
  }))
  graphics::mtext(names.first.out, side = 1, line = 1.5, at = lab.xx, cex = 0.8, col = cols)

  graphics::mtext(xlab, side = 1, line = 3.0, at = base::mean(xpos[,1]), cex = 1)
  # done! No return, just a plot!
}

#' Run a Preliminary Process Assess analysis.
#'
#' This function is an attempt to analyze the relationship between error and k. In other words,
#' the goal of prelimProcessAssess is to visualize the reduction in the error/residuals
#'
#' @param input a mutationCounts-class object
#' @param maxProcess integer, maximum k to test
#' @param approach sting, "counts" or "freq"
#' @param plot logical, shall a plot be printed to the active device
#' @param verbose logical, info about the ongoing analysis be messaged/printed to console
#'
#' @return a data.frame showing the estimated total error with respect to the range of k values
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#'
#' @importFrom stats median
#' @importFrom graphics plot axis lines points title box
#'
#' @export
prelimProcessAssess <- function(input,
                                maxProcess = 6,
                                approach = "counts",
                                plot = TRUE,
                                verbose = TRUE)
{
  tmpParams <- setMutClusterParams(num_processesToExtract = 1,
                                   approach = approach,
                                   stopconv = 10000,
                                   niter = 100000,
                                   thresh_removeWeakMutTypes = 0.01,
                                   thresh_removeLastPercent =0.025,
                                   num_totIterations = 10,
                                   num_parallelCores = 1,
                                   debug = FALSE,
                                   logIterations = "full")

  tmpParams <- coerceObj(x = tmpParams, to = "list")
  input_mutCounts <- coerceObj(x = input, to = "matrix")

  if (ncol(input@counts) < 2) {
    stop("It is not possible to analyze mtational signatures in a single sample.
         Please, include two or more samples before running this analysis.")
  }

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

  # Calc initial (max) error
  medianMaxErr <- sum(base::t(sapply(1:nrow(bckgrnd.removed.mutCounts), (function(i){
    ((stats::median(bckgrnd.removed.mutCounts[i,]) - bckgrnd.removed.mutCounts[i,]))^2
  }))))

  # Message
  if(verbose)
    message("Preliminary Mutational Process Assessment: ", appendLF = FALSE)

  outRes <- lapply(1:maxProcess, (function(i){
    tmpParams$num_processesToExtract <- i
    tmpRes <- suppressWarnings(suppressMessages(
      extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                        params = tmpParams,
                        bootStrap = FALSE)))

    #apply(tmpRes$mutCounts.reconstructed, 2, sum)
    #mutCounts.reconstructed
    #tmpErr <- sum((tmpRes$mutCounts.errors)^2)
    tmpErr <- sum((bckgrnd.removed.mutCounts - tmpRes$mutCounts.reconstructed) ^ 2)

    if(verbose)
      message(".", appendLF = FALSE)
    tmpErr
  }))
  if(verbose)
    message("", appendLF = TRUE)

  err.points <- c(medianMaxErr, do.call(c, outRes))
  err.points <- (-1) * (err.points - max(err.points))
  err.points <- err.points / max(err.points, na.rm = TRUE)

  if (plot) {
    graphics::plot(err.points ~ c(0:(length(err.points) - 1)),
                   type = "n", las = 1, axes = FALSE,
                   ylab = "", xlab = "", main = "Preliminary Mutational Process Assessment")
    graphics::axis(side = 1, at = NULL, cex = 0.65, font = 3)
    graphics::axis(side = 2, at = seq(0, 1, by = 0.2), cex = 0.65, font = 3,
                   labels = format(seq(1, 0, by = -0.2), digits = 2, nsmall = 2), las = 1)
    graphics::lines(c(0:(length(err.points) - 1)), err.points, lwd = 1.5, col = "red2")
    graphics::points(c(0:(length(err.points) - 1)), err.points, pch = 19, col = "gray30")
    graphics::title(xlab="Num of Signatures", line=2.2, cex.lab=1.2, font = 2)
    graphics::title(ylab="Error (% vs. Median)", line=3.1, cex.lab=1.2, font = 2)
    graphics::box()
  }

  return(data.frame(numProcess = c(0:(length(err.points) - 1)),
                    percentErr = (1 - err.points),
                    stringsAsFactors = TRUE))
}


#' Match Mutational Signatures.
#'
#' Analyze the similarity between mutational signatures from different analyses/runs.
#' THis function can be helpful to match de novo extracted signatures with previously
#' described signatures (such as COSMIC), or to reveal signatures that can be
#' identified with alternative NMF algorithms, or that may be due to an algorithm bias.
#'
#' @param mutSign a mutationSignatures object
#' @param reference a mutationSignatures object. If NULL, COSMIC signatures will be retrieved
#' @param method distance method used to compute similarity (1 - distance)
#' @param threshold signal (similarity) upper threshold for maxing the signal
#' @param plot logical, shall a heatmap be plotted
#'
#' @return list, including distance matrix and a heatmap plot
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @importFrom proxy dist
#' @importFrom ggplot2 ggplot aes geom_tile scale_y_discrete scale_x_discrete theme_minimal labs scale_fill_gradient2 theme element_text element_blank margin
#' @importFrom stats quantile
#'
#' @export
matchSignatures <- function(mutSign,
                            reference = NULL,
                            method = 'cosine',
                            threshold = 0.5,
                            plot = "TRUE")
{
  # avoid NOTEs
  newSign <- NULL
  refSign <- NULL

  if(is.null(reference)){
    my.ref <- getCosmicSignatures()
  } else if (class(reference) == "mutationSignatures") {
    my.ref <- reference
  } else {
    stop("reference: wrong data type. Provide NULL or a 'mutationSignatures'-class object")
  }
  if(class(mutSign) != "mutationSignatures") {
    stop("mutSign: wrong data type. Provide a 'mutationSignatures'-class object")
  }
  if (sum(!getMutationTypes(mutSign) %in% getMutationTypes(my.ref)) != 0 |
      sum(!getMutationTypes(my.ref) %in% getMutationTypes(mutSign)) != 0)
    stop("non-compatible mutation types")

  # Coerce to tabular objects
  my.ref2 <- my.ref@mutationFreq
  if (ncol(my.ref2) < 2) { stop("You need to provide at least to reference signatures")}
  dimnames(my.ref2) <- list(my.ref@mutTypes[,1],
                            my.ref@signatureId[,1])

  mutSign2 <- mutSign@mutationFreq
  if (ncol(mutSign2) < 2) { stop("You need to provide at least two signatures to compare to the references")}
  dimnames(mutSign2) <- list(mutSign@mutTypes[,1],
                             mutSign@signatureId[,1])

  my.ref <- my.ref2
  mutSign <- mutSign2[rownames(my.ref), ]

  #Debug
  #message("Debug - proceeded! :-)")

  distMat <- sapply(1:ncol(mutSign), function(i){
    TMP <- cbind(tmpSig=mutSign[,i], my.ref)
    TMPdist <- proxy::dist(TMP, method = method, by_rows = FALSE)
    TMPdist[1:ncol(my.ref)]
  })

  dimnames(distMat) <- list(colnames(my.ref), colnames(mutSign))
  p <- NULL
  nuDF <- NULL

  if(plot) {

    nuDF <- do.call(rbind, lapply(1:nrow(distMat), function(i){
      do.call(rbind, lapply(1:ncol(distMat), function(j){
        data.frame(newSign = colnames(distMat)[j],
                   refSign = rownames(distMat)[i],
                   dist = distMat[i,j],
                   stringsAsFactors = FALSE)
      }))
    }))

    fillLims <- c(round(min(nuDF$dist)), round(max(nuDF$dist)))
    if (fillLims[1] != fillLims[2]) {
      nuDF$dist[nuDF$dist >= fillLims[2]] <- fillLims[2]
      nuDF$dist[nuDF$dist <= fillLims[1]] <- fillLims[1]
    } else if (fillLims[1] == 0 && fillLims[2] == 0) {
      fillLims <- c(0, 1)
    } else {
      msg00 <- paste0("Ouch! An error has occurred, this is embarassing! \n",
                      "Please, try selecting more signatures to compare, and run the function again. \n",
                      "If the error persists, please email the maintainer ",
                      "with a reproducible example of what has happened. Thank you.")

      stop(msg00)
    }

    # Quick cosine fix
    nuDF$dist[nuDF$dist < fillLims[1]] <- fillLims[1]
    nuDF$dist[nuDF$dist > fillLims[2]] <- fillLims[2]


    p <- ggplot2::ggplot(nuDF, ggplot2::aes(y=newSign, x=refSign))
    p <- p + ggplot2::geom_tile(ggplot2::aes(fill=dist), width=.875, height=.875)
    p <- p + ggplot2::scale_y_discrete(limits=base::rev(colnames(distMat)))
    p <- p + ggplot2::scale_x_discrete(limits=rownames(distMat))

    p <- p + ggplot2::theme_minimal(base_size = 11) + ggplot2::labs(x = "", y = "")
    p <- p + ggplot2::labs(title = "Signature Comparison")
    #p <- p + scale_fill_gradient2(low = "#f40000", mid = "#ffe9e9", high = "white",
    #                              midpoint = as.numeric(quantile(distMat, probs = 0.1, na.rm = TRUE)),
    #                              guide = "colourbar",
    #                              name = paste(method, "\ndistance", sep = ""))
    #p <- p + scale_fill_gradient(low = "red2", high = "white",
    #                             guide = "colourbar",
    #                             name = paste(method, "\ndistance", sep = ""))

    nuMid <- base::as.numeric(stats::quantile(distMat, probs = 0.2, na.rm = TRUE))
    #nuMid2 <- base::as.numeric(mean(distMat, na.rm = TRUE))
    nuMax <- fillLims
    # Try adjust
    #nuMid <- max(c(nuMid, 0.5), na.rm = TRUE)

    if (!is.null(threshold)) {
      if (!is.na(base::as.numeric(threshold[1]))){
        nuMid <- base::as.numeric(threshold[1])
      }
    } else if (fillLims[1] == 0 && fillLims[2] == 1) {
      if (nuMid < 0.425) { nuMid <- 0.425}
    }

    p <- p + ggplot2::scale_fill_gradient2(low = "#f40000", mid = "#ffe9e9", high = "white",
                                           midpoint = nuMid,
                                           guide = "colourbar",
                                           name = paste(method, "\ndistance", sep = ""),
                                           limits = fillLims)
    p <- p + ggplot2::theme(legend.position = c('right'), # position the legend in the upper left
                            legend.justification = 0, # anchor point for legend.position.
                            legend.text = ggplot2::element_text(size = 9, color = 'gray10'),
                            legend.title = ggplot2::element_text(size = 11),
                            plot.title = ggplot2::element_text(size = 13, face = 'bold', color = 'gray10', hjust = 0.5),
                            axis.text = ggplot2::element_text(face = 'bold'),
                            panel.grid.major.y = ggplot2::element_blank(),
                            panel.grid.major.x = ggplot2::element_blank(),
                            axis.text.x=ggplot2::element_text(angle = 90, hjust = 1 , vjust = 0.5,
                                                              margin=ggplot2::margin(-5,0,-15,0)),
                            axis.text.y=ggplot2::element_text(hjust = 1, vjust = 0.5,
                                                              margin = ggplot2::margin(0,3,0,0)))
    #print(p)
  }
  out <- list()
  out[["distanceMatrix"]] <- distMat
  out[["distanceDataFrame"]] <- nuDF
  out[["plot"]] <- p

  return(out)
}


## -----------------------

#' Process VCF data.
#'
#' Check, annotate, and process variants imported from a list of VCF files, so that it can be used
#' to run a mutational signature analysis
#'
#'
#' @param vcfData data.frame, includes mutation data from 2 or more samples
#' @param BSGenomeDb a BSGenomeDb-class object storing the genomic sequences and coordinates
#' @param chr_colName string, name of the column including the chromosome (seq) name. Defaults to "CHROM"
#' @param pos_colName string, name of the column including the genomic coordinates/position. Defaults to "POS"
#' @param ref_colName string, name of the column including the reference nucleotide. Defaults to "REF"
#' @param alt_colName string, name of the column including the variant nucleotide. Defaults to "ALT"
#' @param sample_colName string, name of the column including the sample ID. Can be NULL
#' @param nucl_contextN integer, span (in nucelotides) of the context around the variants. Defaults to 3
#' @param verbose logical, shall information about the ongoing analysis be printed to console
#'
#' @return a data.frame including processed variants from VCF files
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @export
processVCFdata <- function(vcfData,
                           BSGenomeDb,
                           chr_colName = "CHROM",
                           pos_colName = "POS",
                           ref_colName = "REF",
                           alt_colName = "ALT",
                           sample_colName = NULL,
                           nucl_contextN = 3,
                           verbose = TRUE)
{
  # make sure
  if (verbose)
    message("Processing VCF data: ", appendLF = FALSE)

  #Initialize
  reprepInput <- NULL

  # Fill reprepInput
  if (!is.null(sample_colName)){
    if (length(sample_colName) == 1 & sample_colName[1] %in% colnames(vcfData)){

      # compute how many unique samples
      unique.SID <- unique(vcfData[,sample_colName])
      reprepInput <- lapply(unique.SID, function(iID) {
        TMP <- vcfData[vcfData[,sample_colName] == iID, ]
        rownames(TMP)
        TMP
      })
    } else {
      stop("Bad input")
    }
  } else {
    #initialize input as list
    reprepInput <- list(vcfData)
  }

  # Triple check
  if(is.null(reprepInput))
    stop("Bad input")

  # Now loop
  bigOUT <- lapply(1:length(reprepInput), (function(jj){
    finalOut <- tryCatch({
      if(verbose)
        message(paste("[", jj, "/", length(reprepInput), "]", sep = ""), appendLF = F)

      tmpVCF <- suppressMessages({filterSNV(dataSet = reprepInput[[jj]],
                                            seq_colNames = c(ref_colName, alt_colName))})
      if (verbose)
        message(".", appendLF = FALSE)

      tmpVCF <- suppressMessages({attachContext(mutData = tmpVCF,
                                                chr_colName = chr_colName,
                                                start_colName = pos_colName,
                                                end_colName = pos_colName,
                                                nucl_contextN = 3,
                                                BSGenomeDb = BSGenomeDb)})
      if (verbose)
        message(".", appendLF = F)

      tmpVCF <- suppressMessages({removeMismatchMut(mutData = tmpVCF,
                                                    refMut_colName = ref_colName,
                                                    context_colName = "context",
                                                    refMut_format = "N")})
      if (verbose)
        message(".", appendLF = F)

      tmpVCF <- suppressMessages({attachMutType(mutData = tmpVCF,
                                                ref_colName = ref_colName,
                                                var_colName = alt_colName,
                                                var2_colName = alt_colName,
                                                context_colName = "context")})
      if (verbose)
        message(".", appendLF = F)

      tmpVCF
    }, error = function(e) NA)
  }))

  bigOUT <- do.call(rbind, bigOUT)

  if(verbose) {
    message("", appendLF = T)
    message("Done!", appendLF = T)
  }
  return(bigOUT)
}

#' Remove Mismatched Mutations.
#'
#' Remove mutation types that do not match the expected nucleotidic context.
#'
#' @param mutData data.frame including mutation data, as well as the nucleotide
#' context around the mutated position
#' @param refMut_colName string, name of the column storing REF and VAR data. Defaults to "N>N"
#' @param context_colName string, name of the column storing nucleotide context around the variant.
#' @param refMut_format string, format of mutation types. Defaults to "N>N"
#'
#' @return filtered data.frame
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("removeMismatchMut")
#' y <- mutSignatures:::removeMismatchMut(x,
#'                                        refMut_colName = "REF",
#'                                        context_colName = "context",
#'                                        refMut_format = "N")
#' y
#'
#' @export
removeMismatchMut <- function(mutData,
                              refMut_colName = "mutation",
                              context_colName = "context",
                              refMut_format = "N>N")
{
  # Validate input data and fields in mutData
  if(!((is.data.frame(mutData) | is.matrix(mutData) ) &
       sum(c(refMut_colName, context_colName) %in% colnames(mutData)) == 2))
    stop ("Issue with the input dataset. Make sure to feed in a data.frame or
          a matrix and double check the name of the fields pointing to chromosome
          name, start and end positions")
  if (! refMut_format %in% c("N>N", "N"))
    stop("Unrecognized refMut_format")

  output <- data.frame(mutData, stringsAsFactors = FALSE, row.names = NULL)
  colnames(output) <- colnames(mutData)
  if (refMut_format == "N>N") {
    output$mutSite.dnaSeq.Ref <- substr(mutData[,refMut_colName], 1, 1)
    output$mutSite.dnaSeq.Mut <- substr(mutData[,refMut_colName], 3, 3)
    refMut_colName <- "mutSite.dnaSeq.Ref"
  }

  pos.to.extr <- (0.5 * (nchar(output[,context_colName]) - 1)) + 1
  keep.id <- output[,refMut_colName] == substr(output[,context_colName], pos.to.extr, pos.to.extr)
  if (sum( keep.id == FALSE) > 0)
    message(paste("Removing", sum( keep.id == FALSE), "rows because of mismatch"))

  output <- output[keep.id, ]
  rownames(output) <- NULL
  return(output)
}

#' Resolve Mutation Signatures.
#'
#' If Mutation signatures are known (such as COSMIC signatures), we can
#' estimate the contribution of each signature in different samples.
#' This functions used a matrix of mutation counts and a matrix of mutation
#' signatures, and estimates Exposures to Mutational Signature of each sample.
#'
#' @param mutCountData object storing mutation counts
#' @param signFreqData object storing mutation signatures
#' @param byFreq logical, shall exposures be estimated on per_mille normalized counts
#'
#' @return a list of objects including data about exposures to mutational signatures
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#'
#' @examples
#'     x <- mutSignatures:::getTestRunArgs("resolveMutSignatures")
#'     y <- mutSignatures::resolveMutSignatures(mutCountData = x$muts, signFreqData = x$sigs)
#'     y
#'
#'
#' @export
resolveMutSignatures <- function(mutCountData,
                                 signFreqData,
                                 byFreq = TRUE )
{
  # process expected objects
  if (class(mutCountData) == "mutationCounts")
    mutCountData <- coerceObj(x = mutCountData, to = "data.frame")

  if (class(signFreqData) == "mutationSignatures")
    signFreqData <- coerceObj(x = signFreqData, to = "data.frame")

  if (!(sum(!rownames(mutCountData) %in% rownames(signFreqData)) == 0 &
        sum(! rownames(signFreqData) %in% rownames(mutCountData)) == 0)) {
    stop (paste("There is an issue with the mutTypes.",
                "MutTypes in the mutType Count Matrix",
                "have to match those in the signature",
                "Matrix... check rownames()"))
  }
  mutCountData <- mutCountData[rownames(signFreqData),]

  if (byFreq) {
    full.Y <- apply(mutCountData, 2, (function(clmn){
      1000 * clmn / sum(clmn)
    }))
  } else {
    full.Y <- base::as.matrix(mutCountData)
  }


  # record sum counts
  mutSums <- apply(mutCountData, 2, sum)

  # make sure we have frequencies
  my.signatures <- apply(signFreqData, 2, (function(clmn){
    clmn / sum(clmn)
  }))
  X <- base::as.matrix(my.signatures)

  # initialize results collector
  out <- list()

  res <- custom_fcnnls(mutCounts = full.Y, signatures = X)

  beta.hat <- data.frame(base::t(res$coef / ifelse (byFreq, 1000, 1)), stringsAsFactors = FALSE)

  colnames(beta.hat) <- colnames(signFreqData)
  rownames(beta.hat) <- colnames(mutCountData)
  out$results <- list()
  out$coeffs <- list()
  out$coeffs$beta.hat <- beta.hat
  out$coeffs$unexplained.mutNum <- round((1 - apply(beta.hat, 1, sum)) * mutSums, digits = 0)
  out$coeffs$unexplained.mutFrac <- out$coeffs$unexplained.mutNum  / mutSums

  if (byFreq) {
    for (i in 1:nrow(beta.hat)) {
      beta.hat[i,] <- beta.hat[i,] * mutSums[i]
    }
  }

  out$results$count.result <- mutSignatures::as.mutsign.exposures(x = beta.hat, samplesAsCols = FALSE)
  out$results$freq.result <- mutSignatures::as.mutsign.exposures(do.call(rbind, lapply(1:nrow(beta.hat), (function(jjj){
    beta.hat[jjj,] / sum(beta.hat[jjj,] )
  }))), samplesAsCols = FALSE)

  out$results$fitted <- res$fitted
  out$results$residuals <- res$residuals

  return(out)
}

#' Compute Reverse Complement sequences.
#'
#' Transform a DNA sequence into its reverse-complement sequence.
#' ALternatively, only the reverse sequence (or only the complement) can be returned.
#'
#' @param DNAseq character vector of DNA sequences
#' @param rev logical, shall the reverse sequence be computed
#' @param compl logical, shall the complementary sequence be computed
#'
#' @return a character vector including transformed DNA sequences
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' A <- c("TAACCG", "CTCGA", "CNNA")
#' mutSignatures::revCompl(A)
#'
#' @export
revCompl <- function (DNAseq,
                      rev = TRUE,
                      compl = TRUE)
{
  if (is.character(DNAseq)) {
    resultVect <- sapply(DNAseq, (function(seqString) {
      mySeq <- toupper(seqString)
      if (rev == TRUE) {
        mySeq <- paste(sapply(1:nchar(mySeq), (function(i) {
          pos <- nchar(mySeq) + 1 - i
          substr(mySeq, pos, pos)
        })), collapse = "")
      }
      if (compl == TRUE) {
        mySeq <- paste(sapply(1:nchar(mySeq), (function(i) {
          aBase <- substr(mySeq, i, i)
          complBase <- c("A", "T", "C", "G", "N", "A")
          names(complBase) <- c("T", "A", "G", "C", "N",
                                "U")
          returnBase <- complBase[aBase]
          if (is.na(returnBase)) {
            stop("Bad input")
          }
          returnBase
        })), collapse = "")
      }
      return(mySeq)
    }))
    names(resultVect) <- NULL
    return(resultVect)
  }
  else {
    warning("Bad input")
  }
}

#' Table Mutation Types by Sample.
#'
#' Prepare a molten data.frame starting from a mutation count matrix.
#' Mutation types (rows) are countes for each sample (cols). The results are returned in a
#' 3-column data.frame.
#'
#' @param dataMatrix a numeric matrix including mutation counts
#' @param rowLab string, name for the column that will be storing row IDs, typically sample IDs
#' @param colLab string, name for the column that will be storing column IDs, typically sample IDs
#' @param valueLab string, name for the column that will be storing mutation count values
#'
#' @return data.frame storing mutation counts by sample
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' A <- cbind(`A>G`=c(5,10),`A>T`=c(3,20),`A>C`=c(15,0))
#' rownames(A) = c("Smpl1", "Smpl2")
#' mutSignatures::table2df(A)
#'
#' @export
table2df <- function(dataMatrix,
                     rowLab = "sample",
                     colLab = "feature",
                     valueLab = "count")
{
  #
  if(!is.null(colnames(dataMatrix))) {
    names.X <- colnames(dataMatrix)
  } else {
    names.X <- 1:ncol(dataMatrix)
  }
  if(!is.null(rownames(dataMatrix))) {
    names.Y <- rownames(dataMatrix)
  } else {
    names.Y <- 1:nrow(dataMatrix)
  }
  #
  tmp <- do.call(rbind, lapply(1:nrow(dataMatrix), (function(i){
    t(sapply(1:ncol(dataMatrix), (function(j){
      c(names.Y[i], names.X[j], dataMatrix[i,j])
    })))
  })))
  tmp <- data.frame(tmp, stringsAsFactors = FALSE)
  if (is.numeric(dataMatrix[1,1]))
    tmp[,3] <- base::as.numeric(base::as.character(tmp[,3]))
  rownames(tmp) <- NULL
  colnames(tmp) <- c(rowLab, colLab, valueLab)
  return(tmp)
}

#' Attach Nucleotide Context.
#'
#' Retrieve the nucleotide context around each DNA variant based on the
#' genomic coordinates of the variant and a reference BSGenome database.
#'
#' @param mutData data.frame storing mutation data
#' @param BSGenomeDb a BSGenomeDb-class object, storing info about the genome of interest
#' @param chr_colName string, name of the column storing seqNames. Defaults to "chr"
#' @param start_colName string, name of the column storing start positions. Defaults to "start_position"
#' @param end_colName string, name of the column storing end positions. Defaults to "end_position"
#' @param nucl_contextN integer, the span of nucleotides to be retrieved around the variant. Defaults to 3
#' @param context_colName string, name of the column that will be storing the
#' nucleotide context. Defaults to "context"
#' @param skip_seqName_check logical, shall seqNames be checked to remove non-official chromosomes. 
#' Defaults to FALSE
#'
#'
#' @return a modified data.frame including the nucleotide context in a new column
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{GitHub Repo}: \url{https://github.com/dami82/mutSignatures/}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Sci Rep paper}, introducing mutS: \url{https://www.nature.com/articles/s41598-020-75062-0/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#'
#' @importFrom methods slotNames
#' @importFrom utils installed.packages data
#'
#' @export
attachContext <- function(mutData,
                          BSGenomeDb,
                          chr_colName = "chr",
                          start_colName = "start_position",
                          end_colName = "end_position",
                          nucl_contextN = 3,
                          context_colName = "context",
                          skip_seqName_check = FALSE)
{
  # Init Exec
  exec <- TRUE

  # Check if dependencies from BIOC are installed
  BiocX <- c("GenomicRanges", "IRanges", "BSgenome", "GenomeInfoDb")
  all.packs <- rownames(utils::installed.packages())
  if (sum(BiocX %in% all.packs) != length(BiocX)) {
    exec <- FALSE
  }

  # Check if BSGenomeDb dependencies from BIOC are installed
  if (!"BSgenome" %in% class(BSGenomeDb)) {
    stop("BSGenomeDb is not a BSgenome-class object")
  }

  if (exec) {

    attachContext.addON <- mutSignatures::mutSigData$.addON$attachContext.addON

    YY <- tryCatch({
      attachContext.addON(mutData = mutData,
                          BSGenomeDb = BSGenomeDb,
                          chr_colName = chr_colName,
                          start_colName = start_colName,
                          end_colName = end_colName,
                          nucl_contextN = nucl_contextN,
                          context_colName = context_colName, 
                          skip_seqName_check = skip_seqName_check)},

      error = function(e) {
        message("An error has occurred!")
        NULL})

    return(YY)

  } else {
    message("Sorry, this operation could not be executed!")
    message("In order to run the `attachContext()` function, the following libraries have to be installed from Bioconductor:")
    message("")
    for (xbi in BiocX) {
      message(paste0("  --> ", xbi))
    }
    message("")
    message("Missing Bioconductor libs:")
    message("  --> ", appendLF = FALSE)
    bxx <- paste(BiocX[!BiocX %in% all.packs], collapse = ", ")
    message(bxx, appendLF = FALSE)
    message("")

    # no return
  }
}


#' Count Mutation Types.
#'
#' Analyze a table (data.frame) including mutation counts. Count and aggregate Count Mutation Types.
#' If multiple samples are included in the same table, results are aggregated by samples.
#'
#' @param mutTable data.frame including mutation types and an optional sample ID column
#' @param mutType_colName string, name of the column storing mutTypes
#' @param sample_colName string, name of the column storing sample identifiers. Can be NULL
#'
#' @return a mutationCounts-class object
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' x <- mutSignatures:::getTestRunArgs("countMutTypes")
#' x
#' y <- mutSignatures::countMutTypes(mutTable = x,
#'                                   mutType_colName = "mutation",
#'                                   sample_colName = "sample")
#' y
#'
#'
#' @export
countMutTypes <- function (mutTable, mutType_colName = "mutType", sample_colName = NULL)
{
  if (!(is.data.frame(mutTable) | is.matrix(mutTable)))
    stop("Bad input")
  if (!mutType_colName %in% colnames(mutTable))
    stop("mutType field not found")
  if (!is.null(sample_colName)) {
    if (!sample_colName %in% colnames(mutTable))
      stop("sample_colName field not found")
  }
  mutType.labels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
                      "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "A[C>T]A",
                      "A[C>T]C", "A[C>T]G", "A[C>T]T", "A[T>A]A", "A[T>A]C",
                      "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G",
                      "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A",
                      "C[C>G]C", "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C",
                      "C[C>T]G", "C[C>T]T", "C[T>A]A", "C[T>A]C", "C[T>A]G",
                      "C[T>A]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
                      "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[C>A]A",
                      "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G",
                      "G[C>T]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
                      "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "G[T>G]A",
                      "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[C>A]A", "T[C>A]C",
                      "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C", "T[C>G]G",
                      "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A",
                      "T[T>C]C", "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C",
                      "T[T>G]G", "T[T>G]T")
  custPatt01 <- "^(A|C|G|T)\\[(A|C|G|T)>(A|C|G|T)\\](A|C|G|T)$"
  if (sum(regexpr(custPatt01, mutTable[, mutType_colName]) >
          0) != length(mutTable[, mutType_colName])){
    message("Non-standard format... MutTypes will be inferred.")

    mutType.labels <- sort(unique(mutTable[, mutType_colName]))

  } else {
    # fix reverse-complement if needed
    # only works for the standard three-nucletide system
    idx.to.fix <- which(!mutTable[, mutType_colName] %in% mutType.labels)
    if (length(idx.to.fix) > 0) {
      tmp <- mutTable[, mutType_colName][idx.to.fix]
      corrected.mutTypes <- sapply(1:length(tmp), (function(i) {
        out <- c(revCompl(substr(tmp[i], 7, 7)), "[", revCompl(substr(tmp[i],
                                                                      3, 3)), ">", revCompl(substr(tmp[i], 5, 5)),
                 "]", revCompl(substr(tmp[i], 1, 1)))
        paste(out, collapse = "", sep = "")
      }))
      mutTable[, mutType_colName][idx.to.fix] <- corrected.mutTypes
    }
  }

  if (sum(!mutTable[, mutType_colName] %in% mutType.labels) > 0)
    stop("Problem with the mutType... Please check the input")
  if (is.null(sample_colName)) {
    my.mutTypes <- mutTable[, mutType_colName]
    out.1 <- sapply(mutType.labels, (function(mtt) {
      sum(my.mutTypes == mtt)
    }))
    out.1 <- data.frame(cbind(sample = out.1))
  } else {
    unique.cases <- unique(mutTable[, sample_colName])
    out.1 <- lapply(unique.cases, (function(csid) {
      tmp.df <- mutTable[mutTable[, sample_colName] ==
                           csid, ]
      out.3 <- sapply(mutType.labels, (function(mtt) {
        sum(tmp.df[, mutType_colName] == mtt)
      }))
      out.3 <- data.frame(cbind(out.3))
      colnames(out.3) <- csid
      out.3
    }))
    out.1 <- data.frame(do.call(cbind, out.1), stringsAsFactors = FALSE)
    tryCatch({
      colnames(out.1) <- unique.cases
    }, error = function(e) NULL)
  }
  fOUT <- new(Class = "mutationCounts",
              x = out.1,
              muts = data.frame(mutTypes = rownames(out.1),
                                stringsAsFactors = FALSE),
              samples = data.frame(ID = colnames(out.1),
                                   stringsAsFactors = FALSE))
  return(fOUT)
}

#' Sort Data by Mutation Type.
#'
#' Reorder a mutationSignatures, mutationCounts, data.frame, or matrix object by
#' sorting entries by mutation type.
#'
#' @param x an object storing mutation count data
#'
#' @return an object of the same class as x, with entries sorted according to mutation types.
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#'
#' @examples
#' A <- data.frame(S1=1:5, S2=5:1, S3=1:5)
#' rownames(A) <- c("A[A>T]G", "A[C>G]G", "T[A>T]G", "T[C>G]T", "T[C>G]G")
#' mutSignatures::sortByMutations(A)
#'
#' @export
sortByMutations <- function(x) {

  if (class(x) == "mutationSignatures") {
    outClass <- "mutationSignatures"
    x1 <- coerceObj(x = x, to = "data.frame")
  } else if (class(x) == "mutationCounts") {
    outClass <- "mutationCounts"
    x1 <- coerceObj(x = x, to = "data.frame")
  } else if (class(x) %in% c("matrix","data.frame")){
    outClass <- "data.frame"
    x1 <- x
  } else {
    stop ("Bad input")
  }

  # Get rownames == mutation names
  if (is.null(rownames(x1)))
    stop("Bad input")

  aRN <- rownames(x1)
  aRN <- sapply(sort(unique(sub("^.*\\[(.*)\\].*$", "\\1", aRN))),
                function(mt) {
                  sort(aRN[grepl(mt, aRN)])
                }, USE.NAMES = FALSE, simplify = FALSE)
  aRN <- do.call(c, aRN)
  ox <- x1[aRN,]

  # Format and return
  if (outClass == "data.frame") {
    return(ox)
  } else if (outClass == "mutationSignatures") {
    return(mutSignatures::as.mutation.signatures(ox))
  } else if (outClass == "mutationCounts") {
    return(mutSignatures::as.mutation.counts(ox))
  }
}

#' Simplify Mutational Signatures.
#'
#' This function is useufl when working with non-standard muation types, such as
#' tetra-nnucleotide mutation types or mutation types with long/complex context.
#' THe goal of this function is to aggregated together mutations that can be simplified
#' because of a common mutation core.
#' For example, mutation types AA[A>T]A, TA[A>T]A, CA[A>T]A, and GA[A>T]A can be simplified to the
#' core tri-nucleotide mutation A[A>T]A. THis function identifies mergeable mutation types,
#' and aggregates the corresponding counts/freqs.
#'
#' @param x a mutationSignatures-class object
#' @param asMutationSignatures logical, shall the results be returned as a mutationSignatures-class object
#'
#' @return  object including simplified mutational signatures data
#'
#' @details This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' More information and examples about mutational signature analysis can be found here:
#' \enumerate{
#'   \item \bold{Official website}: \url{http://www.mutsignatures.org}
#'   \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'   \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
#'  }
#'
#' @examples
#' A <- data.frame(Sig1=1:5, Sig2=5:1, Sig3=1:5)
#' A <- A/apply(A, 2, sum)
#' rownames(A) <- c("AA[C>A]A", "CA[C>A]A", "TA[C>A]A", "TA[C>G]A", "A[C>G]AT")
#' A <- mutSignatures::as.mutation.signatures(A)
#' mutSignatures::simplifySignatures(x = A)
#'
#'
#' @export
simplifySignatures <- function(x, asMutationSignatures = TRUE) {
  if(class(x) != "mutationSignatures")
    stop("Bad input")
  curMT <- getMutationTypes(x)
  if(sum(!grepl("[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]", curMT)) > 0)
    stop("Mutation Types are not compatible with simplifySignatures")

  # Start by defined tri-nucleotides
  mutType.labels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
                      "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
                      "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
                      "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                      "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
                      "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
                      "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
                      "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
                      "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")

  # Needs sorting
  input <- coerceObj(x = x, to = "data.frame")

  # sigNames <- x@signatureId$ID
  sigNames <- getSignatureIdentifiers(x)

  EXPLD <- sapply(mutType.labels, (function(pat){
    tpat <- sub("\\[", "\\\\[",
                sub("\\]", "\\\\]", pat))
    input[grepl(tpat, rownames(input)),]
  }), simplify = FALSE, USE.NAMES = TRUE)

  OUT <- base::as.data.frame(base::t(sapply(EXPLD, function(xi) {
    base::as.numeric(apply(xi, 2, sum, na.rm = TRUE))
  })), row.names = mutType.labels)
  colnames(OUT) <- sigNames

  if (asMutationSignatures){
    return(mutSignatures::as.mutation.signatures(OUT))
  }

  altOut <- lapply(1:ncol(OUT), (function(j){
    curSGN <- base::as.numeric(OUT[,j])
    iSig <- data.frame(curSGN, row.names = mutType.labels)
    colnames(iSig) <- sigNames[j]

    AttDF <- do.call(rbind, sapply(mutType.labels, function(mmd) {
      zi <- EXPLD[[mmd]]
      nuNames <- gsub("[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]", "...", rownames(zi))
      nuNums <- base::as.data.frame(rbind(zi[,j]))
      names(nuNums) <- nuNames
      nuNums
    }, simplify = FALSE, USE.NAMES = TRUE))

    # out
    cbind(iSig, AttDF)
  }))
  names(altOut) <- sigNames
  return(altOut)
}






