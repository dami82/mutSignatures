\name{getTestRunArgs}
\alias{getTestRunArgs}
\title{
Build mutSignature Objects for Examples
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getTestRunArgs(testN = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{testN}{
%%     ~~Describe \code{testN} here~~
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
function (testN = 1) 
{
    out <- list()
    if (testN == 3) {
        out$v <- cbind(c(142, 133, 1, 24, 53, 55, 4, 4, 100), 
            c(132, 113, 0, 34, 50, 52, 3, 3, 17), c(155, 139, 
                10, 14, 53, 45, 2, 5, 13), c(124, 156, 22, 21, 
                52, 45, 2, 7, 100))
        out$r <- 2
        out$params <- setMutClusterParams(num_processesToExtract = 2)
        out$params$num_parallelCores <- 1
        out$params$stopconv <- 800
        out$params$niter <- 8000
    }
    else if (testN == 4) {
        out$data <- cbind(c(runif(4, 10, 25), runif(6, 20, 50), 
            runif(3, 0, 5)), c(runif(4, 50, 60), runif(6, 45, 
            55), runif(3, 30, 40)), c(runif(4, 12, 15), runif(6, 
            10, 15), runif(3, 10, 12)), c(runif(4, 5, 20), runif(6, 
            16, 26), runif(3, 24, 29)))
        out$fac <- c(rep(1, 4), rep(2, 6), rep(3, 3))
    }
    else if (testN == 5) {
        out$W <- do.call(cbind, lapply(1:20, (function(i) {
            cbind(c(runif(4, 0.05, 0.15), c(1e-15 * runif(1, 
                1, 9)), runif(3, 0.003, 0.007), runif(9, 0.04, 
                0.09)), c(runif(3, 0.08, 0.18), c(1e-15 * runif(1, 
                1, 9)), runif(1, 0.03, 0.07), runif(3, 0.004, 
                0.009), runif(9, 0.04, 0.09)))
        })))
        out$H <- do.call(rbind, lapply(1:40, (function(i) {
            c(runif(1, 0.005, 0.099), runif(2, 50, 800), runif(2, 
                0.005, 0.099), runif(1, 50, 750), runif(2, 0.005, 
                0.099), runif(1, 20, 800))
        })))
        out$params <- setMutClusterParams(num_processesToExtract = 2, 
            num_totIterations = 7, num_parallelCores = 3)
        out$params$analyticApproach <- "denovo"
    }
    else if (testN == 6) {
        tmut <- runif(10, 150, 1350)
        eff1 <- runif(10, 0.45, 0.89)
        eff2 <- 1 - eff1
        out$exposures <- sapply(1:10, (function(i) {
            c(eff1[i], eff2[i]) * tmut[i]
        }))
    }
    else {
        set.seed(999)
        my.mat <- sapply(1:10, (function(i) {
            c(as.integer(runif(3, 80, 150)), as.integer(runif(7, 
                0, 10)), as.integer(runif(9, 40, 80)), as.integer(runif(1, 
                0, 3)), as.integer(runif(10, 60, 120)))
        }))
        rownames(my.mat) <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", 
            "A[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", 
            "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "A[T>A]A", 
            "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C", 
            "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", 
            "A[T>G]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", 
            "C[C>G]A", "C[C>G]C")
        out$mutCount.obj <- setMutCountObject(mutCountMatrix = my.mat)
        out$params <- setMutClusterParams(num_processesToExtract = 2)
        out$params$num_totIterations <- 1
        out$params$num_parallelCores <- 1
        out$params$stopconv <- 800
        out$params$niter <- 8000
        if (testN == 2) {
            out$params$stopconv <- 2000
            out$params$niter <- 15000
            out$params$num_totIterations <- 3
        }
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
