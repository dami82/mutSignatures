\name{alexaNMF}
\alias{alexaNMF}
\title{
Perform NMF According to Brunet Algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alexaNMF(v, r, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
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
function (v, r, params) 
{
    debug <- params$debug
    chk.step <- 50
    dot.eachSteps <- 2000
    eps <- params$eps
    num.processes <- r
    stopconv <- params$stopconv
    niter <- params$niter
    err.threshold <- 1e-10
    stopRule <- ifelse(params$stopRule == "LA", "LA", "DF")
    v <- as.matrix(v)
    rownames(v) <- NULL
    colnames(v) <- NULL
    if (min(v) < 0) 
        stop("Matrix entries cannot be negative")
    if (min(apply(v, 1, sum)) == 0) 
        stop("Entries cannot all be equal to 0")
    W.k <- do.call(cbind, lapply(1:num.processes, (function(i) {
        out.o <- runif(n = nrow(v), min = eps, max = 100)
        out.o/sum(out.o)
    })))
    H.k <- matrix((1/num.processes), nrow = num.processes, ncol = ncol(v))
    itr <- 1
    chk.j <- 1
    stationary.chk <- 0
    force.out <- 1
    if (debug) 
        graphics::plot(-10, xlim = c(1000, niter), ylim = c(ifelse(stopRule == 
            "DF", (0.1 * err.threshold), 0.001), ifelse(stopRule == 
            "DF", 10, ncol(H.k))), log = "xy", xlab = "Iteration", 
            ylab = "Variation", main = "Convergence")
    cons.old <- as.vector(W.k)
    consold <- matrix(0, nrow = ncol(H.k), ncol = ncol(H.k))
    while (itr < niter) {
        if (itr\%\%dot.eachSteps == 0) {
            if (stationary.chk > chk.step) {
                message(":", appendLF = FALSE)
            }
            else {
                message(".", appendLF = FALSE)
            }
        }
        delta.01 <- apply(W.k, 2, sum)
        H.k <- H.k * (t(W.k) \%*\% (v/(W.k \%*\% H.k)))/delta.01
        H.k[H.k < eps] <- eps
        W.tmp <- W.k * ((v/(W.k \%*\% H.k)) \%*\% t(H.k))
        W.k <- do.call(cbind, lapply(1:ncol(W.tmp), (function(ci) {
            W.tmp[, ci]/sum(H.k[ci, ])
        })))
        W.k[W.k < eps] <- eps
        if (itr > stopconv & itr\%\%chk.step == 0 & stopRule == 
            "DF") {
            chk.j <- chk.j + 1
            H.k[H.k < eps] <- eps
            W.k[W.k < eps] <- eps
            cons <- as.vector(W.k)
            dist.measure <- proxy::dist(rbind(cons, cons.old), 
                method = "cosine")[1]
            cons.old <- cons
            if (debug) 
                points(itr, (dist.measure + (err.threshold * 
                  0.1)), pch = 19, col = "red2")
            if (dist.measure < err.threshold) {
                stationary.chk <- stationary.chk + 1
            }
            else {
                stationary.chk <- 0
            }
            if (stationary.chk > (stopconv/chk.step)) {
                force.out <- 0
                message("$", appendLF = FALSE)
                (break)()
            }
        }
        else if (itr > stopconv & itr\%\%chk.step == 0 & stopRule == 
            "LA") {
            chk.j <- chk.j + 1
            H.k[H.k < eps] <- eps
            W.k[W.k < eps] <- eps
            y <- apply(H.k, 2, max)
            index <- apply(H.k, 2, (function(dt) {
                which.max(dt)[1]
            }))
            mat1 = t(sapply(1:ncol(H.k), (function(ii) {
                index
            })))
            mat2 = sapply(1:ncol(H.k), (function(ii) {
                index
            }))
            cons <- mat1 == mat2
            if (sum(cons != consold) == 0) {
                stationary.chk <- stationary.chk + 1
            }
            else {
                stationary.chk <- 0
            }
            consold <- cons
            if (debug) 
                points(itr, (sum(cons != consold)/100), pch = 19, 
                  col = "red2")
            if (stationary.chk > (stopconv/chk.step)) {
                force.out <- 0
                message("$", appendLF = FALSE)
                (break)()
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
