\name{chihJenNMF}
\alias{chihJenNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alternative NMF Algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
chihJenNMF(v, r, params)
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
    v <- as.matrix(v)
    rownames(v) <- NULL
    colnames(v) <- NULL
    if (min(v) < 0) 
        stop("Matrix entries cannot be negative")
    if (min(apply(v, 1, sum)) == 0) 
        stop("Entries cannot all be equal to 0")
    if (debug) 
        graphics::plot(-10, xlim = c(1000, niter), ylim = c((0.1 * 
            err.threshold), 10), log = "xy", xlab = "Iteration", 
            ylab = "Variation", main = "Convergence")
    W.k <- do.call(cbind, lapply(1:num.processes, (function(i) {
        out.o <- runif(n = nrow(v), min = eps, max = 100)
        out.o/sum(out.o)
    })))
    H.k <- matrix((1/num.processes), nrow = num.processes, ncol = ncol(v))
    itr <- 1
    chk.j <- 1
    cons.old <- as.vector(W.k)
    stationary.chk <- 0
    force.out <- 1
    while (itr < niter) {
        if (itr\%\%dot.eachSteps == 0) {
            if (stationary.chk > chk.step) {
                message(":", appendLF = FALSE)
            }
            else {
                message(".", appendLF = FALSE)
            }
        }
        WtW <- t(W.k) \%*\% W.k
        gradH <- ((WtW \%*\% H.k) - (t(W.k) \%*\% v))
        H.b <- H.k
        H.b[H.b < eps] <- eps
        H.k <- H.k - (H.b/((WtW \%*\% H.b) + eps)) * gradH
        HHt <- H.k \%*\% t(H.k)
        gradW <- (W.k \%*\% HHt) - (v \%*\% t(H.k))
        W.b <- W.k
        W.b[W.b < eps] <- eps
        W.k <- W.k - (W.b/((W.b \%*\% HHt) + eps)) * gradW
        S <- apply(W.k, 2, sum)
        W.k <- do.call(cbind, lapply(1:ncol(W.k), (function(ci) {
            W.k[, ci]/S[ci]
        })))
        H.k <- do.call(rbind, lapply(1:nrow(H.k), (function(ri) {
            H.k[ri, ] * S[ri]
        })))
        H.k <- do.call(cbind, lapply(1:ncol(H.k), (function(ci) {
            H.k[, ci]/sum(H.k[, ci])
        })))
        H.k[H.k < eps] <- eps
        W.k[W.k < eps] <- eps
        if (itr > stopconv & itr\%\%chk.step == 0) {
            chk.j <- chk.j + 1
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
