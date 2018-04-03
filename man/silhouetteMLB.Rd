\name{silhouetteMLB}
\alias{silhouetteMLB}
\title{
Run a Silhouette Analysis Similar to Matlab
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
silhouetteMLB(data, fac, method = "cosine", plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{fac}{
%%     ~~Describe \code{fac} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (data, fac, method = "cosine", plot = TRUE) 
{
    if (nrow(data) != length(fac)) 
        stop("Bad input!")
    dist.matrix <- as.matrix(proxy::dist(x = data, method = method))
    sil.check <- cluster::silhouette(x = as.numeric(as.factor(fac)), 
        dist = dist.matrix)
    if (plot == TRUE) {
        tmp <- lapply(unique(sil.check[, 1]), (function(clid) {
            part.out <- sil.check[sil.check[, 1] == clid, ]
            part.out[order(part.out[, 3], decreasing = TRUE), 
                ]
        }))
        tmp <- do.call(rbind, tmp)
        xrange <- c(min(tmp[, 3]), max(tmp[, 3]))
        xrange[1] <- ifelse(xrange[1] > 0, 0, (-1.15) * abs(xrange[1]))
        xrange[2] <- 1.15
        graphics::barplot(tmp[nrow(tmp):1, 3], col = as.factor(tmp[nrow(tmp):1, 
            1]), xlim = xrange, horiz = TRUE, xlab = "Silhouette Value", 
            ylab = "", main = "Silhouette Plot", border = as.factor(tmp[nrow(tmp):1, 
                1]))
        graphics::abline(v = 0)
        graphics::title(ylab = "Iter. Results (by Group)", line = 1, 
            cex.lab = 1, font = 2)
    }
    return(as.vector(sil.check[, 3]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
