\name{setMutClusterParams}
\alias{setMutClusterParams}
\title{
Set Mutation Cluster Parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setMutClusterParams(num_processesToExtract = 2, num_totIterations = 10, num_parallelCores = 1, thresh_removeWeakMutTypes = 0.01, thresh_removeLastPercent = 0.07, distanceFunction = "cosine", num_totReplicates = 100, eps = 2.2204e-16, stopconv = 20000, niter = 1e+06, guided = TRUE, debug = FALSE, approach = "freq", stopRule = "DF", algorithm = "brunet", logIterations = "lite")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num_processesToExtract}{
%%     ~~Describe \code{num_processesToExtract} here~~
}
  \item{num_totIterations}{
%%     ~~Describe \code{num_totIterations} here~~
}
  \item{num_parallelCores}{
%%     ~~Describe \code{num_parallelCores} here~~
}
  \item{thresh_removeWeakMutTypes}{
%%     ~~Describe \code{thresh_removeWeakMutTypes} here~~
}
  \item{thresh_removeLastPercent}{
%%     ~~Describe \code{thresh_removeLastPercent} here~~
}
  \item{distanceFunction}{
%%     ~~Describe \code{distanceFunction} here~~
}
  \item{num_totReplicates}{
%%     ~~Describe \code{num_totReplicates} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{stopconv}{
%%     ~~Describe \code{stopconv} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{guided}{
%%     ~~Describe \code{guided} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{approach}{
%%     ~~Describe \code{approach} here~~
}
  \item{stopRule}{
%%     ~~Describe \code{stopRule} here~~
}
  \item{algorithm}{
%%     ~~Describe \code{algorithm} here~~
}
  \item{logIterations}{
%%     ~~Describe \code{logIterations} here~~
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
function (num_processesToExtract = 2, num_totIterations = 10, 
    num_parallelCores = 1, thresh_removeWeakMutTypes = 0.01, 
    thresh_removeLastPercent = 0.07, distanceFunction = "cosine", 
    num_totReplicates = 100, eps = 2.2204e-16, stopconv = 20000, 
    niter = 1e+06, guided = TRUE, debug = FALSE, approach = "freq", 
    stopRule = "DF", algorithm = "brunet", logIterations = "lite") 
{
    paramList <- list()
    if (!((is.numeric(num_processesToExtract[1]) & num_processesToExtract[1] > 
        0))) 
        stop("Provide a reasonable number of signatures/processes to extract")
    paramList$num_processesToExtract <- round(num_processesToExtract[1])
    if (!(is.numeric(num_totIterations[1]) & num_totIterations[1] > 
        0)) 
        stop("Provide a reasonable number of iterations to run (Bootstrapping)")
    paramList$num_totIterations <- round(num_totIterations[1])
    if (!(is.numeric(num_parallelCores[1]) & num_parallelCores[1] > 
        0)) 
        stop("Provide a reasonable number of CPU cores to use for the analysis")
    paramList$num_parallelCores <- round(num_parallelCores[1])
    if (!(is.numeric(thresh_removeWeakMutTypes[1]) & thresh_removeWeakMutTypes[1] >= 
        0 & thresh_removeWeakMutTypes[1] < 1)) 
        stop("Provide a reasonable (0.00-0.99) number of low-occurring mutation types to remove from the input before starting the analysis")
    paramList$thresh_removeWeakMutTypes <- thresh_removeWeakMutTypes[1]
    if (!(is.numeric(thresh_removeLastPercent[1]) & thresh_removeLastPercent[1] >= 
        0 & thresh_removeLastPercent[1] < 1)) 
        stop("Provide a reasonable (0.00-0.99) number for filtering out poor iterations")
    paramList$thresh_removeLastPercent <- thresh_removeLastPercent[1]
    allowed.dist.methods <- c("Braun-Blanquet", "Chi-squared", 
        "correlation", "cosine", "Cramer", "Dice", "eDice", "eJaccard", 
        "Fager", "Faith", "Gower", "Hamman", "Jaccard", "Kulczynski1", 
        "Kulczynski2", "Michael", "Mountford", "Mozley", "Ochiai", 
        "Pearson", "Phi", "Phi-squared", "Russel", "simple matching", 
        "Simpson", "Stiles", "Tanimoto", "Tschuprow", "Yule", 
        "Yule2", "Bhjattacharyya", "Bray", "Canberra", "Chord", 
        "divergence", "Euclidean", "fJaccard", "Geodesic", "Hellinger", 
        "Kullback", "Levenshtein", "Mahalanobis", "Manhattan", 
        "Minkowski", "Podani", "Soergel", "supremum", "Wave", 
        "Whittaker")
    if (!(is.character(distanceFunction[1]) & distanceFunction[1] \%in\% 
        allowed.dist.methods)) 
        stop("Unknown method for calculating distances. For options, run: <<summary(proxy::pr_DB)>>")
    paramList$distanceFunction <- distanceFunction[1]
    if (!(is.numeric(num_totReplicates[1]) & num_totReplicates[1] > 
        99)) 
        stop("Provide a reasonable number of replicates for stability evaluation of the results")
    paramList$num_totReplicates <- round(num_totReplicates[1])
    if (!(is.numeric(eps[1]) & eps[1] > 0 & eps[1] < 1e-04)) 
        stop("Provide a reasonably small number (0 < n < 0.0001) for data overflow prevention")
    paramList$eps <- eps[1]
    if (!(is.numeric(stopconv[1]) & stopconv[1] >= 500)) 
        stop("Provide a reasonable large number: number of 'conn-matrix-stable' iterations before stopping NMF")
    paramList$stopconv <- round(stopconv[1])
    if (!(is.numeric(niter[1]) & niter[1] >= 20000)) 
        stop("Provide a reasonable large number: total NMF iterations")
    paramList$niter <- round(niter[1])
    paramList$guided <- ifelse(guided, TRUE, FALSE)
    paramList$debug <- ifelse(debug, TRUE, FALSE)
    paramList$approach <- ifelse(approach == "counts", "counts", 
        "freq")
    paramList$stopRule <- ifelse(stopRule == "LA", "LA", "DF")
    paramList$algorithm <- ifelse(tolower(algorithm) \%in\% c("brunet", 
        "alexa"), "brunet", "chihjen")
    paramList$logIterations <- ifelse(tolower(logIterations) \%in\% 
        c("lite", "light", "li"), "lite", "full")
    return(new(Class = "mutFrameworkParams", params = paramList))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
