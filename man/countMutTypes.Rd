\name{countMutTypes}
\alias{countMutTypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count Mutation Types
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
countMutTypes(mutTable, mutType_colName = "mutType", sample_colName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutTable}{
%%     ~~Describe \code{mutTable} here~~
}
  \item{mutType_colName}{
%%     ~~Describe \code{mutType_colName} here~~
}
  \item{sample_colName}{
%%     ~~Describe \code{sample_colName} here~~
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
function (mutTable, mutType_colName = "mutType", sample_colName = NULL) 
{
    if (!(is.data.frame(mutTable) | is.matrix(mutTable))) 
        stop("Bad input")
    if (!mutType_colName \%in\% colnames(mutTable)) 
        stop("mutType field not found")
    if (!is.null(sample_colName)) {
        if (!sample_colName \%in\% colnames(mutTable)) 
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
        0) != length(mutTable[, mutType_colName])) 
        stop("Problem with the mutation type format... Please use the following format: A[C>G]T")
    idx.to.fix <- which(!mutTable[, mutType_colName] \%in\% mutType.labels)
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
    if (sum(!mutTable[, mutType_colName] \%in\% mutType.labels) > 
        0) 
        stop("Problem with the mutType... Please check the input")
    if (is.null(sample_colName)) {
        my.mutTypes <- mutTable[, mutType_colName]
        out.1 <- sapply(mutType.labels, (function(mtt) {
            sum(my.mutTypes == mtt)
        }))
        out.1 <- data.frame(cbind(sample = out.1))
    }
    else {
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
    fOUT <- new(Class = "mutationCounts", x = out.1, muts = data.frame(mutTypes = rownames(out.1), 
        stringsAsFactors = FALSE), samples = data.frame(ID = colnames(out.1), 
        stringsAsFactors = FALSE))
    return(fOUT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
