\name{filterSNV}
\alias{filterSNV}
\title{
Filter Single Nucleotide Variants
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filterSNV(dataSet, seq_colNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
%%     ~~Describe \code{dataSet} here~~
}
  \item{seq_colNames}{
%%     ~~Describe \code{seq_colNames} here~~
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
function (dataSet, seq_colNames) 
{
    if (!(is.data.frame(dataSet) | is.matrix(dataSet)) | sum(!seq_colNames \%in\% 
        colnames(dataSet)) > 0 | length(seq_colNames) < 2) {
        stop("Bad input or seq_colNames not found")
    }
    check.tab <- sapply(1:length(seq_colNames), (function(i) {
        tmp <- gsub("[[:space:]]", "", dataSet[, seq_colNames[i]])
        toupper(tmp) \%in\% c("A", "C", "G", "T")
    }))
    toKeep <- apply(check.tab, 1, (function(rw) {
        sum(rw) == length(rw)
    }))
    out <- dataSet[toKeep, ]
    rownames(out) <- NULL
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
