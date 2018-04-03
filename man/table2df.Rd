\name{table2df}
\alias{table2df}
\title{
Convert Matrix to Data Frame}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
table2df(dataMatrix, rowLab = "sample", colLab = "feature", valueLab = "count")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrix}{
%%     ~~Describe \code{dataMatrix} here~~
}
  \item{rowLab}{
%%     ~~Describe \code{rowLab} here~~
}
  \item{colLab}{
%%     ~~Describe \code{colLab} here~~
}
  \item{valueLab}{
%%     ~~Describe \code{valueLab} here~~
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
function (dataMatrix, rowLab = "sample", colLab = "feature", 
    valueLab = "count") 
{
    if (!is.null(colnames(dataMatrix))) {
        names.X <- colnames(dataMatrix)
    }
    else {
        names.X <- 1:ncol(dataMatrix)
    }
    if (!is.null(rownames(dataMatrix))) {
        names.Y <- rownames(dataMatrix)
    }
    else {
        names.Y <- 1:nrow(dataMatrix)
    }
    tmp <- do.call(rbind, lapply(1:nrow(dataMatrix), (function(i) {
        t(sapply(1:ncol(dataMatrix), (function(j) {
            c(names.Y[i], names.X[j], dataMatrix[i, j])
        })))
    })))
    tmp <- data.frame(tmp, stringsAsFactors = FALSE)
    if (is.numeric(dataMatrix[1, 1])) 
        tmp[, 3] <- as.numeric(as.character(tmp[, 3]))
    rownames(tmp) <- NULL
    colnames(tmp) <- c(rowLab, colLab, valueLab)
    return(tmp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
