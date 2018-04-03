\name{removeMismatchMut}
\alias{removeMismatchMut}
\title{
Remove Mutations Corresponding to Mismatched Reference Positions 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
removeMismatchMut(mutData, refMut_colName = "mutation", context_colName = "context", refMut_format = "N>N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutData}{
%%     ~~Describe \code{mutData} here~~
}
  \item{refMut_colName}{
%%     ~~Describe \code{refMut_colName} here~~
}
  \item{context_colName}{
%%     ~~Describe \code{context_colName} here~~
}
  \item{refMut_format}{
%%     ~~Describe \code{refMut_format} here~~
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
function (mutData, refMut_colName = "mutation", context_colName = "context", 
    refMut_format = "N>N") 
{
    if (!((is.data.frame(mutData) | is.matrix(mutData)) & sum(c(refMut_colName, 
        context_colName) \%in\% colnames(mutData)) == 2)) 
        stop("Issue with the input dataset. Make sure to feed in a data.frame or\n          a matrix and double check the name of the fields pointing to chromosome\n          name, start and end positions")
    if (!refMut_format \%in\% c("N>N", "N")) 
        stop("Unrecognized refMut_format")
    output <- data.frame(mutData, stringsAsFactors = FALSE, row.names = NULL)
    colnames(output) <- colnames(mutData)
    if (refMut_format == "N>N") {
        output$mutSite.dnaSeq.Ref <- substr(mutData[, refMut_colName], 
            1, 1)
        output$mutSite.dnaSeq.Mut <- substr(mutData[, refMut_colName], 
            3, 3)
        refMut_colName <- "mutSite.dnaSeq.Ref"
    }
    pos.to.extr <- (0.5 * (nchar(output[, context_colName]) - 
        1)) + 1
    keep.id <- output[, refMut_colName] == substr(output[, context_colName], 
        pos.to.extr, pos.to.extr)
    if (sum(keep.id == FALSE) > 0) 
        message(paste("Removing", sum(keep.id == FALSE), "rows because of mismatch"))
    output <- output[keep.id, ]
    rownames(output) <- NULL
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
