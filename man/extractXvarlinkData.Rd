\name{extractXvarlinkData}
\alias{extractXvarlinkData}
\title{
Extract Data From XVarLink Fields
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractXvarlinkData(xvarLink_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xvarLink_data}{
%%     ~~Describe \code{xvarLink_data} here~~
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
function (xvarLink_data) 
{
    tmpVars <- sub("^.*&var=[[:alnum:]]+(,.*)&.*$", "\\1", xvarLink_data)
    tmpVars[!grepl("^,.+", tmpVars)] <- NA
    tmpVars <- strsplit(tmpVars, ",")
    tmpVars <- do.call(rbind, lapply(1:length(tmpVars), function(i) {
        if (tmpVars[[i]][1] == "" & length(tmpVars[[i]]) == 5) {
            tmpVars[[i]][2:5]
        }
        else {
            c(NA, NA, NA, NA)
        }
    }))
    out <- data.frame(chrXvar = as.character(tmpVars[, 1]), posXvar = as.numeric(tmpVars[, 
        2]), refXvar = as.character(tmpVars[, 3]), mutXvar = as.character(tmpVars[, 
        4]), stringsAsFactors = FALSE)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
