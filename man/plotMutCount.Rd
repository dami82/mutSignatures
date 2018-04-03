\name{plotMutCount}
\alias{plotMutCount}
\title{
Plot Mutation Counts
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotMutCount(mutCount, top = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutCount}{
%%     ~~Describe \code{mutCount} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
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
function (mutCount, top = 50) 
{
    if (class(mutCount) == "mutationExposure") {
        mutCount <- as.data.frame(mutCount)
    }
    rownames(mutCount) <- NULL
    mutCount <- mutCount[order(apply(mutCount, 1, sum), decreasing = TRUE), 
        ]
    rownames(mutCount) <- NULL
    if (is.null(top)) {
        mutDF <- table2df(dataMatrix = mutCount)
    }
    else if (is.na(top[1]) | top[1] > nrow(mutCount) | top[1] < 
        2) {
        mutDF <- table2df(dataMatrix = mutCount)
    }
    else {
        mutDF <- table2df(dataMatrix = mutCount[1:top, ])
    }
    mutDF$sample <- 1e+05 + as.numeric(as.character(mutDF$sample))
    mutDF$sample <- as.character(mutDF$sample)
    mutDF$feature <- factor(mutDF$feature, levels = rev(colnames(mutCount)))
    bp <- ggplot(data = mutDF, aes(x = sample, y = count, fill = feature)) + 
        geom_bar(stat = "identity")
    bp <- bp + theme_minimal() + theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + theme(axis.line.y = element_line(colour = "black", 
        size = 0.75), axis.line.x = element_line(colour = "black", 
        size = 0.75), axis.ticks.y = element_line(colour = "black", 
        size = 1), axis.ticks.length = unit(x = 6, "points"), 
        plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 
        0), limits = c(0, 1.2 * max(apply(mutCount, 1, function(rx) sum(rx, 
        na.rm = TRUE)))))
    return(bp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
