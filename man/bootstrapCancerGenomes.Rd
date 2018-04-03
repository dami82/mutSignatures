\name{bootstrapCancerGenomes}
\alias{bootstrapCancerGenomes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Cancer Genomes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapCancerGenomes(genomes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomes}{
%%     ~~Describe \code{genomes} here~~
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
function (genomes) 
{
    genome.col.sums <- apply(genomes, 2, sum)
    norm.genomes <- genomes/matrix(genome.col.sums, ncol = ncol(genomes), 
        nrow = nrow(genomes), byrow = TRUE)
    bootstrapGenomes <- sapply(1:length(genome.col.sums), (function(i) {
        stats::rmultinom(1, size = genome.col.sums[i], prob = norm.genomes[, 
            i])
    }))
    return(bootstrapGenomes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
