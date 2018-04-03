\name{revCompl}
\alias{revCompl}
\title{
Reverse Complement of a Nucleic Acid Sequence
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
revCompl(DNAseq, rev = TRUE, compl = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DNAseq}{
%%     ~~Describe \code{DNAseq} here~~
}
  \item{rev}{
%%     ~~Describe \code{rev} here~~
}
  \item{compl}{
%%     ~~Describe \code{compl} here~~
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
function (DNAseq, rev = TRUE, compl = TRUE) 
{
    if (is.character(DNAseq)) {
        resultVect <- sapply(DNAseq, (function(seqString) {
            mySeq <- toupper(seqString)
            if (rev == TRUE) {
                mySeq <- paste(sapply(1:nchar(mySeq), (function(i) {
                  pos <- nchar(mySeq) + 1 - i
                  substr(mySeq, pos, pos)
                })), collapse = "")
            }
            if (compl == TRUE) {
                mySeq <- paste(sapply(1:nchar(mySeq), (function(i) {
                  aBase <- substr(mySeq, i, i)
                  complBase <- c("A", "T", "C", "G", "N", "A")
                  names(complBase) <- c("T", "A", "G", "C", "N", 
                    "U")
                  returnBase <- complBase[aBase]
                  if (is.na(returnBase)) {
                    stop("Bad input")
                  }
                  returnBase
                })), collapse = "")
            }
            return(mySeq)
        }))
        names(resultVect) <- NULL
        return(resultVect)
    }
    else {
        warning("Bad input")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
