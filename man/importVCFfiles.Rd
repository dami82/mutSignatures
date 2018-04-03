\name{importVCFfiles}
\alias{importVCFfiles}
\title{
Import VCF Files
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
importVCFfiles(vcfFiles, sampleNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcfFiles}{
%%     ~~Describe \code{vcfFiles} here~~
}
  \item{sampleNames}{
%%     ~~Describe \code{sampleNames} here~~
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
function (vcfFiles, sampleNames = NULL) 
{
    my.colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
        "FILTER", "INFO", "FORMAT", "XTR1", "XTR2", "XTR3")
    if (is.null(sampleNames) | length(sampleNames) != length(vcfFiles)) {
        bypassNames <- paste("sample", 1:length(vcfFiles), sep = ".")
        sampleNames <- vcfFiles
        sampleNames <- sub("\\.vcf$", "", sub("^.*(\\\\|/)", 
            "", tolower(sampleNames)))
        sampleNames[sampleNames == ""] <- bypassNames[sampleNames == 
            ""]
    }
    out <- sapply(1:length(vcfFiles), (function(j) {
        x <- vcfFiles[j]
        if (!is.file(x)) {
            NULL
        }
        else {
            tmpVCF <- read.delim(x, comment.char = "#", header = F, 
                stringsAsFactors = F)
            for (i in 1:ncol(tmpVCF)) {
                colnames(tmpVCF)[i] <- my.colnames[i]
            }
            tmpVCF$SAMPLEID <- sampleNames[j]
            tmpVCF
        }
    }), simplify = FALSE, USE.NAMES = TRUE)
    out <- do.call(rbind, out)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
