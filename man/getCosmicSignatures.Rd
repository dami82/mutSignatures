\name{getCosmicSignatures}
\alias{getCosmicSignatures}
\title{
Retrieve COSMIC Mutational Signatures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCosmicSignatures(forceUseMirror = FALSE, asMutSign = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forceUseMirror}{
%%     ~~Describe \code{forceUseMirror} here~~
}
  \item{asMutSign}{
%%     ~~Describe \code{asMutSign} here~~
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
function (forceUseMirror = FALSE, asMutSign = TRUE) 
{
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
    cosmic.url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
    my_fullW <- NULL
    if (!forceUseMirror) {
        my_fullW <- tryCatch({
            TMP <- suppressWarnings(read.delim(cosmic.url, header = TRUE))
            rownames(TMP) <- TMP$Somatic.Mutation.Type
            TMP <- TMP[, grep("Signature", colnames(TMP))]
            TMP
        }, error = function(e) NULL)
    }
    if (is.null(my_fullW)) {
        private.mirror.url <- "http://www.labwizards.com/rlib/mutSignatures/cosmic.signatures.csv"
        my_fullW <- tryCatch({
            suppressWarnings(read.csv(private.mirror.url, header = TRUE, 
                as.is = TRUE, row.names = 1))
        }, error = function(e2) {
            NULL
        })
    }
    if (is.null(my_fullW)) {
        message("An error occurred!")
        return(NULL)
    }
    else {
        if (sum(mutType.labels \%in\% rownames(my_fullW)) == length(mutType.labels)) {
            obj2rt <- my_fullW[mutType.labels, ]
            if (asMutSign) 
                obj2rt <- as.mutation.signatures(obj2rt)
            return(obj2rt)
        }
        else {
            message("An error occurred!")
            return(NULL)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
