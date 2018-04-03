\name{attachMutType}
\alias{attachMutType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Attach Mutation Type
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
attachMutType(mutData, ref_colName = "reference_allele", var_colName = "variant_allele", var2_colName = NULL, context_colName = "context", format = 1, mutType_dict = "alexa", mutType_colName = "mutType")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutData}{
%%     ~~Describe \code{mutData} here~~
}
  \item{ref_colName}{
%%     ~~Describe \code{ref_colName} here~~
}
  \item{var_colName}{
%%     ~~Describe \code{var_colName} here~~
}
  \item{var2_colName}{
%%     ~~Describe \code{var2_colName} here~~
}
  \item{context_colName}{
%%     ~~Describe \code{context_colName} here~~
}
  \item{format}{
%%     ~~Describe \code{format} here~~
}
  \item{mutType_dict}{
%%     ~~Describe \code{mutType_dict} here~~
}
  \item{mutType_colName}{
%%     ~~Describe \code{mutType_colName} here~~
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
function (mutData, ref_colName = "reference_allele", var_colName = "variant_allele", 
    var2_colName = NULL, context_colName = "context", format = 1, 
    mutType_dict = "alexa", mutType_colName = "mutType") 
{
    if (!((is.data.frame(mutData) | is.matrix(mutData)) & sum(c(ref_colName, 
        var_colName, context_colName) \%in\% colnames(mutData)) == 
        3)) 
        stop("Issue with the input dataset. Make sure to feed in a data.frame or\n          a matrix and double check the name of the fields pointing to chromosome\n          name, start and end positions")
    if (!(format \%in\% c(1, 2))) 
        stop("Please, specify a valid format number (example: 1)")
    if (!(is.null(var2_colName))) {
        if (!(var2_colName \%in\% colnames(mutData))) 
            stop("Invalid var2 column")
    }
    if (!is.character(mutType_colName) | length(mutType_colName) > 
        1) 
        stop("Bad mutType_colName")
    if (mutType_colName \%in\% colnames(mutData)) 
        stop("mutType_colName already exists as column name in the current dataset")
    mutData <- data.frame(mutData, stringsAsFactors = FALSE, 
        row.names = NULL)
    my.key.cols <- c(ref_colName, var_colName, var2_colName, 
        context_colName)
    my.key.cols <- my.key.cols[!is.na(my.key.cols)]
    for (clmn in my.key.cols) {
        mutData[, clmn] <- as.character(as.vector(mutData[, clmn]))
    }
    message("Assigning mutation types ", appendLF = FALSE)
    mutData[, mutType_colName] <- sapply(1:nrow(mutData), (function(i) {
        if (nrow(mutData) > 1000 & i \%in\% as.integer(seq(1, nrow(mutData), 
            length.out = 20))) 
            message(".", appendLF = FALSE)
        ctx.len <- nchar(mutData[i, context_colName])
        half.ln <- (ctx.len - 1)/2
        mid.seq <- substr(mutData[i, context_colName], (half.ln + 
            1), (half.ln + 1))
        pre.seq <- substr(mutData[i, context_colName], 1, half.ln)
        post.seq <- substr(mutData[i, context_colName], (half.ln + 
            2), ctx.len)
        if (mid.seq != mutData[i, ref_colName] | (is.null(var2_colName) & 
            mid.seq == mutData[i, var_colName]) | (tryCatch({
            (mid.seq == mutData[i, var_colName] & mid.seq == 
                mutData[i, var2_colName])
        }, error = function(e) {
            FALSE
        }))) {
            mut.base <- NA
        }
        else {
            if (mutData[i, ref_colName] != mutData[i, var_colName]) {
                mut.base <- mutData[i, var_colName]
            }
            else if (!is.null(var2_colName)) {
                if (mutData[i, ref_colName] != mutData[i, var2_colName]) {
                  mut.base <- mutData[i, var2_colName]
                }
                else {
                  mut.base <- NA
                }
            }
            else {
                mut.base <- NA
            }
            if (is.na(mut.base)) {
                NA
            }
            else {
                paste(mid.seq, ".", mut.base, "[", pre.seq, mid.seq, 
                  post.seq, "][", pre.seq, mut.base, post.seq, 
                  "]", sep = "", collapse = "")
            }
        }
    }))
    message(". Done!", appendLF = TRUE)
    if (sum(is.na(mutData[, mutType_colName])) > 0) {
        message(paste("Removing", sum(is.na(mutData[, mutType_colName])), 
            "positions."))
        mutData <- mutData[!is.na(mutData[, mutType_colName]), 
            ]
    }
    message("Now applying RevCompl transformation", appendLF = FALSE)
    if (mutType_dict == "alexa") {
        idx <- grep("^((G|A)\\.)", mutData[, mutType_colName])
        mutData[idx, mutType_colName] <- sapply(mutData[idx, 
            mutType_colName], (function(seq) {
            base.wt <- revCompl(gsub("\\..+$", "", seq))
            base.mut <- revCompl(gsub("^.+\\.", "", gsub("\\[.+$", 
                "", seq)))
            seq.wt <- revCompl(gsub("^.+\\[", "", gsub("\\]\\[.+$", 
                "", seq)))
            seq.mut <- revCompl(gsub("^.+\\]\\[", "", gsub("\\]$", 
                "", seq)))
            paste(base.wt, ".", base.mut, "[", seq.wt, "][", 
                seq.mut, "]", sep = "", collapse = "")
        }))
    }
    else if (mutType_dict == "custom") {
        idx <- grep("^((G|T)\\.)", mutData[, mutType_colName])
        mutData[idx, mutType_colName] <- sapply(mutData[idx, 
            mutType_colName], (function(seq) {
            base.wt <- revCompl(gsub("\\..+$", "", seq))
            base.mut <- revCompl(gsub("^.+\\.", "", gsub("\\[.+$", 
                "", seq)))
            seq.wt <- revCompl(gsub("^.+\\[", "", gsub("\\]\\[.+$", 
                "", seq)))
            seq.mut <- revCompl(gsub("^.+\\]\\[", "", gsub("\\]$", 
                "", seq)))
            paste(base.wt, ".", base.mut, "[", seq.wt, "][", 
                seq.mut, "]", sep = "", collapse = "")
        }))
    }
    message(". Done!", appendLF = TRUE)
    message("Final formatting", appendLF = FALSE)
    mutData[, mutType_colName] <- sapply(mutData[, mutType_colName], 
        (function(seq) {
            base.wt <- gsub("\\..+$", "", seq)
            base.mut <- gsub("^.+\\.", "", gsub("\\[.+$", "", 
                seq))
            seq.wt <- gsub("^.+\\[", "", gsub("\\]\\[.+$", "", 
                seq))
            seq.mut <- gsub("^.+\\]\\[", "", gsub("\\]$", "", 
                seq))
            half.len <- (nchar(seq.wt) - 1)/2
            pre.seq <- substr(seq.wt, 1, half.len)
            post.seq <- substr(seq.wt, half.len + 2, nchar(seq))
            if (format == 1) {
                paste(pre.seq, "[", base.wt, ">", base.mut, "]", 
                  post.seq, sep = "", collapse = "")
            }
            else if (format == 2) {
                paste(pre.seq, post.seq, ".", base.wt, ">", base.mut, 
                  sep = "", collapse = "")
            }
            else {
                paste(base.wt, ".", base.mut, "[", seq.wt, "][", 
                  seq.mut, "]", sep = "", collapse = "")
            }
        }))
    message(". Done!", appendLF = TRUE)
    return(mutData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
