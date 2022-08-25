#
##  ~~~~~~~~~~~~~~~~~~~~~~~~~
### ~~~~~ MutSignatures ~~~~~
##  ...all methods...
#   ~~~~~~~~~~~~~~~~~~~~~~~~~

# by Damiano Fantini, Ph.D.
# 2020-Mar-05
#
# This file includes class methods, as well
# as series of generics.
# Entries are ordered by 
# method types (getters, setters, casters, ...) 
# and then by method/generic. 


#
# 1. Getters
#

#' Method getFwkParam.
#' 
#' Retrieve the list of parameters used for running a Mutation Signature Analysis.
#' 
#' @param x a mutFrameworkParams-class object
#' @param label string, corresponding to the parameter name to extract
#' 
#' @rdname getFwkParam
#' @exportMethod getFwkParam
setGeneric("getFwkParam", function(x, label) {
  standardGeneric("getFwkParam")
})

#' @rdname getFwkParam
#' @aliases getFwkParam,mutFrameworkParams,character-method
setMethod("getFwkParam", signature(x = "mutFrameworkParams", label = "character"),
          function(x, label) {
            out <- x@params
            out[[label]]})

## ~~~~~~~~~~~~~~~~~~


#' Method getMutationTypes.
#' 
#' Retrieve the list of mutation types from an object.
#' 
#' @param x an object to extract Mutation types from, i.e. a mutationSignatures-class
#' or a mutationCounts-class object
#' 
#' @rdname getMutationTypes
#' @exportMethod getMutationTypes
setGeneric("getMutationTypes", function(x) {
  standardGeneric("getMutationTypes")
})

#' @rdname getMutationTypes
#' @aliases getMutationTypes,mutationSignatures-method
setMethod("getMutationTypes", "mutationSignatures",
          function(x){
            base::as.vector(x@mutTypes[,1])
          })

#' @rdname getMutationTypes
#' @aliases getMutationTypes,mutationCounts-method
setMethod("getMutationTypes", "mutationCounts",
          function(x){
            base::as.vector(x@mutTypes[,1])
          })

## ~~~~~~~~~~~~~~~~~~


#' Method getSampleIdentifiers.
#' 
#' Retrieve the list of sample identifiers from an object.
#' 
#' @param x an object to extract Mutation types from, i.e. a mutationCounts-class
#' or a mutSignExposures-class object
#' 
#' @rdname getSampleIdentifiers
#' @exportMethod getSampleIdentifiers
setGeneric("getSampleIdentifiers", function(x) {
  standardGeneric("getSampleIdentifiers")
})

#' @rdname getSampleIdentifiers
#' @aliases getSampleIdentifiers,mutationCounts-method
setMethod("getSampleIdentifiers", "mutationCounts",
          function(x = "mutationCounts"){
            base::as.vector(x@sampleId[,1])
          })

#' @rdname getSampleIdentifiers
#' @aliases getSampleIdentifiers,mutSignExposures-method
setMethod("getSampleIdentifiers", "mutSignExposures",
          function(x){
            base::as.vector(x@sampleId[,1])
          })

## ~~~~~~~~~~~~~~~~~~


#' Method getCounts.
#' 
#' Retrieve mutation counts from an object.
#' 
#' @param x an object to extract Mutation counts from, i.e. a mutationCounts-class object
#' 
#' @rdname getCounts
#' @exportMethod getCounts  
setGeneric("getCounts", function(x) {
  standardGeneric("getCounts")
})

#' @rdname getCounts
#' @aliases getCounts,mutationCounts-method
setMethod("getCounts", signature(x="mutationCounts"),
          function(x) {
            
            out <- x@counts
            rownames(out) <- x@mutTypes[,1]
            colnames(out) <- x@sampleId[,1]
            out
          })

## ~~~~~~~~~~~~~~~~~~


#' Method getSignatureIdentifiers.
#' 
#' Retrieve the list of signature identifiers from an object.
#' 
#' @param x an object to extract Signature Identifiers from, i.e. a mutSignExposures-class
#' or a mutationSignatures-class object
#' 
#' @rdname getSignatureIdentifiers
#' @exportMethod getSignatureIdentifiers
setGeneric("getSignatureIdentifiers", function(x) {
  standardGeneric("getSignatureIdentifiers")
})

#' @rdname getSignatureIdentifiers
#' @aliases getSignatureIdentifiers,mutSignExposures-method
setMethod("getSignatureIdentifiers", signature(x="mutSignExposures"),
          function(x){
            base::as.vector(x@signatureId[,1])
          })

#' @rdname getSignatureIdentifiers
#' @aliases getSignatureIdentifiers,mutationSignatures-method
setMethod("getSignatureIdentifiers", signature(x="mutationSignatures"),
          function(x){
            base::as.vector(x@signatureId[,1])
          })

## ~~~~~~~~~~~~~~~~~~


#
# 2. Setters
#

#' Method setFwkParam.
#' 
#' Set or update one of the parameters in a mutFrameworkParams-class object. 
#' Individual paramaters can be set or updated, by passing the parameter label, and the corresponding
#' parameter value.
#' 
#' @param x an object to extract Signature Identifiers from, i.e. a mutSignExposures-class
#' @param label string corresponding to the parameter label to be updated
#' @param value new value (string or numeric) of the parameter to be updated
#' 
#' @rdname setFwkParam
#' @exportMethod setFwkParam
setGeneric("setFwkParam", function(x, label, value) {
  standardGeneric("setFwkParam")
})

#' @rdname setFwkParam
#' @aliases setFwkParam,mutFrameworkParams,character,ANY-method
setMethod("setFwkParam", signature(x="mutFrameworkParams", 
                                   label = "character", 
                                   value = "ANY"),
          function(x, label, value) {
            if (length(label) != 1 | length(value) != 1)
              stop("Bad input")
            if(!is.character(label))
              stop("Bad input")
            
            x@params[[ label ]] <- value
            x
          })

## ~~~~~~~~~~~~~~~~~~


#' Method setSignatureNames.
#' 
#' Update signature names of a Mutation Signatures object.
#' 
#' @param x a mutationSignatures object
#' @param names character vector, these are the new names that will be assigned 
#' to the signatures.
#' 
#' @rdname setSignatureNames
#' 
#' @exportMethod setSignatureNames
setGeneric("setSignatureNames", function(x, names) {
  standardGeneric("setSignatureNames")
})

#' @rdname setSignatureNames
#' 
#' @aliases setSignatureNames,mutationSignatures,character-method
setMethod("setSignatureNames", 
          signature(x="mutationSignatures", names = "character"),
          function(x, names) {
            
            old_sigs <- as.character(x@signatureId$ID)
            
            # Checks
            if (is.null(names) || length(names) < 1)
              stop("No names were provided")
            
            if (sum(duplicated(names)) > 0 || sum(is.na(names)) > 0)
              stop('N/A or duplicated signature names can NOT be accepted')
            
            if (length(names) != length(old_sigs))
              stop(paste0("Wrong number of signature names was provided. ", 
                          "Please provide n=", length(old_sigs), 
                          " signature names"))
            
            x@signatureId$ID <- names
            x
          })

## ~~~~~~~~~~~~~~~~~~



#
# 3. Coercers
#

#' Method coerceObj.
#' 
#' Cast an object to a different format, by extracting and returning the
#' most appropriate information. Note that data.frames can be coerced to one
#' of the classes defined in the mutSignatures package using coerceObj.
#' 
#' @param x an object to coerce to a different format
#' @param to string, indicates the expected format (such as list or data.frame)
#' @param ... additional parameters passed to the functions used for the coercion
#' 
#' @rdname coerceObj-methods
#' @exportMethod coerceObj
setGeneric("coerceObj", function(x, to, ...) {
  standardGeneric("coerceObj")
})

#' @rdname coerceObj-methods
#' @aliases coerceObj,mutFrameworkParams,character-method
setMethod("coerceObj", signature(x = "mutFrameworkParams", to = "character"),
          function(x, to) {
            if (to[1] == "list") {
              x@params
            } else if (to[1] == class(x)[1]){
              x
            } else {
              stop("Object could not be coerced to the desired data type")
            }
          })


#' @rdname coerceObj-methods
#' @aliases coerceObj,mutationSignatures,character-method
setMethod("coerceObj", signature(x = "mutationSignatures", to = "character"),
          function(x, to) {
            if (to[1] == "data.frame") {
              out <- x@mutationFreq
              rownames(out) <- x@mutTypes[,1]
              colnames(out) <- x@signatureId[,1]
              out
            } else if (to[1] == "matrix") {
              out <- x@mutationFreq
              rownames(out) <- x@mutTypes[,1]
              colnames(out) <- x@signatureId[,1]
              base::as.matrix(out)
            } else if (to[1] == "list") {
              list(mutationFreq=x@mutationFreq,
                   mutTypes=x@mutTypes,
                   signatureId=x@signatureId)
            } else if (to[1] == class(x)[1]){
              x
            } else {
              stop("Object could not be coerced to the desired data type")
            }         
          })

#' @rdname coerceObj-methods
#' @aliases coerceObj,mutationCounts,character-method
setMethod("coerceObj", signature(x = "mutationCounts", to = "character"),
          function(x, to, ...) {
            
            mDots <- list(...)
            if ("keepNames" %in% names(mDots)) {
              keepNames <- tryCatch({as.logical(mDots[["keepNames"]])[1]}, error = function(e) {FALSE})
            } else {
              keepNames <- FALSE
            }
            
            
            if (to[1] == "data.frame") {
              out <- getCounts(x)
            } else if (to[1] == "matrix") {
              out <- getCounts(x)
              out <- base::as.matrix(out)
              if (is.logical(keepNames[1])) {
                if (!keepNames[1]) {
                  dimnames(out) <- NULL
                }
              }
              out
            } else if (to[1] == class(x)[1]){
              x
            } else {
              stop("Object could not be coerced to the desired data type")
            }
          })

#' @rdname coerceObj-methods
#' @aliases coerceObj,mutSignExposures,character-method
setMethod("coerceObj", signature(x = "mutSignExposures", to = "character"),
          function(x, to, ...) {
            
            mDots <- list(...)
            if ("transpose" %in% names(mDots)) {
              transpose <- tryCatch({as.logical(mDots[["transpose"]])[1]}, error = function(e) {FALSE})
            } else {
              transpose <- FALSE
            }
            
            if (to[1] == "data.frame" & is.logical(transpose) & transpose[1]) {
              out <- base::data.frame(t(x@exposures))
              rownames(out) <- x@sampleId[,1]
              colnames(out) <- x@signatureId[,1]
              out
            } else if (to[1] == "matrix" & is.logical(transpose) & transpose[1]) {
              out <- t(x@exposures)
              rownames(out) <- x@sampleId[,1]
              colnames(out) <- x@signatureId[,1]
              base::as.matrix(out)
            } else if (to[1] == "data.frame" & is.logical(transpose) & !transpose[1]) {
              out <- x@exposures
              rownames(out) <- x@signatureId[,1]
              colnames(out) <- x@sampleId[,1]
              out
            } else if (to[1] == "matrix" & is.logical(transpose) & !transpose[1]) {
              out <- x@exposures
              rownames(out) <- x@signatureId[,1]
              colnames(out) <- x@sampleId[,1]
              base::as.matrix(out)
            } else if (to[1] == class(x)[1]){
              x
            } else {
              stop("Object could not be coerced to the desired data type")
            }         
          })

#' @rdname coerceObj-methods
#' @aliases coerceObj,data.frame,character-method
setMethod("coerceObj", signature(x = "data.frame", to = "character"),
          function(x, to, ...) {
            if (to[1] == "mutationCounts" & nargs() %in% 2:4) {
              as.mutation.counts(x, ...)
            } else if (to[1] == "mutationSignatures" & nargs() == 2){
              as.mutation.signatures(x)
            } else if (to[1] == "mutSignExposures" & nargs() %in% 2:3) {
              as.mutsign.exposures(x, ...)
            } else if (to[1] == class(x)[1]){
              x
            } else {
              stop("Object could not be coerced to the desired data type")
            }
          })

## ~~~~~~~~~~~~~~~~~~
# as.list(), as.data.frame() and as.matrix()- methods


#' Convert a mutFrameworkParams object to list.
#' 
#' Coerce a mutFrameworkParams-class object to list by applying the coerceObj method. 
#'
#' @docType methods
#' @param x a mutFrameworkParams object
#' 
#' @aliases as.list,mutFrameworkParams-method
#' @export
setMethod("as.list", signature(x = "mutFrameworkParams"),
          function(x) {
            coerceObj(x, to = "list")
          })

#' Convert a mutationSignatures object to list.
#' 
#' Coerce a mutationSignatures-class object to list by applying the coerceObj method. 
#'
#' @docType methods
#' @param x a mutationSignatures object
#' 
#' @aliases as.list,mutationSignatures-method
#' @export
setMethod("as.list", "mutationSignatures", 
          function(x) {
            coerceObj(x, to = "list")
          })

#' Convert a mutationSignatures object to data.frame.
#' 
#' Coerce a mutationSignatures-class object to data.frame by applying the coerceObj method. 
#'
#' @docType methods
#' @param x a mutationSignatures object
#' 
#' @aliases as.data.frame,mutationSignatures-method
#' @export
setMethod("as.data.frame", signature(x ="mutationSignatures"),
          function(x){
            coerceObj(x, to = "data.frame")
          })


#' Convert a mutationCounts object to data.frame.
#' 
#' Coerce a mutationCounts-class object to data.frame by applying the coerceObj method. 
#'
#' @docType methods
#' @param x a mutationCounts object
#' 
#' @aliases as.data.frame,mutationCounts-method
#' @export
setMethod("as.data.frame", signature(x ="mutationCounts"),
          function(x){
            coerceObj(x, to = "data.frame")
          })


#' Convert and/or transpose a mutSignExposures object to data.frame.
#' 
#' Coerce a mutSignExposures-class object to data.frame by applying the coerceObj method.
#' The data.frame can be returned in a transposed or non-transposed format. 
#'
#' @docType methods
#' 
#' @param x a mutSignExposures object
#' @param row.names NULL, not used
#' @param optional NULL, not used
#' @param ... additional parameters to be passed to coerceObj, such as transpose (logical)
#' 
#' @aliases as.data.frame,mutSignExposures-method
#' @export
setMethod("as.data.frame", signature(x ="mutSignExposures"),
          function(x, row.names = NULL, optional = NULL, ...){
            
            mDots <- list(...)
            if ("transpose" %in% names(mDots)) {
              transpose <- tryCatch({as.logical(mDots[["transpose"]])[1]}, error = function(e) {FALSE})
            } else {
              transpose <- FALSE
            }
            
            if (is.logical(transpose)) {
              coerceObj(x, to = "data.frame", transpose = transpose[1])  
            } else {
              coerceObj(x, to = "data.frame")
            }
          })

#' Convert a mutationCounts object to matrix.
#' 
#' Coerce a mutationCounts-class object to matrix by applying the coerceObj method. 
#'
#' @docType methods
#' @param x a mutationCounts object
#' 
#' @aliases as.matrix,mutationCounts-method
#' @export
setMethod("as.matrix", signature(x = "mutationCounts"),
          function(x) {
            coerceObj(x, to = "matrix", keepNames = FALSE)
          })

## ~~~~~~~~~~~~~~~~~~



#' Method as.mutation.counts.
#' 
#' Cast a data.frame into a mutationCounts-class object.
#' 
#' @param x an object to extract Signature Identifiers from, i.e. a mutSignExposures-class
#' @param rownames character vector to overwrite data row names. Can be NULL if rownames(x) is not NULL.
#' @param colnames character vector to overwrite data column names. Can be NULL if colnames(x) is not NULL.
#' 
#' @rdname as.mutation.counts
#' 
#' @exportMethod as.mutation.counts
setGeneric("as.mutation.counts", function(x, rownames=NULL, colnames=NULL) {
  standardGeneric("as.mutation.counts")
})

#' @rdname as.mutation.counts
#' @importFrom methods new
#' 
#' @aliases as.mutation.counts,data.frame,ANY,ANY-method
setMethod("as.mutation.counts", signature(x="data.frame",rownames = "ANY", colnames = "ANY"),
          function(x, rownames = NULL, colnames = NULL) {
            
            # Checks
            if (is.null(rownames) & (is.null(rownames(x))))
              stop("Incomplete Mutation Type info")
            
            if (is.null(colnames) & is.null(colnames(x)))
              stop("Incomplete Mutation Type info")
            
            if (!is.null(rownames)){
              if (!is.character(rownames) |
                  length(rownames) != nrow(x) |
                  sum(duplicated(rownames) > 0))
                stop("Unsuitable Mutation Types")
            }
            
            if (!is.null(colnames)){
              if (!is.character(colnames) |
                  length(colnames) != ncol(x) |
                  sum(duplicated(colnames) > 0))
                stop("Unsuitable Sample Identifiers")
            }
            
            if (is.null(rownames))
              rownames <- rownames(x)
            
            if (is.null(colnames))
              colnames <- colnames(x)
            
            fOUT <- new(Class = "mutationCounts", 
                        x = x, 
                        muts = data.frame(mutTypes = rownames, stringsAsFactors = FALSE), 
                        samples = data.frame(ID = colnames, stringsAsFactors = FALSE))            
            
            fOUT
          })



## ~~~~~~~~~~~~~~~~~~

#' Method as.mutation.signatures.
#' 
#' Cast a data.frame into a mutationSignatures-class object.
#' 
#' @param x a data.frame to be converted to a mutationSignatures-class object.
#' 
#' @rdname as.mutation.signatures
#' @exportMethod as.mutation.signatures
setGeneric("as.mutation.signatures", function(x) {
  standardGeneric("as.mutation.signatures")
})

#' @rdname as.mutation.signatures
#' @importFrom methods new
#' @aliases as.mutation.signatures,data.frame-method
setMethod("as.mutation.signatures", signature(x ="data.frame"),
          function(x){
            
            data <- x
            rw1 <- data.frame(type = rownames(x), stringsAsFactors = FALSE)
            cl1 <- data.frame(ID = colnames(x), stringsAsFactors = FALSE)
            
            out <- new(Class = "mutationSignatures", 
                       x= data, 
                       muts=rw1, 
                       signNames=cl1)
            out
          })

## ~~~~~~~~~~~~~~~~~~


#' Method as.mutsign.exposures.
#' 
#' Cast a data.frame into a mutSignExposures-class object.
#' 
#' @param x a data.frame to be converted to a mutSignExposures-class object.
#' @param samplesAsCols logical, are samples listed as columns in the input data.frame. If FALSE, 
#' samples are expected to be listed as rows in the input data.frame 
#' 
#' @rdname as.mutsign.exposures
#' @exportMethod as.mutsign.exposures
setGeneric("as.mutsign.exposures", function(x, samplesAsCols = TRUE) {
  standardGeneric("as.mutsign.exposures")
})

#' @rdname as.mutsign.exposures
#' 
#' @importFrom methods new
#' @aliases as.mutsign.exposures,data.frame,logical-method
setMethod("as.mutsign.exposures", signature(x ="data.frame", samplesAsCols = "logical"),
          function(x, samplesAsCols = TRUE){
            if (is.null(rownames(x)) | is.null(colnames(x))) 
              stop("Bad input")
            
            data <- x
            rw1 <- data.frame(type = rownames(x), stringsAsFactors = FALSE)
            cl1 <- data.frame(ID = colnames(x), stringsAsFactors = FALSE)
            
            if(samplesAsCols) {
              out <- new(Class = "mutSignExposures", 
                         x = data, 
                         samples=cl1, 
                         signNames=rw1)
            } else {
              out <- new(Class = "mutSignExposures", 
                         x = base::as.data.frame(base::t(data)), 
                         samples = rw1, 
                         signNames = cl1)
              
            }
            out
          })

## ~~~~~~~~~~~~~~~~~~

#
# 4. Combine and Subset
#

#' Subset a mutationSignatures-class object.
#'
#' @aliases [,mutationSignatures,numeric-method
#' @docType methods
#' 
#' @param x a mutationSignatures-class object to subset
#' @param i numeric, indeces of the elements to be extracted
#' 
#' @importFrom methods new
#' @export 
setMethod("[", signature(x = "mutationSignatures", i = "numeric"),
          function (x, i){
            
            if(nargs() > 2 )
              stop("object of type 'S4' is not subsettable using more than 1 dimension")
            
            new(Class = "mutationSignatures", 
                x=base::as.data.frame(x@mutationFreq[,i]), 
                muts = x@mutTypes,
                signNames=base::data.frame(x@signatureId[i,], 
                                           stringsAsFactors = FALSE,
                                           row.names = NULL))
          })

#' Subset a mutationCounts-class object.
#'
#' @aliases [,mutationCounts,numeric-method
#' @docType methods
#' 
#' @param x a mutationCounts-class object to subset
#' @param i numeric, indeces of the elements to be extracted
#'
#' @importFrom methods new 
#' @export
setMethod("[", signature(x = "mutationCounts", i = "numeric"),
          function (x, i){
            
            if(nargs() > 2 )
              stop("object of type 'S4' is not subsettable using more than 1 dimension")
            
            new(Class = "mutationCounts", 
                x=base::as.data.frame(x@counts[,i]), 
                muts = x@mutTypes,
                samples=base::data.frame(x@sampleId[i,], 
                                         stringsAsFactors = FALSE,
                                         row.names = NULL))
          })

#' Subset a mutSignExposures-class object.
#'
#' @aliases [,mutSignExposures,numeric-method
#' @docType methods
#' 
#' @param x a mutSignExposures-class object to subset
#' @param i numeric, indeces of the elements to be extracted
#'
#' @importFrom methods new
#' @export 
setMethod("[", signature(x = "mutSignExposures", i = "numeric"),
          function (x, i){
            
            if(nargs() > 2 )
              stop("object of type 'S4' is not subsettable using more than 1 dimension")
            
            new(Class = "mutSignExposures", 
                x = base::as.data.frame(x@exposures[,i]), 
                signNames = x@signatureId,
                samples=base::data.frame(x@sampleId[i,], 
                                         stringsAsFactors = FALSE, 
                                         row.names = NULL))
          })


#' Combine two mutationSignatures-class objects.
#'
#' @name cbind2
#' @aliases cbind2,mutationSignatures,mutationSignatures-method
#' @docType methods
#' 
#' @param x the first mutSignExposures-class object to combine
#' @param y the first mutSignExposures-class object to combine
#' 
#' @details a variant of this method accepting more than 2 object to combine together is 
#' under preparation and be available soon...
#' 
#' @importFrom methods new cbind2
#' @export
setMethod("cbind2", signature(x="mutationSignatures", y="mutationSignatures"), 
          
          function(x, y) {         #  <--------------------omitting `...`
            
            xMut <- getMutationTypes(x)
            yMut <- getMutationTypes(y)
            
            if (length(xMut) == length(yMut) &
                sum(!xMut %in% yMut) == 0) {
              
              orderY <- sapply(xMut, function(ix) which(yMut == ix))
              #sortedYmut <- yMut[orderY]
              sortedYexp <- y@mutationFreq[orderY, ]
              
              new(Class = "mutationSignatures", 
                  x=base::as.data.frame(cbind(x@mutationFreq, sortedYexp)), 
                  muts = x@mutTypes,
                  signNames=data.frame(c(x@signatureId[,1], y@signatureId[,1]) , 
                                       stringsAsFactors = FALSE))
              
            } else {
              stop("Incompatible Objects")
            }
          })

## ~~~~~~~~~~~~~~~~~~

#
# 5. Plot 
#

#' Method msigPlot.
#' 
#' Generate standard plots using data from mutsignature-class objects.
#' 
#' @param x a mutSignatures object
#' @param ... additional parameters, including standard graphical parameters, as well as 
#' a set of class-specific arguments:
#' * x is a mutationSignatures object: 
#'     - signature, numeric; numeric index of the signature to display
#'     - main, string; title of the plot
#' * x is a mutationCounts object
#'     - sample, numeric or string, i.e. the identifier or the index of the sample to be plotted
#'     - main, string, title of the plot
#'  * x is a mutSignExposures object
#'     - top, integer, the maximum number of samples to be plotted
#'
#' 
#' @rdname msigPlot
#' @exportMethod msigPlot
setGeneric("msigPlot", function(x, ...) {
  standardGeneric("msigPlot")
})



#' @rdname msigPlot
#' @aliases msigPlot,mutationSignatures-method
setMethod("msigPlot", signature(x = "mutationSignatures"), 
          function(x, ...) {
            
            # Get Arguments and parse out what we want
            all.args <- list(...)
            
            if (length(all.args) < 1) {
              signature <- 1
            } else if ("signature" %in% names(all.args)) {
              signature <- all.args[["signature"]]
            } else {
              signature <- 1
            }
            
            if (length(all.args) < 1) {
              main <- NULL
            } else if ("main" %in% names(all.args)) {
              main <- all.args[["main"]]
            } else {
              main <- NULL
            }
            
            # Get all other plot args
            # Defaulting
            freq <- TRUE; ylim <- "auto"; ylab <- "Fraction of Variants";
            xlab <- "Sequence Motifs"; xaxis_cex <- 0.475;
            cols <- c("#4eb3d3", "#040404", "#b30000", "#bdbdbd", "#41ab5d", "#dd3497")
            
            if ("freq" %in% names(all.args)) {
              freq <- all.args[["freq"]]
            }
            if ("ylim" %in% names(all.args)) {
              ylim <- all.args[["ylim"]]
            }
            if ("ylab" %in% names(all.args)) {
              ylab <- all.args[["ylab"]]
            }
            if ("xlab" %in% names(all.args)) {
              xlab <- all.args[["xlab"]]
            }
            if ("xaxis_cex" %in% names(all.args)) {
              xaxis_cex <- all.args[["xaxis_cex"]]
            }
            if ("cols" %in% names(all.args)) {
              cols <- all.args[["cols"]]
            }
            
            mutCounts <- x
            objOut <- NULL
            if (is.numeric(signature)) {
              if(length(signature) == 1 & signature[1] <= length(getSignatureIdentifiers(mutCounts))) {
                objOut <- coerceObj(x = mutCounts, to = "data.frame")[,signature]
                label <-  getSignatureIdentifiers(mutCounts)[signature]
              }
            } else if (is.character(signature)) {
              if (length(signature) == 1 & signature[1] %in% getSignatureIdentifiers(mutCounts)) {
                objOut <- coerceObj(x = mutCounts, to = "data.frame")[,signature]
                label <- signature
              }
            } 
            
            if(!is.null(objOut)){
              final.label <- ifelse(!is.null(main), main, label)
              plotMutTypeProfile(mutCounts = objOut, 
                                 mutLabs = getMutationTypes(mutCounts),
                                 main = final.label, freq = freq, 
                                 ylim = ylim, ylab = ylab, xlab = xlab, 
                                 xaxis_cex = xaxis_cex, cols = cols)
              
            } else {
              stop("Signature not found")
            }
          })
#---

#' @rdname msigPlot
#' @aliases msigPlot,mutationCounts-method
setMethod("msigPlot", signature(x = "mutationCounts"), 
          function(x, ...) {
            mutCounts <- x
            objOut <- NULL
            
            # Get Arguments and parse out what we want
            all.args <- list(...)
            new.args <- list(freq = TRUE, ylim = "auto", 
                             ylab = "Fraction of Variants", 
                             xlab = "Sequence Motifs", 
                             xaxis_cex = 0.475,
                             cols = c("#4eb3d3","#040404", "#b30000", "#bdbdbd", "#41ab5d", "#dd3497"),
                             main = "MutType Profile", sample = 1)
            
            if (length(all.args) > 0) {
              for (ArG in names(new.args)) {
                if (ArG %in% names(all.args)) {
                  new.args[[ArG]] <- all.args[[ArG]]  
                }
              }
            }
            
            if (is.numeric(new.args$sample)) {
              if(length(new.args$sample) == 1 & 
                 new.args$sample[1] <= length(getSampleIdentifiers(mutCounts))) {
                objOut <- coerceObj(x = mutCounts, to = "data.frame")[,new.args$sample]
                label <-  getSampleIdentifiers(mutCounts)[new.args$sample]
              }
            } else if (is.character(new.args$sample)) {
              if (length(new.args$sample) == 1 & 
                  new.args$sample[1] %in% getSampleIdentifiers(mutCounts)) {
                objOut <- coerceObj(x = mutCounts, to = "data.frame")[,new.args$sample]
                label <- new.args$sample
              }
            } 
            
            if(!is.null(objOut)){
              final.label <- ifelse(!is.null(new.args$main), 
                                    new.args$main, label)
              plotMutTypeProfile(mutCounts = objOut, 
                                 mutLabs = getMutationTypes(mutCounts),
                                 main = final.label,
                                 freq = new.args$freq, 
                                 ylim = new.args$ylim, 
                                 ylab = new.args$ylab, 
                                 xlab = new.args$xlab, 
                                 xaxis_cex = new.args$xaxis_cex, 
                                 cols = new.args$cols)
              
            } else {
              stop("Sample not found")
            }
          })
#---

#' @rdname msigPlot
#' @aliases msigPlot,mutSignExposures-method
setMethod("msigPlot", signature(x = "mutSignExposures"), 
          function(x, ...) {
            
            
            # Get Arguments and parse out what we want
            all.args <- list(...)
            top <- 50
            
            if (length(all.args) > 0 && 
                "top" %in% names(all.args) &&
                is.numeric(all.args$top) &&
                all.args$top > 0) {
              top <- all.args$top
            }
            
            mutCount <- coerceObj(x = x, to = "data.frame", transpose = TRUE)
            plotSignExposures(mutCount, top = top)
          })

## ~~~~~~~~~~~~~~~~~~
