
# mutSignatures
# by Damiano Fantini, Ph.D.
# http://www.mutSignatures.org

###
##### S4 classes and methods
###

## set 'getters' and 'setters' Generics -----------------------------------------------------
setGeneric("getFwkParam", function(x, label) {
  standardGeneric("getFwkParam")
})

setGeneric("getMutationTypes", function(x) {
  standardGeneric("getMutationTypes")
})

setGeneric("getSampleIdentifiers", function(x) {
  standardGeneric("getSampleIdentifiers")
})

setGeneric("getCounts", function(x) {
  standardGeneric("getCounts")
})

setGeneric("getSignatureIdentifiers", function(x) {
  standardGeneric("getSignatureIdentifiers")
})


setGeneric("setFwkParam", function(x, label, value) {
  standardGeneric("setFwkParam")
})

## set Generics for class-to-class Object coercion and other operations -----------------
setGeneric("coerceObj", function(x, to, ...) {
  standardGeneric("coerceObj")
})

setGeneric("as.mutation.counts", function(x, rownames=NULL, colnames=NULL) {
  standardGeneric("as.mutation.counts")
})

setGeneric("as.mutation.signatures", function(x) {
  standardGeneric("as.mutation.signatures")
})

setGeneric("as.mutsign.exposures", function(x, samplesAsCols = TRUE) {
  standardGeneric("as.mutsign.exposures")
})

setGeneric("setSignatureIdentifiers", function(x, names) {
  standardGeneric("setSignatureIdentifiers")
})



## mutFrameworkParams class -------------------------------------------------------------
setClass("mutFrameworkParams",
         slots = list(params = "list"))

# Define an initializer
setMethod("initialize", "mutFrameworkParams",
          function(.Object, params) {
            .Object <- callNextMethod(.Object)
            
            # Check args
            if(!is.list(params))
              stop("Bad Input")
            
            if(sum(sapply(params, length) != 1, na.rm = TRUE) > 0)
              stop("Bad input")
            
            requiredNames <- c("num_processesToExtract","num_totIterations","num_parallelCores",
                               "thresh_removeWeakMutTypes","thresh_removeLastPercent","distanceFunction",
                               "num_totReplicates","eps","stopconv",
                               "niter", "guided", "debug", 
                               "approach", "stopRule", "algorithm",
                               "logIterations")
            
            if(sum(names(params) %in% requiredNames) < 16)
              stop("Missing Params")
            
            .Object@params <- params
            .Object
          })

setMethod("show", signature(object = "mutFrameworkParams"),
          function(object) {
            out <- object@params
            cat(" mutFrameworkParams object - mutSignatures \n")
            cat(paste(" num of params included:", length(out), "\n\n"))
            linex <- lapply(1:length(out), function(i) {
              paste("   - ", names(out)[i], ": ", out[[i]], "\n", sep = "")
            })
            for (li in linex) {cat(li)}
          })

setMethod("print", signature(x = "mutFrameworkParams"),
          function(x) {
            show(object = x)
          })

setMethod("getFwkParam", signature(x = "mutFrameworkParams", label = "character"),
          function(x, label) {
            out <- x@params
            out[[label]]})

setMethod("setFwkParam", signature(x="mutFrameworkParams"),
          function(x, label, value) {
            if (length(label) != 1 | length(value) != 1)
              stop("Bad input")
            if(!is.character(label))
              stop("Bad input")
            
            x@params[[ label ]] <- value
            x
          })

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

setMethod("as.list", signature(x = "mutFrameworkParams"),
          function(x) {
            coerceObj(x, to = "list")
          })


## mutationCounts class ------------------------------------------------------
setClass("mutationCounts",
         slots = list(counts = "data.frame",
                      mutTypes = "data.frame",
                      sampleId = "data.frame"))

# Define an initializer
setMethod("initialize", "mutationCounts",
          function(.Object, x, muts, samples) {
            .Object <- callNextMethod(.Object)
            
            # Check args
            if(!is.data.frame(x) |
               !is.data.frame(muts) |
               !is.data.frame(samples))
              stop("Bad input")
            
            if(sum(duplicated(muts[,1])) > 0)
              stop("Mutation Types cannot be duplicated")
            
            if(sum(duplicated(samples[,1])) > 0)
              stop("Sample Identifiers cannot be duplicated")
            
            if (ncol(x) != nrow(samples) |
                nrow(x) != nrow(muts))
              stop("Data Dimension Issues")
            
            .Object@counts <- x
            rownames(.Object@counts) <- NULL
            colnames(.Object@counts) <- NULL
            
            .Object@mutTypes <- muts
            rownames(.Object@mutTypes) <- NULL
            colnames(.Object@mutTypes)[1] <- "mutTypes"
            
            
            .Object@sampleId <- samples
            rownames(.Object@sampleId) <- NULL
            colnames(.Object@sampleId)[1] <- "ID"
            
            .Object
          })

setMethod("show", "mutationCounts",
          function(object) {
            
            out <- object@counts
            rownames(out) <- object@mutTypes[,1]
            colnames(out) <- object@sampleId[,1]
            
            cat(" Mutation Counts object - mutSignatures")
            cat("\n\n")
            cat(paste(" Total num of MutTypes:", length(object@mutTypes[,1])))
            cat("\n")
            cat(paste(" MutTypes:", paste(head(object@mutTypes[,1], 5), collapse = ", "),
                      ifelse(length(object@mutTypes[,1]) > 5, "...", "")))
            cat("\n\n")
            cat(paste(" Total num of Samples:", length(object@sampleId[,1])))
            cat("\n")
            cat(paste(" Sample Names:", paste(head(object@sampleId[,1], 5), collapse = ", "),
                      ifelse(length(object@sampleId[,1]) > 5, "...", "")))
            
            cat("")
          })

setMethod("print", signature(x="mutationCounts"),
          function(x) {
            show(x)
          })

setMethod("getMutationTypes", "mutationCounts",
          function(x){
            base::as.vector(x@mutTypes[,1])
          })

setMethod("getSampleIdentifiers", "mutationCounts",
          function(x = "mutationCounts"){
            base::as.vector(x@sampleId[,1])
          })

setMethod("getCounts", signature(x="mutationCounts"),
          function(x) {
            
            out <- x@counts
            rownames(out) <- x@mutTypes[,1]
            colnames(out) <- x@sampleId[,1]
            out
          })

setMethod("coerceObj", signature(x = "mutationCounts", to = "character"),
          function(x, to, keepNames = TRUE) {
            if (to[1] == "data.frame") {
              getCounts(x)
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

setMethod("as.data.frame" , signature(x="mutationCounts"),
          function(x) {
            coerceObj(x, to = "data.frame")
          })

setMethod("as.matrix", signature(x = "mutationCounts"),
          function(x) {
            coerceObj(x, to = "matrix", keepNames = FALSE)
          })

setMethod("as.mutation.counts", signature(x="data.frame"),
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


setMethod("[", signature(x="mutationCounts", i="numeric"),
          function(x, i, ...) {
            
            if(nargs() > 2)
              stop("object of type 'S4' is not subsettable")

            new(Class = "mutationCounts", 
                x=base::as.data.frame(x@counts[,i]), 
                muts = x@mutTypes,
                samples=base::data.frame(x@sampleId[i,], 
                                         stringsAsFactors = FALSE,
                                         row.names = NULL))
          })

## mutationSignatures class ------------------------------------------------------
setClass("mutationSignatures",
         slots = list(mutationFreq = "data.frame",
                      mutTypes = "data.frame",
                      signatureId = "data.frame"))

# Define an initializer
setMethod("initialize", "mutationSignatures",
          function(.Object, x, muts, signNames) {
            .Object <- callNextMethod(.Object)
            
            # Check args
            if(!is.data.frame(x) |
               !is.data.frame(muts) |
               !is.data.frame(signNames))
              stop("Bad input")
            
            if(sum(duplicated(muts[,1])) > 0)
              stop("Mutation Types cannot be duplicated")
            
            if(sum(duplicated(signNames[,1])) > 0)
              stop("Sample Identifiers cannot be duplicated")
            
            if (ncol(x) != nrow(signNames) |
                nrow(x) != nrow(muts))
              stop("Data Dimension Issues")
            
            all.Nums <- base::as.numeric(format(apply(x, 2, sum), digits = 2, nsmall = 2))
            if (sum(all.Nums != 1) > 0) 
              stop ("Mutation frequency must always sum up to unity")
            
            .Object@mutationFreq <- x
            rownames(.Object@mutationFreq) <- NULL
            colnames(.Object@mutationFreq) <- NULL
            
            .Object@mutTypes <- muts
            rownames(.Object@mutTypes) <- NULL
            colnames(.Object@mutTypes)[1] <- "mutTypes"
            
            .Object@signatureId <- signNames
            rownames(.Object@signatureId) <- NULL
            colnames(.Object@signatureId)[1] <- "ID"
            
            .Object
          })

setMethod("show", "mutationSignatures",
          function(object) {
            
            out <- object@mutationFreq
            rownames(out) <- object@mutTypes[,1]
            colnames(out) <- object@signatureId[,1]
            
            cat(" Mutation Signatures object - mutSignatures")
            cat("\n\n")
            
            clMax <- ncol(out)
            cat(paste(" Total num of Signatures:", clMax))
            cat("\n")
            
            rowMax <- nrow(out)
            cat(paste(" Total num of MutTypes:", rowMax))
            cat("\n\n")
            
            nuClMax <- ifelse(clMax > 5, 5, clMax)
            colRange <- 1 : nuClMax
            #message(nuClMax)
            nuRwMax <- ifelse(rowMax > 10, 10, rowMax)
            rowRange <- 1 : nuRwMax
            
            # Headings
            cat(paste("    ", paste("Sign.", colRange, sep = "", collapse = "   "), sep = ""))
            cat("\n")
            cat(paste("    ", paste(rep("------", length(colRange)), sep = "", collapse = "   "), sep = ""))
            cat("\n")
            for(jj in rowRange){
              cat(paste("  + ", paste(format(round(out[jj,colRange], digits = 4), nsmall = 4, scientific = FALSE), sep = "", collapse = "   "), sep = ""))
              cat("  +  ")
              cat(rownames(out)[jj])
              cat("\n")
            }
            
            if(rowMax > 10){
              cat(paste("    ", paste(rep("......", length(colRange)), sep = "", collapse = "   "), sep = ""))
            }
            cat("\n")
          })

setMethod("print", signature(x="mutationSignatures"),
          function(x, ...) {
            show(x)
          })

setMethod("getMutationTypes", "mutationSignatures",
          function(x){
            base::as.vector(x@mutTypes[,1])
          })

setMethod("getSignatureIdentifiers", "mutationSignatures",
          function(x){
            base::as.vector(x@signatureId[,1])
          })

setMethod("setSignatureIdentifiers", signature(x="mutationSignatures", names = "character"),
          function(x, names){
            if (sum(duplicated(names)) == 0 &
              length(names) == nrow(x@signatureId)) {
              x@signatureId[,1] <- names
            } else {
              stop("Bad input - incompatible dimensions")
            }
            x
          })



setMethod("as.data.frame", signature(x ="mutationSignatures"),
          function(x){
            coerceObj(x, to = "data.frame")
          })



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

setMethod("[", signature(x="mutationSignatures", i="numeric"),
          function(x, i, ...) {
            
            if(nargs() > 2)
              stop("object of type 'S4' is not subsettable")
            
            new(Class = "mutationSignatures", 
                x=base::as.data.frame(x@mutationFreq[,i]), 
                muts = x@mutTypes,
                signNames=base::data.frame(x@signatureId[i,], 
                                     stringsAsFactors = FALSE,
                                     row.names = NULL))
          })

setMethod("as.list", "mutationSignatures", 
          function(x) {
            coerceObj(x, to = "list")
          })

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

setMethod("cbind2", signature(x="mutationSignatures", y="mutationSignatures"), 
          function(x, y, ...) {
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

## mutSignExposures class ------------------------------------------------------
setClass("mutSignExposures",
         slots = list(exposures = "data.frame",
                      sampleId = "data.frame",
                      signatureId = "data.frame"))

# Define an initializer
setMethod("initialize", "mutSignExposures",
          function(.Object, x, samples, signNames) {
            .Object <- callNextMethod(.Object)
            
            # Check args
            if(!is.data.frame(x) |
               !is.data.frame(samples) |
               !is.data.frame(signNames))
              stop("Bad input")
            
            if(sum(duplicated(samples[,1])) > 0)
              stop("Mutation Types cannot be duplicated")
            
            if(sum(duplicated(signNames[,1])) > 0)
              stop("Sample Identifiers cannot be duplicated")
            
            if (ncol(x) != nrow(samples) |
                nrow(x) != nrow(signNames))
              stop("Data Dimension Issues")
            
            
            .Object@exposures <- x
            rownames(.Object@exposures) <- NULL
            colnames(.Object@exposures) <- NULL
            
            .Object@sampleId <- samples
            rownames(.Object@sampleId) <- NULL
            colnames(.Object@sampleId)[1] <- "ID"
            
            .Object@signatureId <- signNames
            rownames(.Object@signatureId) <- NULL
            colnames(.Object@signatureId)[1] <- "ID"
            
            .Object
          })

setMethod("show", "mutSignExposures",
          function(object) {
            
            out <- t(object@exposures)
            rownames(out) <- object@sampleId[,1]
            colnames(out) <- object@signatureId[,1]
            
            zero.approach <- ifelse (sum(apply(out, 1, sum) > 2) == 0, TRUE, FALSE)
            
            cat(" MutSignature Exposures object - mutSignatures")
            cat("\n\n")
            
            rwMax <- nrow(out)
            cat(paste(" Total num of Samples:", rwMax))
            cat("\n")
            
            clMax <- ncol(out)
            
            nuClMax <- ifelse(clMax > 5, 5, clMax)
            colRange <- 1 : nuClMax
            #message(nuClMax)
            nuRwMax <- ifelse(rwMax > 10, 10, rwMax)
            rowRange <- 1 : nuRwMax
            
            cat(paste(" Total num of Signatures:", clMax, " { first", nuClMax,  "signatures are displayed }"))
            cat("\n")
            cat(paste(" Signature names:  ", 
                      paste(colnames(out)[1:nuClMax], collapse = ", "), 
                      ifelse(clMax > 5, " ...", ""), 
                      sep = "", collapse = " ") )
            
            cat("\n\n")
            
            
            # Headings
            cat(paste("    ", paste("Sign.", colRange, sep = "", collapse = "   "), sep = ""))
            cat("\n")
            cat(paste("    ", paste(rep("------", length(colRange)), sep = "", collapse = "   "), sep = ""))
            cat("\n")
            for(jj in rowRange){
              
              if (zero.approach) {
                prep.nums <- paste(paste(" ", format(round(out[jj,colRange],digits = 2), nsmall = 2), sep = ""), " ", sep = "")
              } else {
                prep.nums <- leadZeros(out[jj,colRange], m = 999999, char = " ", na.value = "999999")
              }
              
              prep.nums[base::as.numeric(prep.nums) > 999] <- "  >999"
              cat(paste("  + ", paste(prep.nums,
                                      sep = "", collapse = "   "), sep = ""))
              cat("  +  ")
              cat(rownames(out)[jj])
              cat("\n")
            }
            
            if(rwMax > 10){
              cat(paste("    ", paste(rep("......", length(colRange)), sep = "", collapse = "   "), sep = ""))
            }
            cat("\n")
          })

setMethod("print", signature(x="mutSignExposures"),
          function(x, ...) {
            show(x)
          })

setMethod("getSampleIdentifiers", "mutSignExposures",
          function(x){
            base::as.vector(x@sampleId[,1])
          })

setMethod("getSignatureIdentifiers", "mutSignExposures",
          function(x){
            base::as.vector(x@signatureId[,1])
          })

setMethod("setSignatureIdentifiers", signature(x="mutSignExposures", names = "character"),
          function(x, names){
            if (sum(duplicated(names)) == 0 &
                length(names) == nrow(x@signatureId)) {
              x@signatureId[,1] <- names
            } else {
              stop("Bad input - incompatible dimensions")
            }
            x
          })

setMethod("[", signature(x="mutSignExposures", i="numeric"),
          function(x, i, ...) {
            
            if(nargs() > 2)
              stop ("object of type 'S4' is not subsettable")

            new(Class = "mutSignExposures", 
                x = base::as.data.frame(x@exposures[,i]), 
                signNames = x@signatureId,
                samples=base::data.frame(x@sampleId[i,], 
                                         stringsAsFactors = FALSE, 
                                         row.names = NULL))
          })

setMethod("coerceObj", signature(x = "mutSignExposures", to = "character"),
          function(x, to, transpose = TRUE) {
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


setMethod("as.data.frame", signature(x ="mutSignExposures"),
          function(x){
            coerceObj(x, to = "data.frame")
          })

setMethod("as.data.frame", signature(x ="mutSignExposures"),
          function(x, transpose=TRUE){
            if (is.logical(transpose)) {
              coerceObj(x, to = "data.frame", transpose = transpose[1])  
            } else {
              coerceObj(x, to = "data.frame")
            }
          })

setMethod("as.mutsign.exposures", signature(x ="data.frame"),
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
                         x = base::as.data.frame(t(data)), 
                         samples = rw1, 
                         signNames = cl1)
              
            }
            out
          })

setMethod("plot", signature(x = "mutationSignatures"), 
          function(x, signature, main = NULL, ...) {
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
                                 main = final.label, ...)
              
            } else {
              stop("Signature not found")
            }
          })


setMethod("plot",signature(x = "mutationCounts"), 
          function(x, sample, main = NULL, ...) {
            mutCounts <- x
            objOut <- NULL
            if (is.numeric(sample)) {
              if(length(sample) == 1 & sample[1] <= length(getSampleIdentifiers(mutCounts))) {
                objOut <- coerceObj(x = mutCounts, to = "data.frame")[,sample]
                label <-  getSampleIdentifiers(mutCounts)[sample]
              }
            } else if (is.character(sample)) {
              if (length(sample) == 1 & sample[1] %in% getSampleIdentifiers(mutCounts)) {
                objOut <- coerceObj(x = mutCounts, to = "data.frame")[,sample]
                label <- sample
              }
            } 
            
            if(!is.null(objOut)){
              final.label <- ifelse(!is.null(main), main, label)
              plotMutTypeProfile(mutCounts = objOut, 
                                 mutLabs = getMutationTypes(mutCounts),
                                 main = final.label,
                                 ...)
              
            } else {
              stop("Sample not found")
            }
          })


setMethod("plot", signature(x = "mutSignExposures"), 
          function(x, top = 100) {
            mutCount <- coerceObj(x = x, to = "data.frame")
            plotMutCount(mutCount, top = top)
          })

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

###
##### Custom Functions - core package
###

addWeak <- function (mutationTypesToAddSet, 
                     processes_I, 
                     processesStd_I, 
                     Wall_I, 
                     genomeErrors_I, 
                     genomesReconstructed_I) 
{
  if (length(mutationTypesToAddSet) > 0 & mutationTypesToAddSet[1] > 0) {
    totalMutTypes <- nrow(Wall_I) + length(mutationTypesToAddSet)
    processes <- matrix(0, nrow = totalMutTypes, ncol = ncol(processes_I))
    processesStd <- matrix(0, nrow = totalMutTypes, ncol = ncol(processesStd_I))
    Wall <- matrix(0, nrow = totalMutTypes, ncol = ncol(Wall_I))
    genomeErrors <- lapply(1:length(genomeErrors_I), (function(i) {
      matrix(0, nrow = totalMutTypes, ncol = ncol(genomeErrors_I[[1]]))
    }))
    genomesReconstructed <- lapply(1:length(genomesReconstructed_I), 
                                   (function(i) {
                                     matrix(0, 
                                            nrow = totalMutTypes, 
                                            ncol = ncol(genomesReconstructed_I[[1]]))
                                   }))
    origArrayIndex <- 1
    for (i in 1:totalMutTypes) {
      #message(i)
      if (!(i %in% mutationTypesToAddSet)) {
        processes[i, ] <- processes_I[origArrayIndex, ]
        processesStd[i, ] <- processesStd_I[origArrayIndex,]
        Wall[i, ] <- Wall_I[origArrayIndex, ]
        for (j in 1:length(genomeErrors_I)) {
          genomeErrors[[j]][i, ] <- genomeErrors_I[[j]][origArrayIndex,  ]
        }
        for (j in 1:length(genomesReconstructed_I)) {
          genomesReconstructed[[j]][i, ] <- genomesReconstructed_I[[j]][origArrayIndex,   ]
        }
        origArrayIndex <- origArrayIndex + 1
      }
    }
  }
  else {
    processes <- processes_I
    processesStd <- processesStd_I
    Wall <- Wall_I
    genomeErrors <- genomeErrors_I
    genomesReconstructed <- genomesReconstructed_I
  }
  weakAdded.list <- list()
  weakAdded.list$processes <- processes
  weakAdded.list$processesStd <- processesStd
  weakAdded.list$Wall <- Wall
  weakAdded.list$mutCountErrors <- genomeErrors
  weakAdded.list$mutCountReconstructed <- genomesReconstructed
  return(weakAdded.list)
}


bootstrapCancerGenomes <- function (genomes) 
{
  genome.col.sums <- apply(genomes, 2, sum)
  norm.genomes <- genomes/matrix(genome.col.sums, 
                                 ncol = ncol(genomes), 
                                 nrow = nrow(genomes), 
                                 byrow = TRUE)
  bootstrapGenomes <- sapply(1:length(genome.col.sums), (function(i) {
    stats::rmultinom(1, size = genome.col.sums[i], prob = norm.genomes[,i])
  }))
  return(bootstrapGenomes)
}


evaluateStability <- function (wall, 
                               hall, 
                               params) 
{
  BIG_NUMBER <- 100
  CONVERG_ITER <- 10
  CONVERG_CUTOFF <- 0.005
  TOTAL_INIT_CONDITIONS <- 5
  num_processesToExtract <- params$num_processesToExtract
  num_totReplicates <- params$num_totReplicates
  distanceFunction <- params$distanceFunction
  minClusterDist <- BIG_NUMBER
  totalIter <- ncol(wall)/num_processesToExtract
  idx = matrix(0, nrow = nrow(hall), ncol = 1)
  clusterCompactness <- matrix(0, 
                               nrow = num_processesToExtract, 
                               ncol = totalIter)
  iStartDataSet = seq(1, ncol(wall), by = num_processesToExtract)
  iStartingDataSet = iStartDataSet[sample(1:totalIter)]
  for (iInitData in 1:min(c(TOTAL_INIT_CONDITIONS, totalIter))) {
    #message(iInitData)
    iStartingData <- iStartingDataSet[iInitData]
    iEnd <- iStartingData + num_processesToExtract - 1
    centroids <- cbind(wall[, iStartingData:iEnd])
    centroidsTest <- sapply(1:ncol(centroids), (function(kk) {
      runif(nrow(centroids))
    }))
    countIRep <- 0
    for (iRep in 1:num_totReplicates) {
      #message(iRep)
      tmp.tab <- t(cbind(centroids, wall))
      tmp.pdist <- base::as.vector(proxy::dist(tmp.tab, distanceFunction))
      if (num_processesToExtract > 1)
        tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
      allDist <- pracma::squareform(tmp.pdist)
      
      cd.colRange <- (ncol(centroids) + 1):ncol(allDist)
      centroidDist = t(allDist[1:ncol(centroids), cd.colRange])
      
      jRange <- sort(1:num_processesToExtract)
      for (jIndex in 1:num_processesToExtract) {
        j <- jRange[jIndex]
        for (i in seq(1, ncol(wall), by = num_processesToExtract)) {
          #message(i)
          iRange = i:(i + num_processesToExtract - 1)
          tmp.min <- min(centroidDist[iRange, j], na.rm = TRUE)
          Ind <- which(centroidDist[iRange, j] == tmp.min)[1]
          centroidDist[iRange[Ind], ] <- BIG_NUMBER
          idx[iRange[Ind], 1] <- j
        }
      }
      maxDistToNewCentroids <- 0
      for (i in 1:ncol(centroids)) {
        tmp.dset <- wall[, base::as.vector(idx == i)]
        centroids[, i] <- apply(tmp.dset, 1, mean)
        tmp.dset <- t(cbind(centroids[, i], 
                            centroidsTest[,i]))
        tmp.pdist <- base::as.vector(proxy::dist(tmp.dset, distanceFunction))
        tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
        maxDistToNewCentroids <- max(maxDistToNewCentroids, 
                                     tmp.pdist, 
                                     na.rm = TRUE)
      }
      if (maxDistToNewCentroids < CONVERG_CUTOFF) {
        countIRep <- countIRep + 1
      }
      else {
        countIRep <- 0
        centroidsTest <- centroids
      }
      if (countIRep == CONVERG_ITER) {
        break
      }
    }
    for (i in 1:ncol(centroids)) {
      tmp.tab <- t(cbind(centroids[, i], 
                         wall[, base::as.vector(idx == i)]))
      tmp.pdist <- base::as.vector(proxy::dist(tmp.tab, distanceFunction))
      tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
      clusterDist <- pracma::squareform(tmp.pdist)
      clusterCompactness[i, ] = clusterDist[1, 2:ncol(clusterDist)]
    }
    dist.test <- apply(clusterCompactness, 2, (function(clm) {
      mean(clm, na.rm = TRUE)
    }))
    if (sum(minClusterDist > dist.test) == length(dist.test)) {
      centroidsFinal <- centroids
      idxFinal <- idx
      clusterCompactnessFinal <- clusterCompactness
    }
  }
  centroids <- t(centroidsFinal)
  idx <- idxFinal
  clusterCompactness <- clusterCompactnessFinal
  centDist <- apply(clusterCompactness, 1, (function(tmprw) {
    mean(tmprw, na.rm = TRUE)
  }))
  centDistInd <- order(centDist, decreasing = FALSE)
  clusterCompactness <- clusterCompactness[centDistInd, ]
  centroids <- centroids[centDistInd, ]
  idxNew <- idx
  for (i in 1:num_processesToExtract) {
    idxNew[base::as.vector(idx == centDistInd[i]), 1] <- i
  }
  idx <- idxNew
  if (num_processesToExtract > 1) {
    processStab <- silhouetteMLB(data = t(wall), fac = idx, 
                                 distanceFunction)
    processStabAvg <- matrix(0, nrow = 1, ncol = num_processesToExtract)
    for (i in 1:num_processesToExtract) {
      processStabAvg[1, i] = mean(processStab[idx == i])
    }
  } else {
    # Adjusted params for a 1-class silhouette!
    tmp.tab <- t(cbind(centroids, wall))
    tmp.pdist <- base::as.vector(proxy::dist(tmp.tab, distanceFunction))
    tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
    allDist <- pracma::squareform(tmp.pdist)
    processStab <- 1 - t(allDist[1, 2:ncol(allDist)])
    
    # Silhouette plot
    xrange <- c(min(processStab), max(processStab))
    xrange[1] <- ifelse(xrange[1] > 0, 0, (-1.15) * abs(xrange[1]))
    xrange[2] <- 1.15
    graphics::barplot(sort(base::as.numeric(processStab), decreasing = TRUE), 
                      col = "gray20",
                      xlim = xrange,
                      horiz = TRUE, xlab = "Silhouette Value", 
                      ylab = "", 
                      main = "Silhouette Plot", 
                      border = "gray20")
    graphics::abline(v=0)
    graphics::title(ylab="Iter. Results (by Group)", line=1, cex.lab=1, font = 2)
    
    processStabAvg <- apply(processStab, 2, (function(clmn) {
      base::mean(clmn, na.rm = TRUE)
    }))
  }
  
  if (num_processesToExtract > 1) {
    centroidStd <- matrix(0, nrow = nrow(centroids), ncol = ncol(centroids))
  } else {
    centroidStd <- matrix(0, ncol = length(centroids), nrow = 1)
  }
  
  for (i in 1:num_processesToExtract) {
    centroidStd[i, ] <- apply(wall[, idx == i], 1, (function(rw) {
      stats::sd(rw, na.rm = TRUE)
    }))
  }
  centroids <- t(cbind(centroids))
  centroidStd <- t(centroidStd)
  idxS <- matrix(0, nrow = length(idx), ncol = 1)
  for (i in seq(1, ncol(wall), by = num_processesToExtract)) {
    iEnd <- i + num_processesToExtract - 1
    idxG <- idx[i:iEnd]
    for (j in 1:num_processesToExtract) {
      idxS[(i + j - 1), ] = which(idxG == j)
    }
  }
  exposure <- matrix(0, nrow = max(idxS), ncol(hall))
  exposureStd <- matrix(0, nrow = max(idxS), ncol(hall))
  for (i in 1:max(idxS)) {
    exposure[i, ] <- apply(hall[idx == i, ], 2, (function(cl) {
      mean(cl, na.rm = TRUE)
    }))
    exposureStd[i, ] <- apply(hall[idx == i, ], 2, (function(cl) {
      sd(cl, na.rm = TRUE)
    }))
  }
  
  # Fix to sign.to.extract.num = 1
  if (num_processesToExtract < 2){
    centroids <- t(centroids)
  }
  
  
  result.list <- list()
  result.list$centroids <- centroids
  result.list$centroidStd <- centroidStd
  result.list$exposure <- exposure
  result.list$exposureStd <- exposureStd
  result.list$idx <- idx
  result.list$idxS <- idxS
  result.list$processStab <- processStab
  result.list$processStabAvg <- processStabAvg
  result.list$clusterCompactness <- clusterCompactness
  return(result.list)
}


filterOutIterations <- function (wall, 
                                 hall, 
                                 cnt_errors, 
                                 cnt_reconstructed, 
                                 params) 
{
  num_processesToExtract <- params$num_processesToExtract
  thresh_removeLastPercent <- params$thresh_removeLastPercent
  num_totIterations <- ncol(wall)/num_processesToExtract
  tot.rm.iterations <- round(thresh_removeLastPercent * num_totIterations)
  if (tot.rm.iterations > 0) {
    closeness.mutCounts <- matrix(0, 
                                  nrow = num_totIterations, 
                                  ncol = 1)
    for (i in 1:num_totIterations) {
      closeness.mutCounts[i, ] <- base::norm(cnt_errors[[i]], "F")
    }
    indexClosenessGenomes <- order(closeness.mutCounts, 
                                   decreasing = TRUE)
    removeIterations <- indexClosenessGenomes[1:tot.rm.iterations]
    removeIterationSets <- matrix(0, 
                                  nrow = (num_processesToExtract * tot.rm.iterations), 
                                  ncol = 1)
    for (i in 1:tot.rm.iterations) {
      iStart <- num_processesToExtract * (removeIterations[i] - 1) + 1
      iEnd <- num_processesToExtract * removeIterations[i]
      tmpRowRange <- (num_processesToExtract * (i - 1) + 1):(num_processesToExtract * i)
      removeIterationSets[tmpRowRange, ] <- iStart:iEnd
    }
    wall <- wall[, -removeIterationSets]
    hall <- hall[-removeIterationSets, ]
    cnt_errors <- cnt_errors[-removeIterations]
    cnt_reconstructed <- cnt_reconstructed[-removeIterations]
  }
  res.list <- list()
  res.list$Wall <- wall
  res.list$Hall <- hall
  res.list$mutCounts.errors <- cnt_errors
  res.list$mutCounts.reconstructed <- cnt_reconstructed
  return(res.list)
}


getTestRunArgs <- function (testN = 1) 
{
  out <- list()
  if (testN == 3) {
    out$v <- cbind(c(142, 133, 1, 24, 53, 55, 4, 4, 100), 
                   c(132, 113, 0, 34, 50, 52, 3, 3, 17), 
                   c(155, 139,10, 14, 53, 45, 2, 5, 13), 
                   c(124, 156, 22, 21,52, 45, 2, 7, 100))
    out$r <- 2
    out$params <- setMutClusterParams(num_processesToExtract = 2)
    out$params$num_parallelCores <- 1
    out$params$stopconv <- 800
    out$params$niter <- 8000
  }
  else if (testN == 4) {
    out$data <- cbind(c(runif(4, 10, 25), 
                        runif(6, 20, 50), 
                        runif(3, 0, 5)), 
                      c(runif(4, 50, 60),
                        runif(6, 45, 55), 
                        runif(3, 30, 40)), 
                      c(runif(4, 12, 15), 
                        runif(6,10, 15), 
                        runif(3, 10, 12)), 
                      c(runif(4, 5, 20), 
                        runif(6, 16, 26), 
                        runif(3, 24, 29)))
    out$fac <- c(rep(1, 4), rep(2, 6), rep(3, 3))
  }
  else if (testN == 5) {
    out$W <- do.call(cbind, lapply(1:20, (function(i) {
      cbind(c(runif(4, 0.05, 0.15), 
              c(1e-15 * runif(1,1, 9)), 
              runif(3, 0.003, 0.007), 
              runif(9, 0.04, 0.09)), 
            c(runif(3, 0.08, 0.18), 
              c(1e-15 * runif(1, 1, 9)), 
              runif(1, 0.03, 0.07), 
              runif(3, 0.004, 0.009), 
              runif(9, 0.04, 0.09)))
    })))
    out$H <- do.call(rbind, lapply(1:40, (function(i) {
      c(runif(1, 0.005, 0.099), 
        runif(2, 50, 800), 
        runif(2, 0.005, 0.099), 
        runif(1, 50, 750), 
        runif(2, 0.005, 0.099), 
        runif(1, 20, 800))
    })))
    out$params$analyticApproach <- "denovo"
  }
  else if (testN == 6) {
    tmut <- runif(10, 150, 1350)
    eff1 <- runif(10, 0.45, 0.89)
    eff2 <- 1 - eff1
    out$exposures <- sapply(1:10, (function(i) {
      c(eff1[i], eff2[i]) * tmut[i]
    }))
  }
  else {
    set.seed(999)
    my.mat <- sapply(1:10, (function(i) {
      c(base::as.integer(runif(3, 80, 150)), 
        base::as.integer(runif(7, 0, 10)), 
        base::as.integer(runif(9, 40, 80)), 
        base::as.integer(runif(1, 0, 3)), 
        base::as.integer(runif(10, 60, 120)))
    }))
    rownames(my.mat) <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", 
                          "A[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", 
                          "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "A[T>A]A", 
                          "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C", 
                          "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", 
                          "A[T>G]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", 
                          "C[C>G]A", "C[C>G]C")
    #out$mutCount.obj <- setMutCountObject(mutCountMatrix = my.mat)
    #out$params <- setMutClusterParams(num_processesToExtract = 2)
    out$params$num_totIterations <- 1
    out$params$num_parallelCores <- 1
    out$params$stopconv <- 800
    out$params$niter <- 8000
    if (testN == 2) {
      out$params$stopconv <- 2000
      out$params$niter <- 15000
      out$params$num_totIterations <- 3
    }
  }
  return(out)
}


leadZeros <- function (n, m, char = "0", na.value = NA) 
{
  max.zeros <- nchar(base::as.character(round(m)))
  tmp.digits <- nchar(base::as.character(round(n)))
  zeros.toAdd <- max.zeros - tmp.digits
  
  returnVect <- sapply(1:length(n), function(i){
    if (zeros.toAdd[i] >= 0) {
      paste(c(rep(char, zeros.toAdd[i]), base::as.character(round(n[i]))), sep = "", collapse = "")    
    } else {
      na.value
    }
  })
  return(returnVect) 
}


removeWeak <- function (input_mutCounts, 
                        params) 
{
  thresh_removeWeakMutTypes <- params$thresh_removeWeakMutTypes
  sum.counts <- apply(input_mutCounts, 1, sum)
  sum.counts.idx <- order(sum.counts, decreasing = FALSE)
  sorted.sum.counts <- sum.counts[sum.counts.idx]
  tot.mut.counts <- sum(input_mutCounts)
  tot.muttypes.toremove <- sum((sapply(1:length(sorted.sum.counts), 
                                       (function(i) {
                                         sum(sorted.sum.counts[1:i])
                                       }))/tot.mut.counts) < thresh_removeWeakMutTypes)
  return.list <- list()
  if (tot.muttypes.toremove > 0) {
    removed.mutset <- sum.counts.idx[c(1:tot.muttypes.toremove)]
    input_mutCounts <- input_mutCounts[-removed.mutset, ]
    return.list$removed.mutset <- removed.mutset
  }
  else {
    return.list$removed.mutset <- (-1)
  }
  return.list$output.mutCounts <- input_mutCounts
  return(return.list)
}


silhouetteMLB <- function (data, 
                           fac, 
                           method = "cosine", 
                           plot = TRUE) 
{
  if (nrow(data) != length(fac)) 
    stop("Bad input!")
  dist.matrix <- base::as.matrix(proxy::dist(x = data, method = method))
  sil.check <- cluster::silhouette(x = base::as.numeric(base::as.factor(fac)), 
                                   dist = dist.matrix)
  if (plot == TRUE) {
    tmp <- lapply(unique(sil.check[, 1]), (function(clid) {
      part.out <- sil.check[sil.check[, 1] == clid, ]
      part.out[order(part.out[, 3], decreasing = TRUE), 
               ]
    }))
    tmp <- do.call(rbind, tmp)
    xrange <- c(min(tmp[,3]), max(tmp[,3]))
    xrange[1] <- ifelse(xrange[1] > 0, 0, (-1.15) * abs(xrange[1]))
    xrange[2] <- 1.15
    graphics::barplot(tmp[nrow(tmp):1, 3], 
                      col = base::as.factor(tmp[nrow(tmp):1,1]),
                      xlim = xrange,
                      horiz = TRUE, xlab = "Silhouette Value", 
                      ylab = "", 
                      main = "Silhouette Plot", 
                      border = base::as.factor(tmp[nrow(tmp):1,1]))
    graphics::abline(v=0)
    graphics::title(ylab="Iter. Results (by Group)", line=1, cex.lab=1, font = 2)
  }
  return(base::as.vector(sil.check[, 3]))
}


alexaNMF <- function(v, r, params)
{
  # Brunet algorithm, from Alexandrov NMF WTSI code
  # define params (hard-set)
  debug <- params$debug
  chk.step <- 50
  dot.eachSteps <- 2000
  
  # retireve user-defined params
  eps <- params$eps
  num.processes <- r
  stopconv <- params$stopconv
  niter <- params$niter
  err.threshold <- 1e-10
  stopRule <- ifelse(params$stopRule == "LA", "LA", "DF")
  
  # set.seed(231082)
  v <- base::as.matrix(v)
  rownames(v) <- NULL
  colnames(v) <- NULL
  
  # Double check input matrix
  if (min(v) < 0)
    stop("Matrix entries cannot be negative")
  if (min(apply(v, 1, sum)) == 0)
    stop("Entries cannot all be equal to 0")
  
  # Initialize W0 and H0
  W.k <- do.call(cbind, lapply(1:num.processes, (function(i){
    out.o <- runif(n = nrow(v), min = eps, max = 100)
    out.o/sum(out.o)
  })))
  H.k <- matrix((1/num.processes), nrow = num.processes, ncol = ncol(v))
  
  # Initialize looping vars
  itr <- 1
  chk.j <- 1
  stationary.chk <- 0
  force.out <- 1
  
  # Debugging plot
  if (debug)
    graphics::plot(-10, xlim = c(1000, niter),
                   ylim = c( ifelse(stopRule == "DF", (0.1*err.threshold), 0.001) ,
                             ifelse(stopRule == "DF", 10, ncol(H.k))), log = "xy",
                   xlab = "Iteration", ylab = "Variation", main = "Convergence")
  
  # Initialize the objects for comparing dissimilarity
  # DF approach
  cons.old <- base::as.vector(W.k)
  # Alexandrov approach
  consold <- matrix(0, nrow = ncol(H.k), ncol = ncol(H.k))
  
  while (itr < niter){
    
    if (itr %% dot.eachSteps == 0) {
      if (stationary.chk > chk.step) {
        message(":", appendLF = FALSE)
      } else {
        message(".", appendLF = FALSE)
      }
    }
    
    delta.01 <- apply(W.k, 2, sum)
    H.k <- H.k * (t(W.k) %*% (v/(W.k %*% H.k))) / delta.01
    H.k[H.k < eps] <- eps
    
    W.tmp <- W.k * ((v/(W.k %*% H.k)) %*% t(H.k))
    W.k <- do.call(cbind, lapply(1:ncol(W.tmp), (function(ci){
      W.tmp[,ci] / sum(H.k[ci,])
    })))
    W.k[W.k<eps] <- eps
    
    # check convergence every 'chk.step' iterations
    if (itr > stopconv & itr %% chk.step == 0 & stopRule == "DF") {
      chk.j <- chk.j + 1
      H.k[H.k < eps] <- eps
      W.k[W.k < eps] <- eps
      cons <- base::as.vector(W.k)
      
      # compare to consold and reorder
      dist.measure <- proxy::dist(rbind(cons, cons.old), method = "cosine") [1]
      cons.old <- cons
      
      if (debug)
        points(itr, (dist.measure + (err.threshold*0.1)), pch = 19, col = "red2")
      
      # evaluate distance
      if (dist.measure < err.threshold) {
        stationary.chk <- stationary.chk + 1
      } else {
        stationary.chk <- 0
      }
      if (stationary.chk > (stopconv / chk.step)) {
        force.out <- 0
        message("$", appendLF = FALSE)
        break()
      }
    } else if (itr > stopconv & itr %% chk.step == 0 & stopRule == "LA") {
      chk.j <- chk.j + 1
      H.k[H.k < eps] <- eps
      W.k[W.k < eps] <- eps
      y <- apply(H.k, 2, max)
      index <- apply(H.k, 2, (function(dt) {
        which.max(dt)[1]
      }))
      mat1 = t(sapply(1:ncol(H.k), (function(ii) {
        index
      })))
      mat2 = sapply(1:ncol(H.k), (function(ii) {
        index
      }))
      cons <- mat1 == mat2
      if (sum(cons != consold) == 0) {
        stationary.chk <- stationary.chk + 1
      } else {
        stationary.chk <- 0
      }
      #
      consold <- cons
      
      if (debug)
        points(itr, (sum(cons != consold) / 100), pch = 19, col = "red2")
      
      if (stationary.chk > (stopconv/chk.step)) {
        force.out <- 0
        message("$", appendLF = FALSE)
        break()
      }
    }
    itr <- itr + 1
  }
  if (force.out == 1) {
    message("!", appendLF = FALSE)
  }
  output <- list()
  output$w <- W.k
  output$h <- H.k
  return(output)
}


attachMutType <- function(mutData,
                          ref_colName = "reference_allele",
                          var_colName = "variant_allele",
                          var2_colName = NULL,
                          context_colName = "context",
                          format = 1,
                          mutType_dict = "alexa",
                          mutType_colName = "mutType")
{
  # Validate input data and fields in mutData
  if(!((is.data.frame(mutData) | is.matrix(mutData) ) &
       sum(c(ref_colName, var_colName, context_colName) %in% colnames(mutData)) == 3 ))
    stop ("Issue with the input dataset. Make sure to feed in a data.frame or
          a matrix and double check the name of the fields pointing to chromosome
          name, start and end positions")
  if (!(format %in% c(1,2)))
    stop ("Please, specify a valid format number (example: 1)")
  if (! (is.null(var2_colName))) {
    if (!(var2_colName %in% colnames(mutData)))
      stop ("Invalid var2 column")
  }
  
  if (!is.character(mutType_colName) |
      length(mutType_colName) > 1)
    stop("Bad mutType_colName")
  if(mutType_colName %in% colnames(mutData))
    stop ("mutType_colName already exists as column name in the current dataset")
  
  # convert factors to chars
  mutData <- data.frame(mutData, stringsAsFactors = FALSE, row.names = NULL)
  my.key.cols <- c(ref_colName, var_colName, var2_colName, context_colName)
  my.key.cols <- my.key.cols[!is.na(my.key.cols)]
  for (clmn in my.key.cols) {
    mutData[,clmn] <- base::as.character(base::as.vector(mutData[,clmn])) 
  }
  
  message("Assigning mutation types ", appendLF = FALSE)
  
  mutData[,mutType_colName] <- sapply(1:nrow(mutData), (function(i){
    
    if (nrow(mutData) > 1000 & i %in% base::as.integer(seq(1, nrow(mutData),length.out = 20)))
      message(".", appendLF = FALSE)
    
    # first, extract elems and check for middle base to match the reference
    ctx.len <- nchar(mutData[i,context_colName])
    half.ln <- (ctx.len - 1) / 2
    mid.seq <- substr(mutData[i,context_colName], (half.ln + 1), (half.ln + 1))
    pre.seq <- substr(mutData[i,context_colName], 1, half.ln)
    post.seq <- substr(mutData[i,context_colName], (half.ln + 2), ctx.len)
    
    if(mid.seq !=  mutData[i,ref_colName] | 
       (is.null(var2_colName) & mid.seq == mutData[i,var_colName]) |
       (tryCatch({(mid.seq == mutData[i,var_colName] & mid.seq == mutData[i,var2_colName])},
                 error = function(e) {FALSE}))){
      # no match means --> NA
      mut.base <- NA
    } else {
      #mut.variant to use in case var2 is specified
      if (mutData[i,ref_colName] != mutData[i,var_colName]) {
        mut.base <- mutData[i,var_colName]
      } else if (!is.null(var2_colName) ){
        if ( mutData[i,ref_colName] != mutData[i,var2_colName]) {
          mut.base <- mutData[i,var2_colName]
        } else {
          mut.base <- NA  
        }
      } else {
        mut.base <- NA  
      }
      
      # match, format and return according to a standard format (for now)
      if (is.na(mut.base)) {
        NA
      } else {
        paste(mid.seq, ".", mut.base, "[", pre.seq, mid.seq, post.seq, "][", pre.seq, mut.base, post.seq, "]", sep = "", collapse = "")
      }
    }
  }))
  message(". Done!", appendLF = TRUE)
  
  if (sum(is.na(mutData[,mutType_colName])) > 0) {
    message(paste("Removing",sum(is.na(mutData[,mutType_colName])), "positions."))
    mutData <- mutData[!is.na(mutData[,mutType_colName]),]
  }
  
  # Now, apply revCompl transformation
  message("Now applying RevCompl transformation", appendLF = FALSE)
  if (mutType_dict == "alexa") {
    idx <- grep("^((G|A)\\.)", mutData[,mutType_colName])
    mutData[idx,mutType_colName] <- sapply(mutData[idx,mutType_colName], (function(seq){
      
      base.wt  <- revCompl(gsub("\\..+$", "", seq))
      base.mut <- revCompl(gsub("^.+\\.", "", gsub("\\[.+$", "", seq)))
      seq.wt   <- revCompl(gsub("^.+\\[", "", gsub("\\]\\[.+$", "", seq)))
      seq.mut  <- revCompl(gsub("^.+\\]\\[", "", gsub("\\]$", "", seq)))
      paste(base.wt,".",base.mut, "[", seq.wt, "][", seq.mut, "]", sep = "", collapse = "")
    }))
  } else if (mutType_dict == "custom") {
    idx <- grep("^((G|T)\\.)", mutData[,mutType_colName])
    mutData[idx,mutType_colName] <- sapply(mutData[idx,mutType_colName], (function(seq){
      base.wt  <- revCompl(gsub("\\..+$", "", seq))
      base.mut <- revCompl(gsub("^.+\\.", "", gsub("\\[.+$", "", seq)))
      seq.wt   <- revCompl(gsub("^.+\\[", "", gsub("\\]\\[.+$", "", seq)))
      seq.mut  <- revCompl(gsub("^.+\\]\\[", "", gsub("\\]$", "", seq)))
      paste(base.wt,".",base.mut, "[", seq.wt, "][", seq.mut, "]", sep = "", collapse = "")
    }))
  }
  
  message(". Done!", appendLF = TRUE)
  message("Final formatting", appendLF = FALSE)
  #Attach the format of interest
  mutData[,mutType_colName] <- sapply(mutData[,mutType_colName], (function(seq){
    base.wt  <- gsub("\\..+$", "", seq)
    base.mut <- gsub("^.+\\.", "", gsub("\\[.+$", "", seq))
    seq.wt   <- gsub("^.+\\[", "", gsub("\\]\\[.+$", "", seq))
    seq.mut  <- gsub("^.+\\]\\[", "", gsub("\\]$", "", seq))
    half.len <- (nchar(seq.wt) - 1 ) / 2
    pre.seq <- substr(seq.wt, 1, half.len)
    post.seq <-substr(seq.wt, half.len + 2, nchar(seq))
    
    if (format == 1) {
      # --> N[N>M]N
      paste(pre.seq, "[", base.wt, ">", base.mut, "]", post.seq, sep = "", collapse = "")
    } else if (format == 2) {
      # --> NN.N>M
      paste(pre.seq, post.seq, ".", base.wt, ">", base.mut, sep = "", collapse = "")
    }  else {
      # --> N.M[NNN][NMN]
      paste(base.wt,".",base.mut, "[", seq.wt, "][", seq.mut, "]", sep = "", collapse = "")
    }
  }))
  
  message(". Done!", appendLF = TRUE)
  
  return(mutData)
}


chihJenNMF <- function(v, r, params) {
  
  # http://ieeexplore.ieee.org/document/4359171/
  # alternative approach for NMF
  # define params (hard-set)
  debug <- params$debug
  chk.step <- 50
  dot.eachSteps <- 2000
  
  # retireve user-defined params
  eps <- params$eps
  num.processes <- r
  stopconv <- params$stopconv
  niter <- params$niter
  err.threshold <- 1e-10
  
  # set.seed(231082)
  v <- base::as.matrix(v)
  rownames(v) <- NULL
  colnames(v) <- NULL
  
  # Double check input matrix
  if (min(v) < 0)
    stop("Matrix entries cannot be negative")
  if (min(apply(v, 1, sum)) == 0)
    stop("Entries cannot all be equal to 0")
  
  # Debugging plot
  if (debug)
    graphics::plot(-10, xlim = c(1000, niter), ylim = c( (0.1*err.threshold), 10), log = "xy", 
                   xlab = "Iteration", ylab = "Variation", main = "Convergence")
  
  # Initialize W0 and H0
  W.k <- do.call(cbind, lapply(1:num.processes, (function(i){
    out.o <- runif(n = nrow(v), min = eps, max = 100)
    out.o/sum(out.o)
  })))
  H.k <- matrix((1/num.processes), nrow = num.processes, ncol = ncol(v))
  
  # Initialize looping vars
  itr <- 1
  chk.j <- 1
  cons.old <- base::as.vector(W.k)
  stationary.chk <- 0
  force.out <- 1
  while (itr < niter){
    
    if (itr %% dot.eachSteps == 0) {
      if (stationary.chk > chk.step) {
        message(":", appendLF = FALSE)
      } else {
        message(".", appendLF = FALSE)
      }
    }
    
    WtW <- t(W.k) %*% W.k
    gradH <- ((WtW %*% H.k) - (t(W.k) %*% v))
    H.b <- H.k; H.b[H.b<eps] <- eps;
    H.k <- H.k - (H.b / ((WtW %*% H.b) + eps)) * gradH
    HHt <- H.k %*% t(H.k)
    gradW <- (W.k %*% HHt) - (v %*% t(H.k))
    
    W.b <- W.k; W.b[W.b < eps] <- eps
    W.k <- W.k - (W.b / ((W.b %*% HHt) + eps)) * gradW
    S <- apply(W.k, 2, sum)
    
    W.k <- do.call(cbind, lapply(1:ncol(W.k), (function(ci){
      W.k[,ci] / S[ci]
    })))
    H.k <- do.call(rbind, lapply(1:nrow(H.k), (function(ri){
      H.k[ri,] * S[ri]
    })))
    
    # optional ? Keep as is for now
    H.k <- do.call(cbind, lapply(1:ncol(H.k), (function(ci){
      H.k[,ci] / sum(H.k[,ci])
    })))
    
    # Final non-negative check
    H.k[H.k < eps] <- eps
    W.k[W.k < eps] <- eps
    
    # check convergence every 'chk.step' iterations
    if (itr > stopconv & itr %% chk.step == 0) {
      chk.j <- chk.j + 1
      W.k[W.k < eps] <- eps
      cons <- base::as.vector(W.k)
      # compare to consold and reorder
      dist.measure <- proxy::dist(rbind(cons, cons.old), method = "cosine") [1]
      cons.old <- cons
      
      if (debug)
        points(itr, (dist.measure + (err.threshold*0.1)), pch = 19, col = "red2")
      
      # evaluate distance
      if (dist.measure < err.threshold) {
        stationary.chk <- stationary.chk + 1
      } else {
        stationary.chk <- 0
      }
      if (stationary.chk > (stopconv / chk.step)) {
        force.out <- 0
        message("$", appendLF = FALSE)
        break()
      }
    }
    itr <- itr + 1
  }
  if (force.out == 1) {
    message("!", appendLF = FALSE)
  }
  output <- list()
  output$w <- W.k
  output$h <- H.k
  return(output)
}


countMutTypes <- function(mutTable, 
                          mutType_colName = "mutType", #mutTypeLab 
                          sample_colName = NULL)  #sampleLab
{
  
  if (!(is.data.frame(mutTable) | is.matrix(mutTable)))
    stop ("Bad input")
  if (! mutType_colName %in% colnames(mutTable))
    stop ("mutType field not found")
  if (!is.null(sample_colName)){
    if (! sample_colName %in% colnames(mutTable))
      stop ("sample_colName field not found")
  } 
  
  mutType.labels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
                      "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
                      "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
                      "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                      "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
                      "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
                      "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
                      "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
                      "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  custPatt01 <- "^(A|C|G|T)\\[(A|C|G|T)>(A|C|G|T)\\](A|C|G|T)$"
  if (sum(regexpr(custPatt01, mutTable[,mutType_colName]) > 0) != length(mutTable[,mutType_colName]))
    stop ("Problem with the mutation type format... Please use the following format: A[C>G]T")
  
  # Force Alexandrov-style mut types
  idx.to.fix <- which(!mutTable[,mutType_colName] %in% mutType.labels)
  if (length(idx.to.fix) > 0) {
    
    tmp <- mutTable[,mutType_colName][idx.to.fix]
    corrected.mutTypes <- sapply(1:length(tmp), (function(i){
      out <- c(revCompl(substr(tmp[i], 7,7)), "[",
               revCompl(substr(tmp[i], 3,3)), ">",
               revCompl(substr(tmp[i], 5,5)), "]",
               revCompl(substr(tmp[i], 1,1)))
      paste(out, collapse = "", sep = "")
    }))
    mutTable[,mutType_colName][idx.to.fix] <- corrected.mutTypes
  }
  
  if (sum(!mutTable[,mutType_colName] %in% mutType.labels) > 0)
    stop("Problem with the mutType... Please check the input")
  
  if(is.null(sample_colName)) {
    my.mutTypes <- mutTable[,mutType_colName]
    out.1 <- sapply(mutType.labels, (function(mtt){
      sum(my.mutTypes == mtt)
    }))
    out.1 <- data.frame(cbind(sample=out.1))
  } else {
    unique.cases <- unique(mutTable[,sample_colName])
    out.1 <- lapply(unique.cases, (function(csid){
      tmp.df <- mutTable[mutTable[,sample_colName] == csid,]
      out.3 <- sapply(mutType.labels, (function(mtt){
        sum(tmp.df[,mutType_colName] == mtt)
      }))
      out.3 <- data.frame(cbind(out.3))
      colnames(out.3) <- csid
      out.3
    }))
    out.1 <- data.frame(do.call(cbind, out.1), stringsAsFactors = FALSE)
    tryCatch({colnames(out.1) <- unique.cases}, error = function(e) NULL)
  }
  
  # Prepare mutationCounts object
  fOUT <- new(Class = "mutationCounts", 
              x = out.1, 
              muts = data.frame(mutTypes = rownames(out.1), stringsAsFactors = FALSE), 
              samples = data.frame(ID = colnames(out.1), stringsAsFactors = FALSE))
  return(fOUT)
}


decipherMutationalProcesses <- function (input, 
                                         params)
{
  if (class(params) == "mutFrameworkParams" & class(input) == "mutationCounts") {
    paramsList <- coerceObj(params, to = "list")
    inputMAT <- coerceObj(input, to = "matrix", keepNames = FALSE) 
  } else {
    stop("Malformed Input")
  }
  
  if (paramsList$approach != "counts") {
    freq.input <- frequencize(inputMAT)
    inputMAT <- freq.input$freqs
    inputColsums <- freq.input$colSums
  }
  currentWarnings <- options()$warn
  options(warn = -1)
  #requiredElements <- c("cancerType", "mutCounts", "sampleNames", "mutTypes")
  #if (!("input" %in% ls()) | !(is.list(input)) | sum(requiredElements %in%
  #                                                   names(input)) != length(requiredElements))
  #  stop("Malformed Input / Input does not include all the required fields!")
  
  if (is.numeric(paramsList$num_processesToExtract)) {
    paramsList$analyticApproach <- "denovo"
  } else {
    stop("An error occurred!")
  }
  if (paramsList$analyticApproach == "denovo") {
    deconvData <- deconvoluteMutCounts(input_mutCounts = inputMAT,
                                       params = paramsList)
  } else {
    stop("An error occurred!")
  }
  
  # Package results and send out
  mutProcesses <- list()
  
  # Results first
  final.proc <- data.frame(deconvData$processes, stringsAsFactors = FALSE)
  colnames(final.proc) <- paste("Sign.",
                                sapply(1:ncol(final.proc), (function(n){
                                  leadZeros(n, (10*ncol(final.proc)))
                                })),  sep = "")
  rownames(final.proc) <- getMutationTypes(input)
  
  # New signatures
  signResult <- mutSignatures::as.mutation.signatures(final.proc)
  
  
  final.expo <- data.frame(deconvData$exposure, stringsAsFactors = FALSE)
  if (getFwkParam(params, "approach") != "counts"){
    final.expo <- data.frame(sapply(1:ncol(final.expo), function(cjj){
      inputColsums[cjj] * final.expo[,cjj] / sum(final.expo[,cjj])
    }), row.names = NULL, stringsAsFactors = FALSE)  
  }
  tryCatch({colnames(final.expo) <- getSampleIdentifiers(input)}, error = function(e) {NULL})
  rownames(final.expo) <- colnames(final.proc)
  
  # New exposures
  expoResult <- mutSignatures::as.mutsign.exposures(x = final.expo) 
  
  mutProcesses$Results <- list()
  mutProcesses$Results$signatures <- signResult
  mutProcesses$Results$exposures <- expoResult
  
  # Attach input and analysis data
  mutProcesses$RunSpecs <- list()
  mutProcesses$RunSpecs$input <- input
  mutProcesses$RunSpecs$params <- params
  
  # Attach Supplementary and Extra Controls stuff
  if (getFwkParam(params, "logIterations") != "lite") {
    mutProcesses$Supplementary <- list()
    mutProcesses$Supplementary$allProcesses <- deconvData$Wall
    mutProcesses$Supplementary$allExposures <- deconvData$Hall
    mutProcesses$Supplementary$idx <- deconvData$idx
    mutProcesses$Supplementary$mutCountErrors <- deconvData$mutCountErrors
    mutProcesses$Supplementary$mutCountReconstructed <- deconvData$mutCountReconstructed
    mutProcesses$Supplementary$processStab <- deconvData$processStab
    mutProcesses$Supplementary$processStabAvg <- deconvData$processStabAvg
  }
  options(warn = currentWarnings)
  return(mutProcesses)
}


deconvoluteMutCounts <- function (input_mutCounts, 
                                  params)
{
  # avoid NOTEs
  j <- NULL
  
  num_totIterations <- params$num_totIterations
  #perCore.iterations <- params$perCore.iterations
  num_processesToExtract <- params$num_processesToExtract
  distanceFunction <- params$distanceFunction
  thresh_removeWeakMutTypes <- params$thresh_removeWeakMutTypes
  num_parallelCores <- params$num_parallelCores
  guided <- params$guided
  debugStatus <- params$debug
  num_totReplicates <- params$num_totReplicates
  thresh_removeLastPercent <- params$thresh_removeLastPercent
  colnames(input_mutCounts) <- NULL
  rownames(input_mutCounts) <- NULL
  input_mutCounts <- base::as.matrix(input_mutCounts)
  bckgrnd.removed.mutCounts <- removeWeak(input_mutCounts, params)
  bckgrnd.removed.mutset <- bckgrnd.removed.mutCounts$removed.mutset
  bckgrnd.removed.mutCounts <- bckgrnd.removed.mutCounts$output.mutCounts
  total.mutationTypes <- nrow(bckgrnd.removed.mutCounts)
  total.samples <- ncol(bckgrnd.removed.mutCounts)
  
  # Start with a sample NMF iteration to guide the rest
  if (guided) {
    if (!debugStatus) {
      guide.W <- suppressMessages(extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                   params = params,
                                   bootStrap = FALSE))
    } else {
      guide.W <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                   params = params,
                                   bootStrap = FALSE)
    }
    guide.W <- guide.W$Wk
  } else {
    # dummy variable
    guide.W <- 0
  }
  
  if (num_parallelCores < 2) {
    muCounts.checkDF <- tryCatch(lapply(1:num_totIterations, (function(j){
      if (debugStatus) {
        if (j %in% base::as.integer(seq(1, num_totIterations, length.out = 100))) {
          message(paste("(",j,")", sep = ""), appendLF = FALSE)
        }
      }
      if (!debugStatus) {
        tmp.out <- suppressMessages(extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                                      bootStrap = TRUE,
                                                      params = params))
        
      } else {
        tmp.out <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                     bootStrap = TRUE,
                                     params = params)
        
      }
      if(guided){
        re.ORD <- rep(0, num_processesToExtract)
        for (ki in 1:num_processesToExtract) {
          my.i <- order(apply(abs(tmp.out$Wk - guide.W[,ki]), 2, sum))
          if (ki > 1) {
            my.i[re.ORD[1:(ki-1)]] <- max(my.i) + 1
          }
          re.ORD[ki] <- which.min(my.i)
        }
      } else {
        re.ORD <- 1:num_processesToExtract
      }
      # compare to ref and sort!
      if(num_processesToExtract > 1) {
        tmp.out$Wk <- tmp.out$Wk[,re.ORD]
        tmp.out$Hk <- tmp.out$Hk[re.ORD,]
      } else {
        tmp.out$Wk <- tmp.out$Wk
        tmp.out$Hk <- rbind(tmp.out$Hk)
      }
      tmp.out
    })), error = (function(e) {
      print(e)
    }))
    if (debugStatus) {
    message("Done!", appendLF = TRUE)
    }
  } else {
    # Initialize cores
    max.cores <- parallel::detectCores()
    max.cores <- max.cores - 1
    max.cores <- ifelse(max.cores < 1, 1, max.cores)
    use.cores <- ifelse(1 <= num_parallelCores & num_parallelCores <= max.cores,
                        num_parallelCores, max.cores)
    
    if (debugStatus) {
      cl <- suppressMessages(parallel::makeCluster(use.cores, outfile = ""))
    } else {
      cl <- suppressMessages(parallel::makeCluster(use.cores))
    }
    
    print(paste("Extracting", num_processesToExtract, "mutational signatures X",
                num_totIterations, "iterations using", use.cores, "cores"))
    suppressMessages(doParallel::registerDoParallel(cl))
    stuffToExp <- c("alexaNMF", "leadZeros", "extractSignatures", "frequencize",
                    "bootstrapCancerGenomes", 
                    #"guide.W", 
                    #"guided", 
                    #"params",
                    "chihJenNMF")
    
    suppressMessages(parallel::clusterExport(cl, stuffToExp))

    # Also, use guide.W to sort results
    muCounts.checkDF <- tryCatch(foreach::foreach(j = (1:num_totIterations),
                                                  .verbose = TRUE, .packages = "stats") %dopar%
                                                  {
                                                    if (j %in% base::as.integer(seq(1, num_totIterations, length.out = 100))) {
                                                      message(paste("(",j,")", sep = ""), appendLF = FALSE)
                                                    }
                                                    tmp.out <- extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                                                                                 params = params)
                                                    if(guided){
                                                      re.ORD <- rep(0, num_processesToExtract)
                                                      for (ki in 1:num_processesToExtract) {
                                                        my.i <- order(apply(abs(tmp.out$Wk - guide.W[,ki]), 2, sum))
                                                        if (ki > 1) {
                                                          my.i[re.ORD[1:(ki-1)]] <- max(my.i) + 1
                                                        }
                                                        re.ORD[ki] <- which.min(my.i)
                                                      }
                                                    } else {
                                                      re.ORD <- 1:num_processesToExtract
                                                    }
                                                    # compare to ref and sort!
                                                    if(num_processesToExtract > 1) {
                                                      tmp.out$Wk <- tmp.out$Wk[,re.ORD]
                                                      tmp.out$Hk <- tmp.out$Hk[re.ORD,]
                                                    } else {
                                                      tmp.out$Wk <- tmp.out$Wk
                                                      tmp.out$Hk <- rbind(tmp.out$Hk)
                                                    }
                                                    tmp.out
                                                  }, error = (function(e) {
                                                    print(e)
                                                  }), finally = (function(f) {
                                                    parallel::stopCluster(cl)
                                                  }))
    message("Done!", appendLF = TRUE)
  }
  
  # Run analysis in parallel
  W.all <- do.call(cbind, lapply(muCounts.checkDF, (function(tmp){
    tmp$Wk
  })))
  H.all <- do.call(rbind, lapply(muCounts.checkDF, (function(tmp){
    tmp$Hk
  })))
  errors.all <- lapply(muCounts.checkDF, (function(tmp){
    tmp$mutCounts.errors
  }))
  reconstruct.all <- lapply(muCounts.checkDF, (function(tmp){
    tmp$mutCounts.reconstructed
  }))
  
  # get rid of pre-processed data 
  # muCounts.checkDF <- NULL
  fltr.mutCounts.data <- filterOutIterations(wall = W.all,
                                             hall = H.all,
                                             cnt_errors = errors.all,
                                             cnt_reconstructed = reconstruct.all,
                                             params)
  stability.check <- evaluateStability(wall = fltr.mutCounts.data$Wall,
                                       hall = fltr.mutCounts.data$Hall, params)
  final.mutCounts.data <- addWeak(mutationTypesToAddSet = bckgrnd.removed.mutset,
                                  processes_I = stability.check$centroids, 
                                  processesStd_I = stability.check$centroidStd,
                                  Wall_I = fltr.mutCounts.data$Wall, 
                                  genomeErrors_I = fltr.mutCounts.data$mutCounts.errors,
                                  genomesReconstructed_I = fltr.mutCounts.data$mutCounts.reconstructed)
  
  # in case of 'freq' approach, adjust 
  deconvoluted.results <- list()
  deconvoluted.results$Wall <- final.mutCounts.data$Wall
  deconvoluted.results$Hall <- fltr.mutCounts.data$Hall
  deconvoluted.results$mutCountErrors <- final.mutCounts.data$mutCountErrors
  deconvoluted.results$mutCountReconstructed <- final.mutCounts.data$mutCountReconstructed
  deconvoluted.results$idx <- stability.check$idx
  deconvoluted.results$idxS <- stability.check$idxS
  deconvoluted.results$processes <- final.mutCounts.data$processes
  deconvoluted.results$processesStd <- final.mutCounts.data$processesStd
  deconvoluted.results$exposure <- stability.check$exposure
  deconvoluted.results$exposureStd <- stability.check$exposureStd
  deconvoluted.results$processStab <- stability.check$processStab
  deconvoluted.results$processStabAvg <- stability.check$processStabAvg
  deconvoluted.results$clusterCompactness <- stability.check$clusterCompactness
  #message("completed")
  return(deconvoluted.results)
}


extractXvarlinkData <- function(xvarLink_data) {
  tmpVars <- sub("^.*&var=[[:alnum:]]+(,.*)&.*$", "\\1", xvarLink_data)
  tmpVars[!grepl("^,.+", tmpVars)] <- NA
  tmpVars <- strsplit(tmpVars, ",")
  tmpVars <- do.call(rbind, lapply(1:length(tmpVars), function(i){
    if (tmpVars[[i]][1] == "" & length(tmpVars[[i]]) == 5) {
      tmpVars[[i]][2:5]
    } else {
      c(NA, NA, NA, NA)
    }
  }))
  out <- data.frame(chrXvar = base::as.character(tmpVars[,1]),
                    posXvar = base::as.numeric(tmpVars[,2]),
                    refXvar = base::as.character(tmpVars[,3]),
                    mutXvar = base::as.character(tmpVars[,4]),
                    stringsAsFactors = FALSE)
  return(out)
}


extractSignatures <- function (mutCountMatrix, 
                               params, 
                               bootStrap = TRUE)
{
  num_processesToExtract <- params$num_processesToExtract
  approach <- params$approach
  algorithm <- params$algorithm
  eps <- params$eps
  
  if (bootStrap) {
    bstrpd.result <- bootstrapCancerGenomes(mutCountMatrix) 
  } else {
    bstrpd.result <- mutCountMatrix
  }
  
  #if (approach != "counts") {
  #  frq.bstrpd <- frequencize(bstrpd.result)
  #  bstrpd.result <-  frq.bstrpd$freqs
  #}
  bstrpd.result[bstrpd.result < eps] <- eps
  
  if (algorithm %in% c("brunet", "alexa")) {
    nmf.results <- alexaNMF(v = bstrpd.result,
                            r = num_processesToExtract,
                            params = params)
  } else {
    nmf.results <- chihJenNMF(v =  bstrpd.result,
                              r = num_processesToExtract,
                              params = params)
  } 
  # nmf.results
  tmp.w <- nmf.results$w
  tmp.h <- nmf.results$h
  
  # This step seems useless for the new approach, as it divides or multiplies by 1
  # Keep for consistency and cause we let the user select what approach to use
  for (jj in 1:num_processesToExtract) {
    tmp.tot <- sum(tmp.w[, jj])
    tmp.w[, jj] <- tmp.w[, jj]/tmp.tot
    tmp.h[jj, ] <- tmp.h[jj, ] * tmp.tot
  }
  
  # modify for frequentized approach 
  #if (approach != "counts"){
  #  if (length(frq.bstrpd$colSums) == ncol(tmp.h)) {
  #    tmp.h <- sapply(1:length(frq.bstrpd$colSums), function(ai) {
  #      frq.bstrpd$colSums[ai] * tmp.h[,ai] / sum(tmp.h[,ai], na.rm = TRUE)
  #    })
  #  }  
  #}
  
  mutCountMatrix.reconstructed <- tmp.w %*% tmp.h
  
  result.list <- list()
  result.list$Wk <- tmp.w
  result.list$Hk <- tmp.h
  result.list$mutCounts.reconstructed <- mutCountMatrix.reconstructed
  result.list$mutCounts.errors <- bstrpd.result - mutCountMatrix.reconstructed
  
  if (params$logIterations != "lite") {
    result.list$inputMatrix <- bstrpd.result
    result.list$cosDist <- proxy::dist(rbind(base::as.vector(bstrpd.result),
                                             base::as.vector(mutCountMatrix.reconstructed)),
                                       method = "cosine")[1]
  }
  
  return(result.list)
}


filterSNV <- function(dataSet, 
                      seq_colNames) 
{
  if (!(is.data.frame(dataSet) | is.matrix(dataSet) ) |
      sum(!seq_colNames %in% colnames(dataSet)) > 0 |
      length(seq_colNames) < 2) {
    stop("Bad input or seq_colNames not found")
  }
  check.tab <- sapply(1:length(seq_colNames), (function(i){
    tmp <- gsub("[[:space:]]", "", dataSet[,seq_colNames[i]])
    toupper(tmp) %in% c("A","C","G","T")
  }))
  toKeep <- apply(check.tab, 1, (function(rw){
    sum(rw) == length(rw)
  }))
  out <- dataSet[toKeep,]
  rownames(out) <- NULL
  return(out)
}


frequencize <- function(countMatrix, 
                        permille = TRUE)
{
  out <- list()
  cf <- ifelse(permille, 1000, 1)
  
  out[["colSums"]] <- apply(countMatrix, 2, sum)
  out[["freqs"]] <- cf * apply(countMatrix, 2, (function(clmn){clmn/sum(clmn)}))
  return(out)
}


getCosmicSignatures <- function(forceUseMirror = FALSE, asMutSign = TRUE)
{
  mutType.labels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
                      "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
                      "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
                      "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                      "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
                      "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
                      "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
                      "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
                      "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  cosmic.url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  
  #initialize variable
  my_fullW <- NULL
  
  if (!forceUseMirror) {
    my_fullW <- tryCatch({TMP <- suppressWarnings(read.delim(cosmic.url, header = TRUE));
    rownames(TMP) <- TMP$Somatic.Mutation.Type;
    TMP <- TMP[,grep("Signature", colnames(TMP))];
    TMP},
    error = function(e) NULL)
  } 
  
  # private mirror
  if (is.null(my_fullW)) {
    private.mirror.url <- "http://www.labwizards.com/rlib/mutSignatures/cosmic.signatures.csv"
    my_fullW <- tryCatch({suppressWarnings(read.csv(private.mirror.url, 
                                                    header = TRUE, as.is = TRUE, row.names = 1))},
                         error = function(e2) { NULL })
  }
  if (is.null(my_fullW)) {
    message("An error occurred!")
    return(NULL)
  } else {
    if (sum(mutType.labels %in% rownames(my_fullW)) == length(mutType.labels)) {
      obj2rt <- my_fullW[mutType.labels,]
      names(obj2rt) <- gsub("Signature", "COSMIC", names(obj2rt), ignore.case = T)
      if(asMutSign)
        obj2rt <- mutSignatures::as.mutation.signatures(obj2rt)
      
      return(obj2rt)
    } else {
      message("An error occurred!")
      return(NULL)
    }   
  }
}


importVCFfiles <- function(vcfFiles, sampleNames = NULL){
  
  my.colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                   "FILTER", "INFO", "FORMAT", "XTR1", "XTR2", "XTR3")
  
  if (is.null(sampleNames) | length(sampleNames) != length(vcfFiles)){
    
    bypassNames <- paste("sample", 1:length(vcfFiles), sep = ".")
    
    sampleNames <- vcfFiles
    sampleNames <- sub("\\.vcf$", "", sub("^.*(\\\\|/)", "", tolower(sampleNames)))
    
    sampleNames[sampleNames == ""] <- bypassNames[sampleNames == ""]
    
  }
  out <- sapply(1:length(vcfFiles), (function(j){
    x <- vcfFiles[j]
    if (!file.exists(x)) {
      NULL
    } else {
      tmpVCF <- read.delim(x, comment.char = '#', header = F, stringsAsFactors = F)
      
      for (i in 1:ncol(tmpVCF)) {
        colnames(tmpVCF)[i] <- my.colnames[i]      
      }
      tmpVCF$SAMPLEID <- sampleNames[j]
      tmpVCF
    }
  }), simplify = FALSE, USE.NAMES = TRUE)
  out <- do.call(rbind, out)
  #names(out) <- sub("\\.vcf$", "", names(out))
  return(out)
}


plotMutCount <- function(mutCount, top = 50) {
  # avoid NOTEs
  count <- NULL
  feature <- NULL
  
  sampleLabs <- rownames(mutCount)
  rownames(mutCount) <- NULL
  
  nuSmplOrder <- order(apply(mutCount, 1, sum), decreasing = TRUE)
  mutCount <- mutCount[nuSmplOrder,]
  if (!is.null(sampleLabs))
    sampleLabs <- sampleLabs[nuSmplOrder]
  
  rownames(mutCount) <- NULL
  
  if (is.null(top)) {
    mutDF <- table2df(dataMatrix = mutCount)
  } else if(is.na(base::as.numeric(top[1]))){
    mutDF <- table2df(dataMatrix = mutCount)
  } else if (top[1] > nrow(mutCount) | top[1] < 2) {
    mutDF <- table2df(dataMatrix = mutCount)
  } else {
    mutDF <- table2df(dataMatrix = mutCount[1:top,])
    if (!is.null(sampleLabs))
      sampleLabs <- sampleLabs[1:top]
  }
  
  mutDF$sample <- 100000 + base::as.numeric(base::as.character(mutDF$sample))
  mutDF$sample <- base::as.character(mutDF$sample)
  
  tryCatch({mutDF$inputLabel <- do.call(c, lapply(sampleLabs, 
                                                  rep, times = ncol(mutCount)))}, 
           error = function(e) {
             message("damn")})
  
  mutDF$feature <- factor(mutDF$feature, levels = rev(colnames(mutCount)))
  bp <- ggplot(data=mutDF, aes(x=sample, y=count, fill=feature)) +
    geom_bar(stat="identity")
  bp <- bp + theme_minimal() + 
    theme(axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.line.y = element_line(colour = "black", size = 0.75),
          axis.line.x = element_line(colour = "black", size = 0.75),
          axis.ticks.y = element_line(colour = "black", size = 1),
          axis.ticks.length = unit(x = 6, "points"),
          plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(expand=c(0,0),
                       limits = c(0, 1.2 * max(apply(mutCount, 1, function(rx) sum(rx, na.rm = TRUE)))))

  return(bp)
}


plotMutTypeProfile <- function(mutCounts,
                               mutLabs,
                               freq = TRUE,
                               ylim = "auto",
                               ylab = "Fraction of Variants",
                               xlab = "Sequence Motifs",
                               xaxis_cex = 0.475,
                               cols = c("#4eb3d3", "#040404", "#b30000", "#bdbdbd", "#41ab5d", "#dd3497"),
                               main = "MutType Profile", 
                               ...)
{
  
  mutLabs <- toupper(mutLabs)
  if (!sum(regexpr("^(A|C|G|T)(\\[)(A|C|G|T)(>)(A|C|G|T)(\\])(A|C|G|T)$", toupper(mutLabs)) > 0) == length(mutLabs)) {
    message("mutTypes are in a non-standard format (ie, not Sanger 'A[G>A]T' style)")
    altPlot <- 1
  } else {
    altPlot <- 0
  }

  
  if (length(mutCounts) != length(mutLabs))
    stop("Mutation Type Labels and Mutation Counts do not match!") 
  
  mutCounts[is.na(mutCounts)] <- 0
  
  if (freq == TRUE) {
    if (ylim[1] == "auto")
      ylim <- c(0,0.2)
    freqs <- mutCounts/sum(mutCounts, na.rm = TRUE)
  } else {
    if (ylim[1] == "auto")
      ylim = c(0,(1.05*max(mutCounts)))
    freqs <- mutCounts
  }
  
  tinyMut <- gsub("\\](.*)$", "",
                  gsub("^(.*)\\[", "", toupper(mutLabs)))
  
  # Order the frequencies
  first.out <- lapply(sort(unique(tinyMut)), (function(mt){
    idxKeep <- which(tinyMut == mt)
    tmp.order <- order(mutLabs[idxKeep])
    out <- freqs[idxKeep][tmp.order]
    names(out) <- mutLabs[idxKeep][tmp.order]
    out
  }))
  names.first.out <- sort(unique(tinyMut))
  
  # Add a spacer
  second.out <- lapply(first.out, (function(vct){
    c(vct, 0)
  }))
  
  # Wrap together
  third.out <- do.call(c, second.out)
  third.out <- third.out[-length(third.out)]
  
  # Define color scheme
  col.out <- lapply(1:length(first.out), (function(ii){
    rep(cols[ii], (length(first.out[[ii]]) + 1))
  }))
  col.out <- do.call(c, col.out)
  col.out <- col.out[-length(col.out)]
  
  third.out[third.out>max(ylim)] <- max(ylim)
  third.shortLab <- gsub("\\[([[:alpha:]]>[[:alpha:]])\\]","-", names(third.out))
  
  xpos <- barplot(third.out,
                  col = col.out,
                  ylim = ylim,
                  axes = FALSE,
                  names.arg = FALSE,
                  ylab = ylab,
                  border = NA,
                  main = main)
  
  xpos.out <- lapply(1:length(first.out), (function(ii){
    c(rep(ii, length(first.out[[ii]])),0)
  }))
  xpos.out <- do.call(c, xpos.out)
  xpos.out <- xpos.out[-length(xpos.out)]
  
  #(max(ylim)*0.0075*0.2)
  box()
  axis(side = 1, tick = FALSE, hadj = 1, cex.axis = xaxis_cex, pos = (max(ylim) * 0.035), family = "mono",
       at = xpos,
       labels = third.shortLab,
       las = 2)
  
  my.y <- seq(min(ylim), max(ylim), length.out = 4)
  if (max(ylim) > 1.75) {
    my.y.labs <- format(round(my.y, digits = 0))
  } else {
    my.y.labs <- format(round(my.y, digits = 2))
  }
  
  axis(side = 2, tick = -0.005, at = my.y, labels = my.y.labs, las = 1, cex.axis = 0.75)
  
  zz <- unique(xpos.out)
  zz <- zz[zz>0]
  lab.xx <- sapply(zz, (function(i){
    mean(xpos[,1][xpos.out == i]) 
  }))
  mtext(names.first.out, side = 1, line = 1.5, at = lab.xx, cex = 0.8, col = cols)
  
  mtext(xlab, side = 1, line = 3.0, at = mean(xpos[,1]), cex = 1)
  # done! No return, just a plot!
}


prelimProcessAssess <- function(input,
                                maxProcess = 6,
                                approach = "counts",
                                plot = TRUE,
                                verbose = TRUE) 
{
  tmpParams <- setMutClusterParams(num_processesToExtract = 1, 
                                   approach = approach,
                                   stopconv = 10000,
                                   niter = 100000,
                                   thresh_removeWeakMutTypes = 0.01,
                                   thresh_removeLastPercent =0.025,
                                   num_totIterations = 10,
                                   num_parallelCores = 1,
                                   debug = FALSE,
                                   logIterations = "full")
  
  tmpParams <- coerceObj(x = tmpParams, to = "list")
  input_mutCounts <- coerceObj(x = input, to = "matrix")
  
  if (tmpParams$approach != "counts") {
    tmpFRQ.input <- frequencize(input_mutCounts)
    tmpParams$approach <- "counts"
    input_mutCounts <- tmpFRQ.input$freqs
  }
  bckgrnd.removed.mutCounts <- removeWeak(input_mutCounts,
                                          tmpParams)
  bckgrnd.removed.mutset <- bckgrnd.removed.mutCounts$removed.mutset
  bckgrnd.removed.mutCounts <- bckgrnd.removed.mutCounts$output.mutCounts
  total.mutationTypes <- nrow(bckgrnd.removed.mutCounts)
  total.samples <- ncol(bckgrnd.removed.mutCounts)
  tmpParams$num_totIterations <- 1
  
  # Calc initial (max) error
  medianMaxErr <- sum(t(sapply(1:nrow(bckgrnd.removed.mutCounts), (function(i){
    ((median(bckgrnd.removed.mutCounts[i,]) - bckgrnd.removed.mutCounts[i,]))^2
  }))))
  
  # Message
  if(verbose)
    message("Preliminary Mutational Process Assessment: ", appendLF = FALSE)
  
  outRes <- lapply(1:maxProcess, (function(i){
    tmpParams$num_processesToExtract <- i
    tmpRes <- suppressWarnings(suppressMessages(
      extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts,
                        params = tmpParams,
                        bootStrap = FALSE)))
    
    #apply(tmpRes$mutCounts.reconstructed, 2, sum)
    #mutCounts.reconstructed
    #tmpErr <- sum((tmpRes$mutCounts.errors)^2)
    tmpErr <- sum((bckgrnd.removed.mutCounts - tmpRes$mutCounts.reconstructed) ^ 2)
    
    if(verbose)
      message(".", appendLF = FALSE)
    tmpErr
  }))
  if(verbose)
    message("", appendLF = TRUE)
  
  err.points <- c(medianMaxErr, do.call(c, outRes))
  err.points <- (-1) * (err.points - max(err.points))
  err.points <- err.points / max(err.points, na.rm = TRUE)
  
  if (plot) {
    plot(err.points ~ c(0:(length(err.points) - 1)), type = "n", las = 1, axes = FALSE,
         ylab = "", xlab = "", main = "Preliminary Mutational Process Assessment")
    axis(side = 1, at = NULL, cex = 0.65, font = 3)
    axis(side = 2, at = seq(0, 1, by = 0.2), cex = 0.65, font = 3, 
         labels = format(seq(1, 0, by = -0.2), digits = 2, nsmall = 2), las = 1)
    lines(c(0:(length(err.points) - 1)), err.points, lwd = 1.5, col = "red2")
    points(c(0:(length(err.points) - 1)), err.points, pch = 19, col = "gray30")
    title(xlab="Num of Signatures", line=2.2, cex.lab=1.2, font = 2)
    title(ylab="Error (% vs. Median)", line=3.1, cex.lab=1.2, font = 2)
    box()
  }
  
  return(data.frame(numProcess = c(0:(length(err.points) - 1)),
                    percentErr = (1 - err.points),
                    stringsAsFactors = TRUE))
}


processVCFdata <- function(vcfData, 
                           BSGenomeDb,
                           chr_colName = "CHROM",
                           pos_colName = "POS",
                           ref_colName = "REF",
                           alt_colName = "ALT",
                           sample_colName = NULL,
                           nucl_contextN = 3,
                           verbose = TRUE)
{
  # make sure
  if (verbose)
    message("Processing VCF data: ", appendLF = FALSE)
  
  #Initialize
  reprepInput <- NULL
  
  # Fill reprepInput
  if (!is.null(sample_colName)){
    if (length(sample_colName) == 1 & sample_colName[1] %in% colnames(vcfData)){
      
      # compute how many unique samples
      unique.SID <- unique(vcfData[,sample_colName])
      reprepInput <- lapply(unique.SID, function(iID) {
        TMP <- vcfData[vcfData[,sample_colName] == iID, ]
        rownames(TMP)
        TMP
      })
    } else {
      stop("Bad input")
    }
  } else {
    #initialize input as list
    reprepInput <- list(vcfData)
  }
  
  # Triple check
  if(is.null(reprepInput))
    stop("Bad input")
  
  # Now loop
  bigOUT <- lapply(1:length(reprepInput), (function(jj){
    finalOut <- tryCatch({
      if(verbose)
        message(paste("[", jj, "/", length(reprepInput), "]", sep = ""), appendLF = F)
      
      tmpVCF <- suppressMessages({filterSNV(dataSet = reprepInput[[jj]],
                                            seq_colNames = c(ref_colName, alt_colName))})
      if (verbose)
        message(".", appendLF = FALSE)
      
      tmpVCF <- suppressMessages({attachContext(mutData = tmpVCF,
                                                chr_colName = chr_colName,
                                                start_colName = pos_colName,
                                                end_colName = pos_colName,
                                                nucl_contextN = 3,
                                                BSGenomeDb = BSGenomeDb)})
      if (verbose)
        message(".", appendLF = F)
      
      tmpVCF <- suppressMessages({removeMismatchMut(mutData = tmpVCF,
                                                    refMut_colName = ref_colName,
                                                    context_colName = "context",
                                                    refMut_format = "N")})
      if (verbose)
        message(".", appendLF = F)
      
      tmpVCF <- suppressMessages({attachMutType(mutData = tmpVCF,
                                                ref_colName = ref_colName,
                                                var_colName = alt_colName,
                                                var2_colName = alt_colName,
                                                context_colName = "context")})
      if (verbose)
        message(".", appendLF = F)
      
      tmpVCF
    }, error = function(e) NA)
  }))
  
  bigOUT <- do.call(rbind, bigOUT)
  
  if(verbose) {
    message("", appendLF = T)
    message("Done!", appendLF = T)
  }
  return(bigOUT)
}


removeMismatchMut <- function(mutData,
                              refMut_colName = "mutation",
                              context_colName = "context",
                              refMut_format = "N>N")
{
  # Validate input data and fields in mutData
  if(!((is.data.frame(mutData) | is.matrix(mutData) ) &
       sum(c(refMut_colName, context_colName) %in% colnames(mutData)) == 2))
    stop ("Issue with the input dataset. Make sure to feed in a data.frame or
          a matrix and double check the name of the fields pointing to chromosome
          name, start and end positions")
  if (! refMut_format %in% c("N>N", "N"))
    stop("Unrecognized refMut_format")
  
  output <- data.frame(mutData, stringsAsFactors = FALSE, row.names = NULL)
  colnames(output) <- colnames(mutData)
  if (refMut_format == "N>N") {
    output$mutSite.dnaSeq.Ref <- substr(mutData[,refMut_colName], 1, 1)
    output$mutSite.dnaSeq.Mut <- substr(mutData[,refMut_colName], 3, 3)
    refMut_colName <- "mutSite.dnaSeq.Ref"
  }
  
  pos.to.extr <- (0.5 * (nchar(output[,context_colName]) - 1)) + 1
  keep.id <- output[,refMut_colName] == substr(output[,context_colName], pos.to.extr, pos.to.extr)
  if (sum( keep.id == FALSE) > 0)
    message(paste("Removing", sum( keep.id == FALSE), "rows because of mismatch"))
  
  output <- output[keep.id, ]
  rownames(output) <- NULL
  return(output)
}


resolveMutSignatures <- function(mutCountData, 
                                 signFreqData, 
                                 byFreq = TRUE )
{
  # process expected objects
  if (class(mutCountData) == "mutationCounts")
    mutCountData <- coerceObj(x = mutCountData, to = "data.frame")
  
  if (class(signFreqData) == "mutationSignatures")
    signFreqData <- coerceObj(x = signFreqData, to = "data.frame")
  
  if (!(sum(!rownames(mutCountData) %in% rownames(signFreqData)) == 0 &
        sum(! rownames(signFreqData) %in% rownames(mutCountData)) == 0)) {
    stop (paste("There is an issue with the mutTypes.",
                "MutTypes in the mutType Count Matrix",
                "have to match those in the signature",
                "Matrix... check rownames()"))
  }
  mutCountData <- mutCountData[rownames(signFreqData),]
  
  if (byFreq) {
    full.Y <- apply(mutCountData, 2, (function(clmn){
      1000 * clmn / sum(clmn)
    }))
  } else {
    full.Y <- base::as.matrix(mutCountData)    
  }
  
  
  # record sum counts
  mutSums <- apply(mutCountData, 2, sum)
  
  # make sure we have frequencies
  my.signatures <- apply(signFreqData, 2, (function(clmn){
    clmn / sum(clmn)
  }))
  X <- base::as.matrix(my.signatures)
  
  # initialize results collector
  out <- list()
  
  res <- NMF::fcnnls(x = X, y = full.Y)
  
  beta.hat <- data.frame(t(res$x / ifelse (byFreq, 1000, 1)), stringsAsFactors = FALSE)
  
  colnames(beta.hat) <- colnames(signFreqData)
  rownames(beta.hat) <- colnames(mutCountData)
  out$results <- list()
  out$coeffs <- list()
  out$coeffs$beta.hat <- beta.hat
  out$coeffs$unexplained.mutNum <- round((1 - apply(beta.hat, 1, sum)) * mutSums, digits = 0)
  out$coeffs$unexplained.mutFrac <- out$coeffs$unexplained.mutNum  / mutSums
  
  if (byFreq) {
    for (i in 1:nrow(beta.hat)) {
      beta.hat[i,] <- beta.hat[i,] * mutSums[i]
    }
  }
  
  out$results$count.result <- mutSignatures::as.mutsign.exposures(beta.hat, samplesAsCols = FALSE)
  out$results$freq.result <- mutSignatures::as.mutsign.exposures(do.call(rbind, lapply(1:nrow(beta.hat), (function(jjj){
    beta.hat[jjj,] / sum(beta.hat[jjj,] )
  }))), samplesAsCols = FALSE)
  
  out$results$fitted <- res$fitted
  out$results$residuals <- res$residuals
  
  return(out)
}


revCompl <- function (DNAseq, 
                      rev = TRUE, 
                      compl = TRUE)
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


setMutClusterParams <- function(num_processesToExtract = 2, # number of de-novo signatures to extract
                                num_totIterations = 10,
                                num_parallelCores = 1,
                                thresh_removeWeakMutTypes = 0.01,
                                thresh_removeLastPercent = 0.07,
                                distanceFunction = "cosine",
                                num_totReplicates = 100,
                                eps = 2.2204e-16,
                                stopconv = 20000,
                                niter = 1000000,
                                guided = TRUE,
                                debug = FALSE,
                                approach = "freq",
                                stopRule = "DF",
                                algorithm = "brunet",
                                logIterations = "lite")
{
  #
  # Step-by-step Parameter Validation and preparation
  paramList <- list()
  #
  if (!((is.numeric(num_processesToExtract[1]) & num_processesToExtract[1] > 0) ))
    stop("Provide a reasonable number of signatures/processes to extract")
  paramList$num_processesToExtract <- round(num_processesToExtract[1])
  #
  if (!(is.numeric(num_totIterations[1]) & num_totIterations[1] > 0))
    stop("Provide a reasonable number of iterations to run (Bootstrapping)")
  paramList$num_totIterations <- round(num_totIterations[1])
  #
  if (!(is.numeric(num_parallelCores[1]) & num_parallelCores[1] > 0))
    stop("Provide a reasonable number of CPU cores to use for the analysis")
  paramList$num_parallelCores <- round(num_parallelCores[1])
  #
  #paramList$perCore.iterations <- as.integer(paramList$num_totIterations / paramList$num_parallelCores)
  #paramList$perCore.iterations <- ifelse(paramList$perCore.iterations < 1, 1, paramList$perCore.iterations)
  #
  if (!(is.numeric(thresh_removeWeakMutTypes[1]) & thresh_removeWeakMutTypes[1] >= 0 & thresh_removeWeakMutTypes[1] < 1))
    stop("Provide a reasonable (0.00-0.99) number of low-occurring mutation types to remove from the input before starting the analysis")
  paramList$thresh_removeWeakMutTypes <- thresh_removeWeakMutTypes[1]
  #
  if (!(is.numeric(thresh_removeLastPercent[1]) & thresh_removeLastPercent[1] >= 0 & thresh_removeLastPercent[1] < 1))
    stop("Provide a reasonable (0.00-0.99) number for filtering out poor iterations")
  paramList$thresh_removeLastPercent <- thresh_removeLastPercent[1]
  #
  allowed.dist.methods <- c("Braun-Blanquet", "Chi-squared", "correlation", "cosine", "Cramer", "Dice",
                            "eDice", "eJaccard", "Fager", "Faith", "Gower", "Hamman", "Jaccard",
                            "Kulczynski1", "Kulczynski2", "Michael", "Mountford", "Mozley", "Ochiai",
                            "Pearson", "Phi", "Phi-squared", "Russel", "simple matching", "Simpson",
                            "Stiles", "Tanimoto", "Tschuprow", "Yule", "Yule2", "Bhjattacharyya",
                            "Bray", "Canberra", "Chord", "divergence", "Euclidean", "fJaccard",
                            "Geodesic", "Hellinger", "Kullback", "Levenshtein", "Mahalanobis",
                            "Manhattan", "Minkowski", "Podani", "Soergel","supremum", "Wave", "Whittaker")
  if (!(is.character(distanceFunction[1]) & distanceFunction[1] %in% allowed.dist.methods))
    stop("Unknown method for calculating distances. For options, run: <<summary(proxy::pr_DB)>>")
  paramList$distanceFunction <- distanceFunction[1]
  #
  if (!(is.numeric(num_totReplicates[1]) & num_totReplicates[1] > 99))
    stop("Provide a reasonable number of replicates for stability evaluation of the results")
  paramList$num_totReplicates <- round(num_totReplicates[1])
  #
  if (!(is.numeric(eps[1]) & eps[1] > 0 & eps[1] < 0.0001))
    stop("Provide a reasonably small number (0 < n < 0.0001) for data overflow prevention")
  paramList$eps <- eps[1]
  #
  if (!(is.numeric(stopconv[1]) & stopconv[1] >= 500))
    stop("Provide a reasonable large number: number of 'conn-matrix-stable' iterations before stopping NMF")
  paramList$stopconv <- round(stopconv[1])
  #
  if (!(is.numeric(niter[1]) & niter[1] >= 20000))
    stop("Provide a reasonable large number: total NMF iterations")
  paramList$niter <- round(niter[1])
  #
  paramList$guided <- ifelse(guided, TRUE, FALSE)
  paramList$debug <- ifelse(debug, TRUE, FALSE)
  paramList$approach <- ifelse (approach == "counts", "counts", "freq")
  paramList$stopRule <- ifelse (stopRule == "LA", "LA", "DF")
  paramList$algorithm <- ifelse (tolower(algorithm) %in% c("brunet", "alexa"), "brunet", "chihjen")
  paramList$logIterations <- ifelse(tolower(logIterations) %in% c("lite", "light", "li"), "lite", "full")
  #
  return(new(Class = "mutFrameworkParams", params = paramList))
}


table2df <- function(dataMatrix, 
                     rowLab = "sample", 
                     colLab = "feature", 
                     valueLab = "count") 
{
  #
  if(!is.null(colnames(dataMatrix))) {
    names.X <- colnames(dataMatrix)
  } else {
    names.X <- 1:ncol(dataMatrix)
  }
  if(!is.null(rownames(dataMatrix))) {
    names.Y <- rownames(dataMatrix)
  } else {
    names.Y <- 1:nrow(dataMatrix)
  }
  #
  tmp <- do.call(rbind, lapply(1:nrow(dataMatrix), (function(i){
    t(sapply(1:ncol(dataMatrix), (function(j){
      c(names.Y[i], names.X[j], dataMatrix[i,j])
    })))
  })))
  tmp <- data.frame(tmp, stringsAsFactors = FALSE)
  if (is.numeric(dataMatrix[1,1]))
    tmp[,3] <- base::as.numeric(base::as.character(tmp[,3]))
  rownames(tmp) <- NULL
  colnames(tmp) <- c(rowLab, colLab, valueLab)
  return(tmp)
}


attachContext <- function(mutData,
                          BSGenomeDb,
                          chr_colName = "chr",
                          start_colName = "start_position",
                          end_colName = "end_position",
                          nucl_contextN = 3,
                          context_colName = "context")
{
  # Make sure about the right dependencies, or kill the process smoothly (no error)
  if (sum(grepl("^GenomicRanges", rownames(installed.packages()))) == 0 &
      sum(grepl("^BSgenome", rownames(installed.packages()))) == 0 ) {
    message("The `attachContext` function depends on the following bioconductor repositories:")
    message(" >  GenomicRanges")
    message(" >  a BSgenome data package, such as BSgenome.Hsapiens.UCSC.hg19")
    message("")
    message("Please install these packages via bioconductor and try again!")
    return(NULL)
  }
  
  # Validate input data and fields in mutData
  if(!((is.data.frame(mutData) | is.matrix(mutData) ) &
       sum(c(chr_colName, start_colName, end_colName) %in% colnames(mutData)) == 3 ))
    stop ("Issue with the input dataset. Make sure to feed in a data.frame or
          a matrix and double check the name of the fields pointing to chromosome
          name, start and end positions")
  if (!(is.numeric(nucl_contextN[1]) &
        nucl_contextN[1] > 2 & ((nucl_contextN[1] - 1) %% 2 == 0)))
    stop ("Please, specify a odd number > 2 as nucl_contextN")
  if (class(BSGenomeDb) != "BSgenome" | 
      sum(c("pkgname", "organism") %in% slotNames(BSGenomeDb)) != 2)
    stop ("Please, provide a valid 'BSgenome'-class object")
    
  # retrieve data, check quality, prepare
  output.data <- data.frame(mutData, stringsAsFactors = FALSE, row.names = NULL)
  rws.toExcl <- which(is.na(output.data[,chr_colName]) | 
                        is.na(output.data[,start_colName]) | 
                        is.na(output.data[,end_colName]))
  if (length(rws.toExcl) != 0) {
    output.data <- output.data[-rws.toExcl,]
  }
  chrPatt <- "(^chr([[:digit:]]{1,2})$)|(^chr(X|Y|M)$)|(^(X|Y|M)$)|(^[[:digit:]]{1,2}$)"
  rws.toKeep <- grep(chrPatt, output.data[,chr_colName])
  output.data <- output.data[rws.toKeep,]
  rws.toKeep <- grep("(^[[:digit:]]+$)", output.data[,start_colName])
  output.data <- output.data[rws.toKeep,]
  rws.toKeep <- grep("(^[[:digit:]]+$)", output.data[,end_colName])
  output.data <- output.data[rws.toKeep,]
  message(paste( (nrow(mutData) - nrow(output.data)), "rows were excluded from analysis"))
  
  tmp.ranges <- data.frame(cbind(base::as.character(gsub("^(C|c)hr", "", base::as.character(base::as.vector(output.data[,chr_colName])))),
                                 base::as.character(base::as.vector(output.data[,start_colName])),
                                 base::as.character(base::as.vector(output.data[,end_colName]))), 
                           stringsAsFactors = FALSE)
  colnames(tmp.ranges) <- c("chr", "start", "end")
  tmp.ranges$chr <- paste("chr", tmp.ranges$chr, sep = "")
  tmp.ranges$start <- base::as.numeric(tmp.ranges$start)
  tmp.ranges$end <- base::as.numeric(tmp.ranges$end)
  
  # Now let's check which are not point mutations
  non.point <- which(abs(tmp.ranges$end - tmp.ranges$start) > 1)
  if (length(non.point) > 0 ) {
    message(paste(length(non.point), "rows removed cause those are not point mutations"))
    tmp.ranges <- tmp.ranges[-non.point,]
    output.data <- output.data[-non.point,]
  }
  
  # Adding this line to fix an issue with human chromosomes X and Y, sometimes referred to as
  # chromosome 23 and 24
  if (BSGenomeDb@organism == "Homo sapiens"){
    tmp.ranges$chr[tmp.ranges$chr == "chr23"] <- "chrX"
    tmp.ranges$chr[tmp.ranges$chr == "chr24"] <- "chrY"
  } else if (BSGenomeDb@organism == "Mus musculus"){
    tmp.ranges$chr[tmp.ranges$chr == "chr20"] <- "chrX"
    tmp.ranges$chr[tmp.ranges$chr == "chr21"] <- "chrY"
  }
  
  # define a GRanges object
  half.context <- round((nucl_contextN[1] - 1) / 2, digits = 0)
  half.context <- ifelse(half.context <1, 1, half.context)
  tmp.granges <- GenomicRanges::GRanges(tmp.ranges$chr,
                                        IRanges::IRanges(start = (tmp.ranges$start - half.context),
                                                         end =  (tmp.ranges$start + half.context)),
                                        strand = "*")

  # make sure to remove seqs out-of-bounds
  message("Removing out-of-bounds positions...  ", appendLF = FALSE)
  tmp.seqNm <- BSgenome::seqnames(BSGenomeDb)
  tmp.seqNm <- tmp.seqNm[tmp.seqNm %in% unique(base::as.vector(BSgenome::seqnames(tmp.granges)))]
  chr.lens <- sapply(tmp.seqNm, (function(sqnm){
    length(BSGenomeDb[[sqnm]])
  }))
  keepId <- sapply(1:nrow(tmp.ranges), (function(jj){
    if (! tmp.ranges[jj, "chr"] %in% tmp.seqNm) {
      FALSE
    } else if ((tmp.ranges[jj, "end"] + half.context) > chr.lens[tmp.ranges[jj, "chr"]] |
               (tmp.ranges[jj, "start"] - 1) < 1) {
      FALSE
    } else {
      TRUE
    }
  }))
  message(sum(!keepId), " records were removed.", appendLF = TRUE)
  
  tmp.ranges <- tmp.ranges[keepId,]
  tmp.granges <- tmp.granges[keepId]
  output.data <- output.data[keepId,]
  
  all.contexts <- BSgenome::getSeq(BSGenomeDb, tmp.granges)
  all.contexts.seq <- sapply(1:length(all.contexts), (function(i){toString(all.contexts[[i]])}))
  
  # Attach sequences and wrap-up
  output.data[,context_colName] <- all.contexts.seq
  rownames(output.data) <- NULL
  message("Done!")
  return(output.data)
}


countMutTypes <- function (mutTable, mutType_colName = "mutType", sample_colName = NULL) 
{
  if (!(is.data.frame(mutTable) | is.matrix(mutTable))) 
    stop("Bad input")
  if (!mutType_colName %in% colnames(mutTable)) 
    stop("mutType field not found")
  if (!is.null(sample_colName)) {
    if (!sample_colName %in% colnames(mutTable)) 
      stop("sample_colName field not found")
  }
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
  custPatt01 <- "^(A|C|G|T)\\[(A|C|G|T)>(A|C|G|T)\\](A|C|G|T)$"
  if (sum(regexpr(custPatt01, mutTable[, mutType_colName]) > 
          0) != length(mutTable[, mutType_colName])){
    message("Non-standard format... MutTypes will be inferred.")
    
    mutType.labels <- sort(unique(mutTable[, mutType_colName]))
    
  } else {
    # fix reverse-complement if needed
    # only works for the standard three-nucletide system
    idx.to.fix <- which(!mutTable[, mutType_colName] %in% mutType.labels)
    if (length(idx.to.fix) > 0) {
      tmp <- mutTable[, mutType_colName][idx.to.fix]
      corrected.mutTypes <- sapply(1:length(tmp), (function(i) {
        out <- c(revCompl(substr(tmp[i], 7, 7)), "[", revCompl(substr(tmp[i], 
                                                                      3, 3)), ">", revCompl(substr(tmp[i], 5, 5)), 
                 "]", revCompl(substr(tmp[i], 1, 1)))
        paste(out, collapse = "", sep = "")
      }))
      mutTable[, mutType_colName][idx.to.fix] <- corrected.mutTypes
    }
  }
  
  if (sum(!mutTable[, mutType_colName] %in% mutType.labels) > 0) 
    stop("Problem with the mutType... Please check the input")
  if (is.null(sample_colName)) {
    my.mutTypes <- mutTable[, mutType_colName]
    out.1 <- sapply(mutType.labels, (function(mtt) {
      sum(my.mutTypes == mtt)
    }))
    out.1 <- data.frame(cbind(sample = out.1))
  } else {
    unique.cases <- unique(mutTable[, sample_colName])
    out.1 <- lapply(unique.cases, (function(csid) {
      tmp.df <- mutTable[mutTable[, sample_colName] == 
                           csid, ]
      out.3 <- sapply(mutType.labels, (function(mtt) {
        sum(tmp.df[, mutType_colName] == mtt)
      }))
      out.3 <- data.frame(cbind(out.3))
      colnames(out.3) <- csid
      out.3
    }))
    out.1 <- data.frame(do.call(cbind, out.1), stringsAsFactors = FALSE)
    tryCatch({
      colnames(out.1) <- unique.cases
    }, error = function(e) NULL)
  }
  fOUT <- new(Class = "mutationCounts", 
              x = out.1, 
              muts = data.frame(mutTypes = rownames(out.1),
                                stringsAsFactors = FALSE), 
              samples = data.frame(ID = colnames(out.1),
                                   stringsAsFactors = FALSE))
  return(fOUT)
}


sortByMutations <- function(x) {
  
  if (class(x) == "mutationSignatures") {
    outClass <- "mutationSignatures"
    x1 <- coerceObj(x = x, to = "data.frame")
  } else if (class(x) == "mutationCounts") {
    outClass <- "mutationCounts"
    x1 <- coerceObj(x = x, to = "data.frame")
  } else if (class(x) %in% c("matrix","data.frame")){
    outClass <- "data.frame"
    x1 <- x
  } else {
    stop ("Bad input")
  }
  
  # Get rownames == mutation names
  if (is.null(rownames(x1)))
    stop("Bad input")
  
  aRN <- rownames(x1)
  aRN <- sapply(sort(unique(sub("^.*\\[(.*)\\].*$", "\\1", aRN))), 
                function(mt) {
                  sort(aRN[grepl(mt, aRN)])
                }, USE.NAMES = FALSE, simplify = FALSE) 
  aRN <- do.call(c, aRN)
  ox <- x1[aRN,]
  
  # Format and return
  if (outClass == "data.frame") {
    return(ox)
  } else if (outClass == "mutationSignatures") {
    return(mutSignatures::as.mutation.signatures(ox))
  } else if (outClass == "mutationCounts") {
    return(mutSignatures::as.mutation.counts(ox))
  }
}


simplifySignatures <- function(x, asMutationSignatures = TRUE) {
  if(class(x) != "mutationSignatures")
    stop("Bad input")
  curMT <- getMutationTypes(x)
  if(sum(!grepl("[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]", curMT)) > 0)
    stop("Mutation Types are not compatible with simplifySignatures")
  
  # Start by defined tri-nucleotides
  mutType.labels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
                      "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
                      "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
                      "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                      "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
                      "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
                      "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
                      "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
                      "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  
  # Needs sorting
  input <- coerceObj(x = x, to = "data.frame")
  sigNames <- getSignatureIdentifiers(x)
  
  EXPLD <- sapply(mutType.labels, (function(pat){
    tpat <- sub("\\[", "\\\\[", 
                sub("\\]", "\\\\]", pat))
    input[grepl(tpat, rownames(input)),]
  }), simplify = FALSE, USE.NAMES = TRUE)
  
  OUT <- base::as.data.frame(t(sapply(EXPLD, function(xi) {
    base::as.numeric(apply(xi, 2, sum, na.rm = TRUE))
  })), row.names = mutType.labels)
  colnames(OUT) <- sigNames
  
  if (asMutationSignatures){
    return(mutSignatures::as.mutation.signatures(OUT))
  }
  
  altOut <- lapply(1:ncol(OUT), (function(j){
    curSGN <- base::as.numeric(OUT[,j])
    iSig <- data.frame(curSGN, row.names = mutType.labels)
    colnames(iSig) <- sigNames[j]
    
    AttDF <- do.call(rbind, sapply(mutType.labels, function(mmd) {
      zi <- EXPLD[[mmd]]
      nuNames <- gsub("[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]", "...", rownames(zi))
      nuNums <- base::as.data.frame(rbind(zi[,j]))
      names(nuNums) <- nuNames
      nuNums
    }, simplify = FALSE, USE.NAMES = TRUE))
    
    # out
    cbind(iSig, AttDF)
  }))
  names(altOut) <- sigNames
  return(altOut)
}


matchSignatures <- function(mutSign, 
                            reference = NULL, 
                            method = 'cosine',
                            threshold = 0.5,
                            plot = "TRUE")
{
  # avoid NOTEs
  newSign <- NULL
  refSign <- NULL
  
  if(is.null(reference)){
    my.ref <- getCosmicSignatures()
  } else if (class(reference) == "mutationSignatures") {
    my.ref <- reference
  } else {
    stop("reference: wrong data type. Provide NULL or a 'mutationSignatures'-class object")
  }
  if(class(mutSign) != "mutationSignatures") {
    stop("mutSign: wrong data type. Provide a 'mutationSignatures'-class object")
  }
  if (sum(!getMutationTypes(mutSign) %in% getMutationTypes(my.ref)) != 0 |
      sum(!getMutationTypes(my.ref) %in% getMutationTypes(mutSign)) != 0)
    stop("non-compatible mutation types")
  
  refSigNames <- getSignatureIdentifiers(my.ref)
  xSigNames <- getSignatureIdentifiers(mutSign)
  
  # Coerce to tabular objects
  my.ref2 <- my.ref@mutationFreq
  dimnames(my.ref2) <- list(my.ref@mutTypes[,1], 
                           my.ref@signatureId[,1])

  mutSign2 <- mutSign@mutationFreq
  dimnames(mutSign2) <- list(mutSign@mutTypes[,1], 
                             mutSign@signatureId[,1])
  
  my.ref <- my.ref2
  mutSign <- cbind(mutSign2[rownames(my.ref), ])
  
  #Debug
  #message("Debug - proceeded! :-)")
  
  distMat <- sapply(1:ncol(mutSign), function(i){
    TMP <- cbind(tmpSig=mutSign[,i], my.ref)
    TMPdist <- proxy::dist(TMP, method = method, by_rows = FALSE)
    TMPdist[1:ncol(my.ref)]
  })
  if (class(distMat) != "matrix") {
    distMat <- rbind(distMat)
  }
  
  dimnames(distMat) <- list(refSigNames, xSigNames)
  p <- NULL
  nuDF <- NULL
  
  if(plot) {
    
    nuDF <- do.call(rbind, lapply(1:nrow(distMat), function(i){
      do.call(rbind, lapply(1:ncol(distMat), function(j){
        data.frame(newSign = colnames(distMat)[j],
                   refSign = rownames(distMat)[i],
                   dist = distMat[i,j],
                   stringsAsFactors = FALSE)
      }))
    }))
    
    fillLims <- c(round(min(nuDF$dist)), round(max(nuDF$dist)))
    nuDF$dist[nuDF$dist >= fillLims[2]] <- fillLims[2]
    nuDF$dist[nuDF$dist <= fillLims[1]] <- fillLims[1]
    
    p <- ggplot(nuDF, aes(y=newSign, x=refSign))
    p <- p + geom_tile(aes(fill=dist), width=.875, height=.875)
    p <- p + scale_y_discrete(limits=rev(colnames(distMat)))
    p <- p + scale_x_discrete(limits=rownames(distMat))
    
    p <- p + theme_minimal(base_size = 11) + labs(x = "", y = "")
    p <- p + labs(title = "Signature Comparison")
    #p <- p + scale_fill_gradient2(low = "#f40000", mid = "#ffe9e9", high = "white",
    #                              midpoint = as.numeric(quantile(distMat, probs = 0.1, na.rm = TRUE)),
    #                              guide = "colourbar", 
    #                              name = paste(method, "\ndistance", sep = ""))
    #p <- p + scale_fill_gradient(low = "red2", high = "white",
    #                             guide = "colourbar", 
    #                             name = paste(method, "\ndistance", sep = ""))
    
    nuMid <- base::as.numeric(quantile(distMat, probs = 0.2, na.rm = TRUE))
    if (!is.null(threshold)) {
      if (!is.na(base::as.numeric(threshold[1]))){
        nuMid <- base::as.numeric(threshold[1])
      }  
    } 
    p <- p + scale_fill_gradient2(low = "#f40000", mid = "#ffe9e9", high = "white",
                                  midpoint = nuMid,
                                  guide = "colourbar", 
                                  name = paste(method, "\ndistance", sep = ""), 
                                  limits = fillLims)
    p <- p + theme(legend.position = c('right'), # position the legend in the upper left
                   legend.justification = 0, # anchor point for legend.position.
                   legend.text = element_text(size = 9, color = 'gray10'),
                   legend.title = element_text(size = 11),
                   plot.title = element_text(size = 13, face = 'bold', color = 'gray10', hjust = 0.5),
                   axis.text = element_text(face = 'bold'),
                   panel.grid.major.y = element_blank(),
                   panel.grid.major.x = element_blank(),
                   axis.text.x=element_text(angle = 90, hjust = 1 , vjust = 0.5,
                                            margin=margin(-5,0,-15,0)),
                   axis.text.y=element_text(hjust = 1, vjust = 0.5,
                                            margin = margin(0,3,0,0)))
    #print(p)
  }
  out <- list()
  out[["distanceMatrix"]] <- distMat
  out[["distanceDataFrame"]] <- nuDF
  out[["plot"]] <- p
  
  return(out)
  
}

#cosmix <- getCosmicSignatures()
#setwd("/media/dami/dami-data/backup/packages_CRAN/mutSignatures/ver_1.3.5/")
#package.skeleton(name = "mutSignatures", code_files = "core_mutSignatures_scr_6.R")

