#
##  ~~~~~~~~~~~~~~~~~~~~~~~~~
### ~~~~~ MutSignatures ~~~~~
##  ...all classes...
#   ~~~~~~~~~~~~~~~~~~~~~~~~~

# by Damiano Fantini, Ph.D.
# 2020-Mar-05
#
# This file includes class definition, as well
# as class initialization methods.
# Print and show methods for each class are also included.
# A total of 4 classes are defined.


#
# 1. mutFrameworkParams-class and related init/show methods
#

#' Class mutFrameworkParams.
#'
#' Class mutFrameworkParams defines objects including the set of 
#' parameters used for running a Mutational Signature Analysis.
#'
#' @slot params list including the set of parameters used for running a Mutational Signature Analysis
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#' 
#' @name mutFrameworkParams-class
#' @rdname mutFrameworkParams-class
#' @exportClass mutFrameworkParams
#' @export
setClass("mutFrameworkParams",
         slots = list(params = "list"))

#' Constructor method of the mutFrameworkParams Class.
#'
#' @rdname mutFrameworkParams-class
#' @importFrom methods callNextMethod
#' 
#' @param .Object the mutFrameworkParams object being built
#' @param params list including values for a set of mutFramework params
#' 
#' @aliases initialize,mutFrameworkParams-method
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
                               "logIterations", "seed")
            
            if(sum(names(params) %in% requiredNames) < 16)
              stop("Missing Params")
            
            .Object@params <- params
            .Object
          })

#' Show method of the mutFrameworkParams Class.
#'
#' @rdname mutFrameworkParams-show
#' 
#' @param object the mutFrameworkParams object being shown
#' 
#' @aliases show,mutFrameworkParams-method
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

#' Print method of the mutFrameworkParams Class.
#'
#' @rdname mutFrameworkParams-show
#' 
#' @param x the mutFrameworkParams object being printed
#' 
#' @importFrom methods show
#' @aliases print,mutFrameworkParams-method
setMethod("print", signature(x = "mutFrameworkParams"),
          function(x) {
            show(object = x)
          })



#
# 2. mutationSignatures-class and related init/show methods
#

#' Class mutationSignatures.
#'
#' Class mutationSignatures defines objects storing Mutational Signatures data.
#'
#' @slot mutationFreq data.frame including information about mutation frequencies
#' @slot mutTypes data.frame including information about mutation types
#' @slot signatureId data.frame including information about mutation signature Identifiers
#'
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#' 
#' @name mutationSignatures-class
#' @rdname mutationSignatures-class
#' @exportClass mutationSignatures
#' @export
setClass("mutationSignatures",
         slots = list(mutationFreq = "data.frame",
                      mutTypes = "data.frame",
                      signatureId = "data.frame"))

#' Constructor method of the mutationSignatures Class.
#'
#' @rdname mutationSignatures-class
#' @importFrom methods callNextMethod
#' 
#' @param .Object the mutationSignatures object being built
#' @param x data.frame including fequency data of multiple mutation signatures
#' @param muts data.frame including information about mutation types
#' @param signNames data.frame including information about mutation 
#' signature names (unique identifiers)
#' 
#' @aliases initialize,mutationSignatures-method
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

#' Show method of the mutationSignatures Class.
#'
#' @rdname mutationSignatures-show
#' 
#' @param object the mutationSignatures object being shown
#' 
#' @aliases show,mutationSignatures-method
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

#' Print method of the mutationSignatures Class.
#'
#' @rdname mutationSignatures-show
#' 
#' @param x the mutationSignatures object being printed
#' 
#' @importFrom methods show
#' @aliases print,mutationSignatures-method
setMethod("print", signature(x="mutationSignatures"),
          function(x) {
            show(x)
          })



#
# 3. mutationCounts-class and related init/show methods
#

#' Class mutationCounts.
#'
#' Class mutationCounts defines objects storing Mutation COunts data.
#'
#' @slot counts data.frame including information about mutation counts
#' @slot mutTypes data.frame including information about mutation types
#' @slot sampleId data.frame including information about sample identifiers
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#' 
#' @name mutationCounts-class
#' @rdname mutationCounts-class
#' @exportClass mutationCounts
#' @export
setClass("mutationCounts",
         slots = list(counts = "data.frame",
                      mutTypes = "data.frame",
                      sampleId = "data.frame"))

#' Constructor method of the mutationCounts Class.
#'
#' @rdname mutationCounts-class
#' @importFrom methods callNextMethod
#' 
#' @param .Object the mutationCounts object being built
#' @param x data.frame including mutation count values for each biological sample
#' @param muts data.frame including information about mutation types
#' @param samples data.frame including information about sample identifiers (unique names)
#' 
#' @aliases initialize,mutationCounts-method
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

#' Show method of the mutationCounts Class.
#'
#' @rdname mutationCounts-show
#' 
#' @param object the mutationCounts object being shown
#' 
#' @importFrom utils head
#' 
#' 
#' @aliases show,mutationCounts-method
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

#' Print method of the mutationCounts Class.
#'
#' @rdname mutationCounts-show
#' 
#' @param x the mutationCounts object being printed
#' 
#' @importFrom methods show
#' @aliases print,mutationCounts-method
setMethod("print", signature(x="mutationCounts"),
          function(x) {
            show(x)
          })




#
# 4. mutSignExposures-class and related init/show methods
#

#' Class mutSignExposures.
#'
#' Class mutSignExposures defines objects storing information about Exposures of 
#' biological samples to Mutational Signatures.
#'
#' @slot exposures data.frame including information about exposures
#' @slot sampleId data.frame including information about sample identifiers
#' @slot signatureId data.frame including information about signature identifiers
#' 
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#' 
#' @name mutSignExposures-class
#' @rdname mutSignExposures-class
#' @exportClass mutSignExposures
#' @export
setClass("mutSignExposures",
         slots = list(exposures = "data.frame",
                      sampleId = "data.frame",
                      signatureId = "data.frame"))

#' Constructor method of the mutSignExposures Class.
#'
#' @rdname mutSignExposures-class
#' @importFrom methods callNextMethod
#' 
#' @param .Object the mutSignExposures object being built
#' @param x data.frame including numeric values of exposures to mutational signatures
#' @param samples data.frame including information about biological sample identifiers (unique names)
#' @param signNames data.frame including information about mutational signature identifiers
#' 
#' @aliases initialize,mutSignExposures-method
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

#' Show method of the mutSignExposures Class.
#'
#' @rdname mutSignExposures-show
#' 
#' @param object the mutSignExposures object being shown
#' 
#' @aliases show,mutSignExposures-method
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

#' Print method of the mutSignExposures Class.
#'
#' @rdname mutSignExposures-show
#' 
#' @param x the mutSignExposures object being printed
#' 
#' @importFrom methods show
#' 
#' @aliases print,mutSignExposures-method
setMethod("print", signature(x="mutSignExposures"),
          function(x) {
            show(x)
          })


