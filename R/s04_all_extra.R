#
##  ~~~~~~~~~~~~~~~~~~~~~~~~~
### ~~~~~ MutSignatures ~~~~~
##  ...all functions...
#   ~~~~~~~~~~~~~~~~~~~~~~~~~

# Damiano Fantini, Ph.D.
# 2020-Mar-05

#' Decipher Mutational Signatures from Somatic Mutational Catalogs.
#'
#' Cancer cells accumulate DNA mutations as result of DNA damage and DNA repair pro-cesses.
#' mutSignatures is a computational framework that is aimed at deciphering DNA mutational
#' signatures oper-ating in cancer. The input is a numeric matrix of DNA mutation counts de-tected
#' in a panel of cancer samples. The framework performs Non-negative Matrix Factorization to extract
#' mutational signatures explaining the observed set of DNA mutations. The framework relies on
#' parallelization and is optimized for use on multi-core systems. This framework was described by Fantini D
#' et al (2020) \url{https://www.nature.com/articles/s41598-020-75062-0/} and is built upon a custom R-based
#' implementation of the original MATLAB WTSI frame-work by Alexandrov LB et al (2013)
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}. The mutSignatures framework has been described in
#' peer-reviewed publications, including Fantini D et al (2018) \url{https://www.nature.com/articles/s41388-017-0099-6/} and
#' Fantini D et al (2019) \url{https://www.sciencedirect.com/science/article/abs/pii/S1078143918303818/}. The framework includes three modules
#' that support raw data import and pre-processing, mutation counts deconvolution,
#' and data visualization.
#'
#' @references
#' More info, examples and vignettes:
#' \enumerate{
#'   \item{\strong{GitHub Repo}: \url{https://github.com/dami82/mutSignatures/} }
#'   \item{\strong{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/} }
#'   \item{\strong{2020 Sci Rep paper} describing the latest version of mutSignatures: \url{https://www.nature.com/articles/s41598-020-75062-0/} }
#'   \item{\strong{Oncogene paper}: Mutational Signatures operative in bladder cancer: \url{https://www.nature.com/articles/s41388-017-0099-6/} }
#'  }
#'
#' @docType package
#' @name mutSignatures-package
NULL




#' Input Data and Examples for Running Mutational Signatures Analyses
#'
#' A series of objects, including collections of DNA mutations from 50 Bladder cancer samples,
#' as well as mutational signatures extracted from the same samples. Mutation catalogs were obtained
#' from a TCGA bladder cancer dataset (data available from the BROAD Institute).
#' Original sample IDs were shuffled and then re-encoded.
#' Data are available in different formats, and can be used as input for running mutational signature analyses.
#'
#' @usage data("mutSigData")
#'
#' @format A list with 6 elements. Each element is a different type of \code{mutSignatures} input/data:
#' \describe{
#'   \item{inputA}{data.frame with 10401 rows and 4 columns. DNA mutation data mimicking a TCGA dataset downloaded using TCGAretriever/cBio}
#'   \item{inputB}{data.frame with 13523 rows and 12 columns. DNA mutation data mimicking a TCGA MAF file}
#'   \item{inputC}{data.frame with 13523 rows and 11 columns. DNA mutation data mimicking a VCF file decorated with a SAMPLEID column}
#'   \item{inputC.ctx}{data.frame with 13523 rows and 11 columns. DNA mutation data mimicking a VCF file decorated with a SAMPLEID column}
#'   \item{inputD}{data.frame with 13487 rows and 56 columns. DNA mutation data mimicking a set of VCF files casted into a 2D matrix (samples as columns)}
#'   \item{inputS}{list including data for silhouette plot generation (used in the vignette)}
#'   \item{blcaMUTS}{data.frame with 96 rows and 50 columns. A table of DNA mutation counts (rows are mutation types; columns are samples) }
#'   \item{blcaSIGS}{data.frame with 96 rows and 8 columns. Set of 8 mutational signatures (rows are mutation types; columns are signatures) }
#'   \item{.addON}{list of add-on functions (executed only upon request, not evaluated; these may require manual installation of external libraries from Bioconductor or GitHUB)}
#' }
#'
#' @source BLCA data were downloaded from \url{http://gdac.broadinstitute.org/} and then further processed, modified, and formatted.
#'
#' @details Examples and more information are available in the vignette, as well as at the following URL: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
#'
#' @examples
#' data(mutSigData)
#' print(mutSigData$input.A[1:6,])
#'
#' @name mutSigData
"mutSigData"




