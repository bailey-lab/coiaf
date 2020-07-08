#' \pkg{coiaf} COI Estimation from Allele Frequency
#'
#' @docType package
#' @name coiaf
#'
#' @importFrom stats rbeta rbinom rgamma runif
#' @importFrom rlang .data
#'
"_PACKAGE"

#' Example Real Data
#'
#' A small example dataset that contains WSAF from a sample of individuals.
#'
#' @format A list of matrices. The list indicates the region of the genome and
#' the VCF file where the samples came from. The two matrices have the same
#' structure. The rows of the matrix indicate the patient and the columns of
#' the matrix indicate the WSAF at each locus.
#'
#' @docType data
#' @keywords datasets
#' @name example_real_data
NULL

# Silence the R CMD check notes on the where function
# This function comes from the tidyselect package but has not yet been exported
# causing R to through notes at the user.
utils::globalVariables("where")
