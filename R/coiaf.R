#' \pkg{coiaf} COI Estimation from Allele Frequency
#'
#' @docType package
#' @name coiaf
#'
#' @importFrom stats rbeta rbinom rgamma runif
#' @importFrom rlang .data
#'
"_PACKAGE"

# Silence the R CMD check notes on the where function
# This function comes from the tidyselect package but has not yet been exported
# causing R to through notes at the user.
utils::globalVariables("where")
