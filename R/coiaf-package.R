#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom stats rbeta rbinom rgamma runif
#' @importFrom rlang abort warn inform
#' @importFrom glue glue
## usethis namespace: end
NULL

#' Example real data
#'
#' A small example dataset that contains within-sample allele frequencies
#' (WSAFs) from a sample of individuals.
#'
#' @format A matrix of data. The rows of the matrix indicate the sample name and
#'   the columns of the matrix indicate the WSAF at each locus.
#'
#' @source ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf6/
#'
#' @docType data
#' @keywords datasets
#' @name example_real_data
NULL

# Silence the R CMD check notes on the where function
# This function comes from the tidyselect package but has not yet been exported
# causing R to throw notes at the user.
utils::globalVariables("where")
