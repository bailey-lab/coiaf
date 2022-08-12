#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom glue glue
#' @import mathjaxr
#' @importFrom stats rbeta rbinom rgamma runif
#' @importFrom tibble tibble as_tibble as_tibble_col as_tibble_row
## usethis namespace: end
NULL

#' Example real data
#'
#' A small example dataset that contains within-sample allele frequencies
#' (WSMAFs) from a sample of individuals.
#'
#' @format A matrix of data. The rows of the matrix indicate the sample name and
#'   the columns of the matrix indicate the WSMAF at each locus.
#'
#' @source ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf6/
#'
#' @docType data
#' @keywords datasets
#' @name example_real_data
NULL

# Silence the R CMD check notes on the where function and the magrittr dot
# The where function comes from the tidyselect package but has not yet been
# exported causing R to throw notes at the user.
# The dot comes from magrittr and is used to pass data into the RHS of an
# expression
utils::globalVariables(c("where", "."))
