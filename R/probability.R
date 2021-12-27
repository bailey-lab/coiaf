#------------------------------------------------
#' Dirichlet distribution
#'
#' Draw from a Dirichlet distribution given a vector of shape parameters. The
#' number of observations is inferred from the length of the shape vector.
#'
#' @param shape Vector of shape parameters.
#'
#' @seealso [rgamma()] for additional details.
#' @family distributions
#' @keywords internal
#' @examples
#' coiaf:::rdirichlet(c(1, 1, 1))
#' coiaf:::rdirichlet(c(1, 3, 2))
rdirichlet <- function(shape = rep(1, 3)) {
  # Check input
  assert_pos(shape)
  assert_vector(shape)

  # Draw from distribution
  x <- rgamma(length(shape), shape = shape, rate = max(shape))
  x / sum(x)
}

#------------------------------------------------
#' Beta-binomial distribution
#'
#' Draw from a Beta-binomial distribution.
#'
#' @param n Number of draws.
#' @param k Number of binomial trials.
#' @param alpha First shape parameter of beta distribution.
#' @param beta Second shape parameter of beta distribution.
#'
#' @seealso [rbeta()] and [rbinom()] for additional details.
#' @family distributions
#' @keywords internal
#' @examples
#' coiaf:::rbetabinom()
#' coiaf:::rbetabinom(n = 10, k = 10)
rbetabinom <- function(n = 1, k = 10, alpha = 1, beta = 1) {
  # Check inputs
  assert_single_pos_int(n)
  assert_pos_int(k)

  # Draw from distribution
  p <- rbeta(n = n, shape1 = alpha, shape2 = beta)
  rbinom(n = n, size = k, prob = p)
}
