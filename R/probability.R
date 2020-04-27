#------------------------------------------------
#' @title Draw from Dirichlet distribution
#'
#' @description Draw from Dirichlet distribution given a vector of shape
#'   parameters.
#'
#' @param shape Vector of shape parameters.
#'
#' @keywords internal

rdirichlet <- function(shape = rep(1,3)) {
  x <- rgamma(length(shape), shape = shape, rate = max(shape))
  return(x/sum(x))
}

#------------------------------------------------
#' @title Draw from Beta-binomial distribution
#'
#' @description Draw from Beta-binomial distribution.
#'
#' @param n Number of draws.
#' @param k Number of binomial trials.
#' @param alpha First shape parameter of beta distribution.
#' @param beta Second shape parameter of beta distribution.
#'
#' @keywords internal

rbetabinom <- function(n = 1, k = 10, alpha = 1, beta = 1) {
  p <- rbeta(n = n, shape1 = alpha, shape2 = beta)
  ret <- rbinom(n = n, size = k, prob = p)
  return(ret)
}
