#------------------------------------------------
#' @title Draw from Dirichlet distribution
#'
#' @description Draw from Dirichlet distribution given a vector of shape
#'   parameters.
#'
#' @param shaper vector of shape parameters.
#'
#' @export

rdirichlet <- function(shaper = rep(1,3)) {
  x <- rgamma(length(shaper), shape = shaper, rate = max(shaper))
  return(x/sum(x))
}

#------------------------------------------------
#' @title Draw from Beta-binomial distribution
#'
#' @description Draw from Beta-binomial distribution.
#'
#' @param n number of draws.
#' @param k number of binomial trials.
#' @param alpha first shape parameter of beta distribution.
#' @param beta second shape parameter of beta distribution.
#'
#' @export

rbetabinom <- function(n = 1, k = 10, alpha = 1, beta = 1) {
  p <- rbeta(n = n, shape1 = alpha, shape2 = beta)
  ret <- rbinom(n = n, size = k, prob = p)
  return(ret)
}

