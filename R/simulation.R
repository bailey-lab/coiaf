#------------------------------------------------
#' @title Simulate biallelic data
#'
#' @description Simulate biallelic data from a simple statistical model. Inputs
#'   include the complexity of infection (COI), population-level allele
#'   frequencies (PLAF) and some parameters dictating skew and error
#'   distributions. Outputs include the phased haplotypes and the unphased read
#'   count and coverage data.
#'
#' @details Simulated data are drawn from a simple statistical model:
#'   \enumerate{
#'     \item Strain proportions are drawn from a symmetric Dirichlet
#'     distribution with shape parameter \code{alpha}.
#'     \item Phased haplotypes are drawn at every locus, one for each
#'     \code{COI}. The allele at each locus is drawn from a Bernoulli
#'     distribution with probability given by the \code{PLAF}.
#'     \item The "true" within-sample allele frequency at every locus is
#'     obtained by multiplying haplotypes by their strain proportions, and
#'     summing over haplotypes. Errors are introduced through the equation
#'     \deqn{wsaf_error = wsaf*(1-e) + (1-wsaf)*e}where \eqn{wsaf} is the WSAF
#'     without error and \eqn{e} is the error parameter \code{epsilon}.
#'     \item Final read counts are drawn from a beta-binomial distribution with
#'     expectation \eqn{w_error}. The raw number of draws is given by the
#'     \code{coverage}, and the skew of the distribution is given by the
#'     \code{overdispersion} parameter. If \code{overdispersion = 0} then the
#'     distribution is binomial, rather than beta-binomial.
#'   }
#'
#' @param COI complexity of infection.
#' @param PLAF vector of population-level allele frequencies at each locus.
#' @param coverage coverage at each locus. If a single value then the same
#'   coverage is applied over all loci.
#' @param alpha shape parameter of the symmetric Dirichlet prior on strain
#'   proportions.
#' @param overdispersion the extent to which counts are over-dispersed relative
#'   to the binomial distribution. Counts are Beta-binomially distributed, with
#'   the beta distribution having shape parameters \code{p/overdispersion} and
#'   \code{(1-p)/overdispersion}.
#' @param epsilon the probability of a single read being miscalled as the other
#'   allele. Applies in both directions.
#'
#' @export

sim_biallelic <- function(COI = 3,
                          PLAF = runif(10, 0, 0.5),
                          coverage = 100,
                          alpha = 1,
                          overdispersion = 0,
                          epsilon = 0) {

  # check inputs
  assert_single_pos_int(COI)
  assert_vector(PLAF)
  assert_bounded(PLAF)

  # if a single value was inputed, then repeat coverage so that coverage is
  # applied over all loci.
  L <- length(PLAF)
  if (length(coverage) == 1) {
    coverage <- rep(coverage, L)
  }

  # contine to check inputs
  assert_vector(coverage)
  assert_pos_int(coverage)
  assert_same_length(PLAF, coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_pos(epsilon, zero_allowed = TRUE)
  assert_bounded(epsilon)

  # generate strain proportions
  w <- rdirichlet(rep(alpha, COI))

  # generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(COI, 1, x), x = PLAF)
  p_levels <- colSums(sweep(m, 1, w, "*"))

  # rounding errors from multiplying w by m can cause numbers greater than 1
  p_levels[p_levels>1] <- 1L

  # add in genotyping error
  p_error <- p_levels*(1-epsilon) + (1-p_levels)*epsilon

  # draw read counts, taking into account overdispersion
  if (overdispersion == 0) {
    counts <- rbinom(L, size = coverage, prob = p_error)
  } else {
    counts <- rbetabinom(L,
                         k = coverage,
                         alpha = p_error/overdispersion,
                         beta = (1-p_error)/overdispersion)
  }

  # return list
  ret <- list(COI = COI,
              strain_proportions = w,
              phased = m,
              data = data.frame(PLAF = PLAF,
                                coverage = coverage,
                                counts = counts,
                                WSAF = counts/coverage))
  return(ret)
}
