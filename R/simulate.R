#------------------------------------------------
#' Simulate biallelic data
#'
#' Simulate biallelic data from a simple statistical model. Inputs
#' include the complexity of infection (COI), population-level allele
#' frequencies (PLAF), and some parameters dictating skew and error
#' distributions. Outputs include the phased haplotypes and the unphased read
#' count and coverage data.
#'
#' \loadmathjax
#' Simulated data are drawn from a simple statistical model:
#' 1. Strain proportions are drawn from a symmetric Dirichlet
#'    distribution with shape parameter `alpha`.
#' 2. Phased haplotypes are drawn at every locus, one for each
#'    `COI`. The allele at each locus is drawn from a Bernoulli
#'    distribution with probability given by the `PLAF`.
#' 3. The "true" within-sample allele frequency at every locus is
#'    obtained by multiplying haplotypes by their strain proportions, and
#'    summing over haplotypes. Errors are introduced through the equation
#'    \mjsdeqn{wsaf_{error} = wsaf(1-e) + (1-wsaf)e}
#'    where \mjseqn{wsaf} is the WSAF without error and \mjseqn{e} is
#'    the error parameter `epsilon`.
#' 4. Final read counts are drawn from a beta-binomial distribution with
#'    expectation \mjseqn{w_{error}}. The raw number of draws is given by the
#'    `coverage`, and the skew of the distribution is given by the
#'    `overdispersion` parameter. If the `overdispersion` is equal to
#'    zero, then the distribution is binomial, rather than beta-binomial.
#'
#' @param COI Complexity of infection.
#' @param PLAF Vector of population-level allele frequencies at each locus.
#' @param coverage Coverage at each locus. If a single value then the same
#'   coverage is applied over all loci.
#' @param alpha Shape parameter of the symmetric Dirichlet prior on strain
#'   proportions.
#' @param overdispersion The extent to which counts are over-dispersed relative
#'   to the binomial distribution. Counts are Beta-binomially distributed, with
#'   the beta distribution having shape parameters
#'   \mjseqn{\frac{p}{overdispersion}} and
#'   \mjseqn{\frac{1-p}{overdispersion}}.
#' @param epsilon The probability of a single read being miscalled as the other
#'   allele. Applies in both directions.
#' @param relatedness The probability that a strain in mixed infections is
#'   related to another. Default = 0 (unrelated).
#'
#' @return A list of:
#' * `COI`: The COI used to simulate the data.
#' * `strain_proportions`: The strain proportion of each strain.
#' * `phased`: The phased haplotype.
#' * `data`: A dataframe of:
#'   + `PLAF`: Population-level allele frequency.
#'   + `coverage`: The coverage at each locus.
#'   + `counts`: The count at each locus.
#'   + `WSAF`: The within-sample allele frequency.
#' * `inputs`: A dataframe of function input arguments:
#'   + `alpha`: Shape parameters of Dirichlet controlling strain proportions.
#'   + `overdispersion`: Overdispersion in count data.
#'   + `relatedness`: Within sample relatedness between strains.
#'   + `epsilon`: Probability of a single read being miscalled.
#' *
#' @family simulated data functions
#' @export

sim_biallelic <- function(COI = 3,
                          PLAF = runif(10, 0, 0.5),
                          coverage = 200,
                          alpha = 1,
                          overdispersion = 0,
                          relatedness = 0,
                          epsilon = 0) {

  # Check inputs
  assert_single_pos_int(COI)
  assert_vector(PLAF)
  assert_bounded(PLAF)

  # If a single value was input, then repeat coverage so that coverage is
  # applied over all loci.
  L <- length(PLAF)
  if (length(coverage) == 1) {
    coverage <- rep(coverage, L)
  }

  # Continue to check inputs
  assert_vector(coverage)
  assert_pos_int(coverage)
  assert_same_length(PLAF, coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_pos(epsilon, zero_allowed = TRUE)
  assert_single_pos(relatedness, zero_allowed = TRUE)
  assert_bounded(epsilon)
  assert_bounded(relatedness)

  # Generate strain proportions
  w <- rdirichlet(rep(alpha, COI))

  # Generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(COI, 1, x), x = PLAF)

  # Handle relatedness
  if (relatedness > 0 && COI > 1) {

    # If there is relatedness we iteratively step through each lineage
    for (i in seq_len(COI - 1)) {

      # For each loci we assume that it is related with probability relatedness
      rel_i <- as.logical(rbinom(L, size = 1, prob = relatedness))

      # And for those sites that related, we draw the other lineages
      if (any(rel_i)) {
        m[i+1, rel_i] <- apply(m[seq_len(i), rel_i, drop = FALSE], 2, sample, size = 1)
      }

    }

  }

  # Draw with the within sample allele frequencies (p_levels) are
  if (COI == 1){
    p_levels = m*w
  } else{
    p_levels <- colSums(sweep(m, 1, w, "*"))
  }

  # Rounding errors from multiplying w by m can cause numbers greater than 1
  p_levels[p_levels > 1] <- 1L

  # Add in genotyping error
  p_error <- p_levels * (1 - epsilon) + (1 - p_levels) *epsilon

  # Draw read counts, taking into account overdispersion
  if (overdispersion == 0) {
    counts <- rbinom(L, size = coverage, prob = p_error)
  } else {
    counts <- rbetabinom(L,
                         k = coverage,
                         alpha = p_error/overdispersion,
                         beta = (1 - p_error)/overdispersion)
  }

  # Return list
  ret <- list(COI = COI,
              strain_proportions = w,
              phased = m,
              data = data.frame(PLAF     = PLAF,
                                coverage = coverage,
                                counts   = counts,
                                WSAF     = counts/coverage),
              inputs = data.frame(alpha = alpha,
                                  overdispersion = overdispersion,
                                  relatedness = relatedness,
                                  epsilon = epsilon))
}
