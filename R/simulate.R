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
#'    `coi`. The allele at each locus is drawn from a Bernoulli
#'    distribution with probability given by the `plaf`.
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
#' @param coi Complexity of infection.
#' @param plaf Vector of population-level allele frequencies at each locus.
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
#'   related to another. Default = 0 (unrelated). The implementation is similar
#'   to relatedness as defined in THE REAL McCOIL simulations. In the original
#'   paper (https://doi.org/10.1371/journal.pcbi.1005348) this is defined as:
#'   "... simulated relatedness (r) among lineages within the same host by
#'   sampling alleles either from an existing lineage within the same host
#'   (with probability r) or from the population (with probability (1-r))."
#'
#' @return A list of:
#' * `coi`: The COI used to simulate the data.
#' * `strain_proportions`: The strain proportion of each strain.
#' * `phased`: The phased haplotype.
#' * `data`: A dataframe of:
#'   + `plaf`: The population-level allele frequency.
#'   + `coverage`: The coverage at each locus.
#'   + `counts`: The count at each locus.
#'   + `wsaf`: The within-sample allele frequency.
#' * `inputs`: A dataframe of function input arguments:
#'   + `alpha`: Shape parameters of Dirichlet controlling strain proportions.
#'   + `overdispersion`: Overdispersion in count data.
#'   + `relatedness`: Within sample relatedness between strains.
#'   + `epsilon`: Probability of a single read being miscalled.
#'
#' @family simulated data functions
#' @export

sim_biallelic <- function(coi = 3,
                          plaf = runif(10, 0, 0.5),
                          coverage = 200,
                          alpha = 1,
                          overdispersion = 0,
                          relatedness = 0,
                          epsilon = 0) {

  # Check inputs
  assert_single_pos_int(coi)
  assert_vector(plaf)
  assert_bounded(plaf)

  # If a single value was input, then repeat coverage so that coverage is
  # applied over all loci.
  L <- length(plaf)
  if (length(coverage) == 1) coverage <- rep(coverage, L)

  # Continue to check inputs
  assert_vector(coverage)
  assert_pos_int(coverage)
  assert_same_length(plaf, coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_bounded(relatedness)
  assert_single_bounded(epsilon)

  # Generate strain proportions
  w <- rdirichlet(rep(alpha, coi))

  # Generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(coi, 1, x), x = plaf)

  ## Handle relatedness
  if (relatedness > 0 && coi > 1) {
    # If there is relatedness, we iteratively step through each lineage
    for (i in seq_len(coi - 1)) {
      # For each loci, we assume that it is related with probability relatedness
      rel_i <- as.logical(rbinom(L, size = 1, prob = relatedness))
      # And for those sites that are related, we draw from the other lineages
      if (any(rel_i)) {
        if (i == 1) {
          m[i + 1, rel_i] <- m[seq_len(i), rel_i]
        } else {
          m[i + 1, rel_i] <- apply(
            m[seq_len(i), rel_i, drop = FALSE],
            2,
            sample,
            size = 1
          )
        }
      }
    }
  }

  # Draw with the within sample allele frequencies (p_levels) are
  if (coi == 1) {
    p_levels <- m * w
  } else {
    p_levels <- colSums(sweep(m, 1, w, "*"))
  }

  # Rounding errors from multiplying w by m can cause numbers greater than 1
  p_levels[p_levels > 1] <- 1L

  # Add in genotyping error
  p_error <- p_levels * (1 - epsilon) + (1 - p_levels) * epsilon

  # Draw read counts, taking into account overdispersion
  if (overdispersion == 0) {
    counts <- rbinom(L, size = coverage, prob = p_error)
  } else {
    counts <- rbetabinom(
      L,
      k = coverage,
      alpha = p_error / overdispersion,
      beta = (1 - p_error) / overdispersion
    )
  }

  # Return list
  list(
    coi = coi,
    strain_proportions = w,
    phased = m,
    data = data.frame(
      plaf     = plaf,
      coverage = coverage,
      counts   = counts,
      wsaf     = counts / coverage
    ),
    inputs = data.frame(
      alpha          = alpha,
      overdispersion = overdispersion,
      relatedness    = relatedness,
      epsilon        = epsilon
    )
  )
}
