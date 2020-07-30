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
#'
#' @export

sim_biallelic <- function(COI = 3,
                          PLAF = runif(10, 0, 0.5),
                          coverage = 200,
                          alpha = 1,
                          overdispersion = 0,
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
  assert_bounded(epsilon)

  # Generate strain proportions
  w <- rdirichlet(rep(alpha, COI))

  # Generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(COI, 1, x), x = PLAF)
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
                                WSAF     = counts/coverage))
}

#------------------------------------------------
#' Process simulated data
#'
#' Generate the simulated COI curve. In order to do this, we
#' utilize the output of [sim_biallelic()], which created simulated
#' data. We keep the PLAF, and compute whether a SNP is a variant or not, based
#' on the simulated WSAF at that SNP -- accounting for potential sequencing
#' error. To check whether our simulated WSAF correctly indicated a variant
#' site or not, we determine whether a site should be variant or not using
#' the phased haplotype of the parasites.
#'
#' @param sim Output of [sim_biallelic()].
#' @param seq_error The level of sequencing error that is assumed.
#' @param cut How often the data is summarized.
#'
#' @return Simulated COI curve.
#'
#' @export

process_sim <- function(sim, seq_error = 0.01, cut = seq(0, 0.5, 0.01)) {
  # Check inputs
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)

  # Extract information from simulation
  df_sim <- data.frame(
    # PLAF
    plaf_cut = cut(sim$data$PLAF, cut, include.lowest = TRUE),

    # Determine if a site is a variant, accounting for sequencing error.
    variant = ifelse(sim$data$WSAF <= seq_error |
                       sim$data$WSAF >= (1 - seq_error), 0, 1)
  )

  # Case where COI is 1
  if (sim$COI == 1){
    df_sim$true_variant = df_sim$variant
  } else {
    df_sim$true_variant = as.integer(!apply(sim$phased, 2,
                                            function(x) {all(x) || all(!x)}))
  }

  # Average over intervals of PLAF
  df_sim_grouped <- df_sim %>%
    dplyr::group_by(.data$plaf_cut, .drop = FALSE) %>%
    dplyr::summarise(m_variant      = mean(.data$variant),
                     m_true_variant = mean(.data$true_variant),
                     bucket_size    = dplyr::n()) %>%
    as.data.frame()

  # Include midpoints
  df_sim_grouped$midpoints <- cut[-length(cut)] + diff(cut)/2
  df_sim_grouped <- stats::na.omit(df_sim_grouped)
}
