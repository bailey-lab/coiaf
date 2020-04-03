#------------------------------------------------
#' @title Compute COI
#'
#' @description Compute the COI of the sample. In order to do this, we utilize
#'  the output of \code{\link{sim_biallelic}}, which created simulated data.
#'  We keep the PLAF, and compute whether a SNP is a variant or not, based
#'  on the simulated WSAF at that SNP -- acounting for potential sequencing
#'  error. To check whether our simulated WSAF correctly indicated a variant
#'  site or not, we determine whether a site should be varaint or not using
#'  the phased haplotype of the parasites.
#'
#' @param sim Output of \code{\link{sim_biallelic}}.
#' @param seq_error The level of sequencing error that is assumed.
#' @param cuts How often the data is summarized.
#' @param theoretical_cois The theoretical COI curves to be examined.
#' @return Dataframe with averaged variant for several PLAF cuts.
#'
#' @export

compute_coi <- function(sim, seq_error, cuts, theoretical_cois){
  # Check inputs
  assert_single_bounded(seq_error)

  # Extract information from simulation
  df_sim <- data.frame(
    # PLAF
    PLAF = sim$data$PLAF,
    PLAF_cut = cut(sim$data$PLAF, cuts, include.lowest = TRUE),

    # Determine if a site is a variant, accounting for sequencing error.
    variant = ifelse(sim$data$WSAF < seq_error | sim$data$WSAF > (1 - seq_error), 0, 1),

    # True variant
    true_variant = as.integer(!apply(sim$phased, 2, function(x) {all(x) || all(!x)}))
  )

  # Average over intervals of PLAF
  df_sim_grouped <- dplyr::group_by(df_sim, PLAF_cut) %>%
    dplyr::summarise(m_variant = mean(variant), m_true_variant = mean(true_variant))
  df_sim_grouped <- as.data.frame(df_sim_grouped)

  # Include midpoints
  df_sim_grouped$midpoints <- cuts[-length(cuts)] + diff(cuts)/2

  return(df_sim_grouped)
}
