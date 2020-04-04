#------------------------------------------------
#' @title Generate Simulated COI Curve
#'
#' @description Generate the simulated COI curve. In order to do this, we utilize
#'  the output of \code{\link{sim_biallelic}}, which created simulated data.
#'  We keep the PLAF, and compute whether a SNP is a variant or not, based
#'  on the simulated WSAF at that SNP -- acounting for potential sequencing
#'  error. To check whether our simulated WSAF correctly indicated a variant
#'  site or not, we determine whether a site should be varaint or not using
#'  the phased haplotype of the parasites.
#'
#' @param sim Output of \code{\link{sim_biallelic}}
#' @param seq_error The level of sequencing error that is assumed
#' @param cuts How often the data is summarized
#' @param theoretical_cois The theoretical COI curves to be examined
#' @return Simulated COI curve
#'
#' @export

simulated_coi <- function(sim, seq_error, cuts, theoretical_cois){
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

#------------------------------------------------
#' @title Compute COI
#'
#' @description Generate the COI of the sample by comparing the simualated
#' COI curve to the theoretical COI curve. To determine the actual COI value,
#' three different methods are utilized:
#' \describe{
#'   \item{\code{end}}{Determines the distance between the theoretical and
#'   simulated curve at a PLAF of 0.5. The COI is whichever theoretical COI curve
#'   is the closest to the simulated data..}
#'   \item{\code{ideal}}{Determines the distance between the theoretical and
#'   simulated curve at the ideal PLAF. The ideal PLAF is calculated by looking at
#'   the change between the COI of \eqn{i} and the COI of \eqn{i-1} and finding
#'   the PLAF for which this distance is maximal. The COI is whichever theoretical
#'   COI curve is the closet to the simulated data at the ideal PLAF.}
#'   \item{\code{end}}{Determines the distance between the theoretical and
#'   simulated curve at for all PLAFs. Computes the distance between the
#'   theoretical curves and the smulated curve. The COI is whichever theoretical
#'   curv has the smallest distance from the simulated curve.}
#'   }
#'
#' @param theory_cois The theoretical COI curves to be examined
#' @param sim_coi The simulated COI curve
#' @param method The method to be employed.
#' @param cuts How often the data is summarized
#' @return COI for the simulation
#'
#' @export

compute_coi <- function(theory_cois, sim_coi, method = "end", cuts){
  # Minus 1 because theory_cois now includes the PLAF at the end
  bound_coi = ncol(theory_cois) - 1

  if (method == "end"){
    ## Method 1: Compare last value
    # Remove first column because it contains COI of 1
    # (needed to compute ideal PLAF)
    last_theory <- theory_cois[, 2:bound_coi]

    # Get last row of theoretical COI curves and simulated date (PLAF of 0.5)
    last_theory <- last_theory[nrow(last_theory),]
    last_sim <- sim_coi$m_variant[nrow(sim_coi)]

    # Find coi by looking at minimum distance
    coi <- stringr::str_sub(names(which.min(abs(last_theory - last_sim))), -1)
  } else if (method == "ideal"){
    ## Method 2: Compute ideal PLAF
    # For each COI, find best PLAF and get theoretical and simulated values at that PLAF
    dist <- list()
    for (i in 2:bound_coi){
      # Want maximum point of the following
      diff = theory_cois[i] - theory_cois[i - 1]

      # Get max value and determine the PLAF
      ideal_PLAF <- theory_cois$PLAF[which.max(diff[[1]])]

      # Get max value and determine the theory WSAF
      theory_WSAF <- theory_cois[which.max(diff[[1]]), i]

      # Using max PLAF, determine which cut it would be part of
      # and then figure out m_variant value at this cut
      m_var <- sim_coi$m_variant[cut(ideal_PLAF, cuts)]

      # Find distance between theoretical and simulated curves
      dist[i-1] <- abs(theory_WSAF - m_var)
    }
    names(dist) <- colnames(theory_cois)[2:bound_coi]

    # Find coi by looking at minimum distance
    coi <- stringr::str_sub(names(which.min(dist)), -1)
  } else if (method == "overall"){
    # Find overlapping PLAFs for theory and simulated data
    match_theory_cois <- dplyr::filter(theory_cois, PLAF %in% sim_coi$midpoints)

    # Remove COI of 1 and PLAF
    match_theory_cois <- match_theory_cois[, 2:bound_coi]

    # Find absolute value of difference
    gap <- abs(colSums(match_theory_cois - sim_coi$m_variant))

    # Find coi by looking at minimum distance
    coi <- stringr::str_sub(names(which.min(gap)), -1)
  }

  ret_str <- paste("COI is", coi, sep = " ")
  return(ret_str)
}
