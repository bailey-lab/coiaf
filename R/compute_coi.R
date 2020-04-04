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
#' coi curve to the theoretical coi curve.
#'
#' @param theory_cois The theoretical COI curves to be examined
#' @param sim_coi The simulated coi curve
#' @param method The method to be employed.
#' @param cuts How often the data is summarized
#' @return COI for the simulation
#'
#' @export

compute_coi <- function(theory_cois, sim_coi, method = "end", cuts){
  # In order to calculate the best points to find COI, need
  # to have COI of one as well. Do COI[i] - COI[i-1].
  # So for method "end",  need to remove the first column

  # Minus 1 because theory_cois now includes the PLAF at the end
  bound_coi = ncol(theory_cois) - 1
  if (method == "end"){
    # Method 1: Compare last value
    last_theory <- theory_cois[, 2:bound_coi]
    last_theory <- last_theory[nrow(last_theory),]
    last_sim <- sim_coi$m_variant[nrow(sim_coi)]

    coi <- stringr::str_sub(names(which.min(abs(last_theory - last_sim))), -1)
  } else if (method == "ideal"){
    # Method 2: Compute ideal PLAF
    dist <- list()
    for (i in 2:bound_coi){
      # COI[i] - COI[i-1]
      diff = theory_cois[i] - theory_cois[i - 1]

      # Determine best PLAF value
      ideal_PLAF <- theory_cois$PLAF[which.max(diff[[1]])]

      theory_WSAF <- theory_cois[which.max(diff[[1]]), i]
      m_var <- sim_coi$m_variant[cut(ideal_PLAF, cuts)]

      dist[i-1] <- abs(theory_WSAF - m_var)
    }
    names(dist) <- colnames(theory_cois)[2:bound_coi]

      # df_grouped$PLAF_cut[cut(0.5, cut)]

      # # Determine what the WSAF value is at the best PLAF is
      # ideal[i-1] <- theory_cois[which.max(diff[[1]]), i]
      #
      # # For each ideal COI curve, find out which point from simulated
      # # curve is closet to that curve
      # dist_to_ideal[i-1] <- min(abs(df_grouped$m_variant - ideal[[i-1]]))
    # }
    # Now have the WSAFs in ideal place for each COI curve
    # names(ideal) <- colnames(theory_cois)[2:bound_coi]
    # names(dist_to_ideal) <- names(ideal)

    coi <- stringr::str_sub(names(which.min(dist)), -1)
  }

  ret_str <- paste("COI is", coi, sep = " ")
  return(ret_str)
}
