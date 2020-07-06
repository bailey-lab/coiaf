#------------------------------------------------
#' @title Compute COI
#'
#' @description Predict the COI of the sample by comparing the simulated
#' COI curve to the theoretical COI curve. To determine the predicted COI value,
#' three different methods are utilized:
#' \describe{
#'   \item{\code{end}}{Determines the distance between the theoretical and
#'   simulated curve at a PLAF of 0.5. The COI is whichever theoretical COI
#'   curve is the closest to the simulated data.}
#'   \item{\code{ideal}}{Determines the distance between the theoretical and
#'   simulated curve at the ideal PLAF. The ideal PLAF is calculated by looking
#'   at the change between the COI of \eqn{i} and the COI of \eqn{i-1} and
#'   finding the PLAF for which this distance is maximal. The COI is whichever
#'   theoretical COI curve is the closet to the simulated data at the ideal
#'   PLAF.}
#'   \item{\code{overall}}{Determines the distance between the theoretical and
#'   simulated curve at for all PLAFs. Computes the distance between the
#'   theoretical curves and the simulated curve. The COI is whichever
#'   theoretical curve has the smallest distance from the simulated curve.
#'   There is an option to choose one of several distance metrics:
#'   \itemize{
#'     \item{\code{abs_sum}:}{ Absolute value of sum of difference.}
#'     \item{\code{sum_abs}:}{ Sum of absolute difference.}
#'     \item{\code{squared}:}{ Sum of squared difference.}
#'     \item{\code{KL}:}{ Kullback-Leibler divergence.}
#'   }}
#'   }
#'
#' @param processed_data The simulated COI curve, which is the output of
#' \link{process_simulated_coi} or \link{process_real_data}.
#' @param theory_coi_range The range of COIs for which theoretical curves
#' will be calculated.
#' @param cut A vector indicating how often the data is summarized.
#' @param method The method to be employed. One of
#' \code{"end", "ideal", "overall"}.
#' @param dist_method The distance method used to determine the distance between
#' the theoretical and simulated curves for the "overall" method. One of
#' \code{"abs_sum", "sum_abs", "squared", "KL"}.
#' @param weighted An indicator indicating whether to compute weighted distance.
#'
#' @return A list of the following:
#' \describe{
#'   \item{\code{coi}}{The predicted COI for the simulation.}
#'   \item{\code{probability}}{A probability density function representing the
#'   probability of each COI.}
#'   }
#'
#' @export

compute_coi <- function(processed_data, theory_coi_range, cut,
                        method = "overall",
                        dist_method = "squared",
                        weighted = TRUE){
  ##Check inputs
  assert_pos_int(theory_coi_range, zero_allowed = FALSE)
  assert_vector(theory_coi_range)
  assert_increasing(theory_coi_range)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_single_string(method)
  assert_in(method, c("end", "ideal", "overall"))
  assert_single_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared", "KL"))
  assert_single_logical(weighted)

  # Calculate theoretical COI curves for the interval specified. Since we want
  # the theoretical curves and the simulated curves to have the PLAF values, we
  # compute the theoretical coi curves at processed_data$midpoints
  theory_cois <- theoretical_coi(theory_coi_range, processed_data$midpoints)

  # To check that PLAFs are the same
  assert_eq(theory_cois$plaf, processed_data$midpoints)

  # Minus 1 because theory_cois now includes the PLAF at the end
  bound_coi = ncol(theory_cois) - 1

  if (method == "end"){
    ## Method 1: Compare last value
    # Remove first column because it contains a COI we do not want tested
    # This COI is included to aid in the calculation of the `ideal` method
    # (needed to compute ideal PLAF)
    last_theory <- theory_cois[, 2:bound_coi]

    # Get last row of theoretical COI curves and simulated date (PLAF of 0.5)
    last_theory <- last_theory[nrow(last_theory),]
    last_sim <- processed_data$m_variant[nrow(processed_data)]

    # Find coi by looking at minimum distance
    dist <- abs(last_theory - last_sim)
    coi <- unlist(stringr::str_split(names(which.min(dist)), "_"))[2]
  } else if (method == "ideal"){
    ## Method 2: Compute ideal PLAF
    # For each COI, find best PLAF and get theoretical and simulated values at
    # that PLAF
    dist <- list()
    for (i in 2:bound_coi){
      # Want maximum point of the following
      diff = theory_cois[i] - theory_cois[i - 1]

      # Get max value and determine the PLAF
      ideal_PLAF <- theory_cois$plaf[which.max(diff[[1]])]

      # Get max value and determine the theory WSAF
      theory_WSAF <- theory_cois[which.max(diff[[1]]), i]

      # Using max PLAF, determine which cut it would be part of
      # and then figure out m_variant value at this cut
      m_var <- processed_data$m_variant[cut(ideal_PLAF, cut)]

      # Find distance between theoretical and simulated curves
      dist[i-1] <- abs(theory_WSAF - m_var)
    }
    # Name the distance vector so can extract COI information
    names(dist) <- colnames(theory_cois)[2:bound_coi]

    # Find coi by looking at minimum distance
    coi <- unlist(stringr::str_split(names(which.min(dist)), "_"))[2]
  } else if (method == "overall"){
    ## Method 3: Find distance between curves
    # Utilize helper function to compute overall distance between two curves
    overall_res <- distance_curves(processed_data, theory_cois, dist_method,
                                   weighted)

    # Extract information from the helper function
    coi  <- overall_res$coi
    dist <- overall_res$dist
  }

  # Distance to probability
  dist <- as.numeric(dist)
  dist <- 1 / (dist + 1e-5)
  dist <- dist / sum(dist)

  # Prepare list to return
  ret <- list(as.numeric(coi),
              dist)
  names(ret) <- c("coi", "probability")
  return(ret)
}

#------------------------------------------------
#' @title Compute distance between two curves
#'
#' @description Compute the distance between two curves using several methods.
#' \describe{
#'   \item{\code{abs_sum}}{Absolute value of sum of difference.}
#'   \item{\code{sum_abs}}{Sum of absolute difference.}
#'   \item{\code{squared}}{Sum of squared difference.}
#'   \item{\code{KL}}{Kullback-Leibler divergence.}
#'   }
#'
#' @inheritParams compute_coi
#' @param theory_cois The theoretical COI curves.
#'
#' @return A list of the following:
#' \describe{
#'   \item{\code{coi}}{The predicted COI for the simulation.}
#'   \item{\code{dist}}{The distance between each theoretical COI and the
#'   simulated COI curve.}
#'   }
#'
#' @keywords internal

distance_curves <- function(processed_data, theory_cois,
                            dist_method = "squared", weighted = TRUE){
  # Check inputs
  assert_single_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared", "KL"))
  assert_single_logical(weighted)

  # Find bound of COIs. Subtract 1 because theory_cois includes the PLAF at
  # the end
  bound_coi = ncol(theory_cois) - 1

  # Remove COI that is used to calculate distance method (1st one) and PLAF
  match_theory_cois <- theory_cois[, 1:bound_coi]

  # First find difference between theoretical and simulate curve. Weigh
  # difference if wanted
  gap <- match_theory_cois - processed_data$m_variant
  if (weighted){
    gap <- (gap * processed_data$bucket_size) / sum(processed_data$bucket_size)
  }

  if (dist_method == "abs_sum"){
    # Find sum of differences
    gap <- abs(colSums(gap))

  } else if (dist_method == "sum_abs"){
    # Find absolute value of differences
    gap <- colSums(abs(gap))

  } else if (dist_method == "squared"){
    # Squared distance
    gap <- colSums(gap ^ 2)

  } else if (dist_method == "KL"){
    # KL divergence
    gap <-  list()
    Q <- processed_data$m_variant
    Q <- Q/sum(Q)
    for (i in 1:ncol(match_theory_cois)){
      P <- match_theory_cois[,i]
      P <- P/sum(P)
      gap[i] <- suppressMessages(philentropy::KL(rbind(P, Q), unit = "log2"))
    }
    names(gap) <- colnames(theory_cois)[1:bound_coi]
    gap <- unlist(gap)
  }

  # Find coi by looking at minimum distance
  coi <- unlist(stringr::str_split(names(which.min(gap)), "_"))[2]

  # Prepare list to return
  ret <- list(coi  = coi,
              dist = gap)
  return(ret)
}
