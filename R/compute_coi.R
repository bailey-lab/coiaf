#------------------------------------------------
#' Predict the COI
#'
#' Predict the COI of the sample.
#'
#' \loadmathjax
#' Compare the within sample allele frequency (WSMAF) and the population level
#' allele frequency (PLMAF) of the sample to what a theoretical WSMAF and PLMAF
#' should look like. By examining the sample's WSMAF and PLMAF to the theoretical
#' WSMAF and PLMAF, an estimation can be made about what the COI of the sample is.
#' We refer to the sample's WSMAF vs PLMAF as the "sample curve" and refer to the
#' theoretical WSMAF vs PLMAF as the "theoretical curve." To determine the
#' predicted COI value, one of three different methods can be selected:
#' \describe{
#'   \item{`end`}{Determines the distance between the theoretical and
#'   sample curve at a PLMAF of 0.5. The COI is whichever theoretical COI
#'   curve has the smallest distance to the simulated data.}
#'   \item{`ideal`}{Determines the distance between the theoretical and
#'   sample curve at the ideal PLMAF. The ideal PLMAF is calculated by looking
#'   at the change between the COI of \mjseqn{i} and the COI of \mjseqn{i-1} and
#'   finding the PLMAF for which this distance is maximized. The COI is whichever
#'   theoretical COI curve has the smallest distance to the simulated data at
#'   the ideal PLMAF.}
#'   \item{`overall`}{Determines the distance between the theoretical and
#'   simulated curve for all PLMAFs. Computes the distance between the
#'   theoretical curves and the simulated curve. The COI is whichever
#'   theoretical curve has the smallest distance to the simulated curve.
#'   There is an option to choose one of several distance metrics:
#'   * `abs_sum`: Absolute value of sum of difference.
#'   * `sum_abs`: Sum of absolute difference.
#'   * `squared`: Sum of squared difference.
#'   }}
#'
#' @param data The data for which the COI will be computed.
#' @param data_type The type of the data to be analyzed. One of
#' `"sim"` or `"real"`.
#' @param max_coi A number indicating the maximum COI to compare the
#' simulated data to.
#' @inheritParams process_real
#' @param comparison The method to be employed. One of `"end"`, `"ideal"`,
#' `"overall"`.
#' @param distance The distance method used to determine the distance between
#' the theoretical and simulated curves for the `"overall"` method. One of
#' `"abs_sum"`, `"sum_abs"`, `"squared"`.
#'
#' @return A list of the following:
#' * `coi`: The predicted COI of the sample.
#' * `probability`: A probability density function representing the probability
#'  of each COI.
#'
#' @export
compute_coi <- function(data,
                        data_type,
                        max_coi = 25,
                        seq_error = NULL,
                        bin_size = 20,
                        comparison = "overall",
                        distance = "squared",
                        coi_method = "variant") {
  ## Check inputs
  assert_in(data_type, c("sim", "real"))
  assert_single_string(data_type)
  assert_single_pos_int(max_coi)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)
  assert_single_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_string(coi_method)
  assert_in(coi_method, c("variant", "frequency"))

  # Warnings
  if (comparison != "overall") {
    message <- glue::glue(
      "Please use the recommended method:",
      '\n\u2139 The recommended method is "overall".',
      '\n\u2716 User specified the "{comparison}" method.'
    )
    warning(message, call. = FALSE)
  }
  if (distance != "squared") {
    message <- glue::glue(
      "Please use the recommended distance metric:",
      '\n\u2139 The recommended distance metric is "squared".',
      '\n\u2716 User specified the "{distance}" metric.'
    )
    warning(message, call. = FALSE)
  }

  # Process data
  if (data_type == "sim") {
    processed <- process_sim(data, seq_error, bin_size, coi_method)
    processed_data <- processed$data
    seq_error <- processed$seq_error
    bin_size <- processed$bin_size
    cuts <- processed$cuts
  } else if (data_type == "real") {
    processed <- process_real(
      data$wsmaf,
      data$plmaf,
      seq_error,
      bin_size,
      coi_method
    )
    processed_data <- processed$data
    seq_error <- processed$seq_error
    bin_size <- processed$bin_size
    cuts <- processed$cuts
  }

  # Special cases for the Frequency Method where COI = 1
  if (coi_method == "frequency") {
    check <- switch(data_type,
      "sim" = check_freq_method(data$data$wsmaf, data$data$plmaf, seq_error),
      "real" = check_freq_method(data$wsmaf, data$plmaf, seq_error)
    )

    # If the check returns FALSE, it means that the COI is likely 1
    if (!check) {
      ret <- list(
        coi = 1,
        probability = c(1, rep(0, max_coi - 1)),
        notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method."
      )
      return(ret)
    }
  }

  # Calculate theoretical COI curves for the interval specified. Since we want
  # the theoretical curves and the simulated curves to have the PLMAF values, we
  # compute the theoretical COI curves at processed_data$midpoints
  theory_cois <- theoretical_coi(
    1:max_coi,
    processed_data$midpoints,
    coi_method
  )

  # To check that PLMAFs are the same
  assert_eq(theory_cois$plmaf, processed_data$midpoints)

  # Minus 1 because theory_cois now includes the PLMAF at the end
  bound_coi <- ncol(theory_cois) - 1

  if (comparison == "end") {
    ## Method 1: Compare last value
    # Get last row of theoretical COI curves and simulated data (PLMAF of 0.5)
    # Last column is removed because it contains the PLMAF
    last_theory <- theory_cois[nrow(theory_cois), 1:bound_coi]
    last_sim <- processed_data$m_variant[nrow(processed_data)]

    # Find COI by looking at minimum distance
    dist <- abs(last_theory - last_sim)
    coi <- unlist(stringr::str_split(names(which.min(dist)), "_"))[2]
  } else if (comparison == "ideal") {
    ## Method 2: Compute ideal PLMAF
    # For each COI, find best PLMAF and get theoretical and simulated values at
    # that PLMAF
    dist <- list()
    for (i in 1:bound_coi) {
      if (i != 1) {
        # Find difference between i and i-1 curve
        diff <- theory_cois[i] - theory_cois[i - 1]

        # Get max value and determine the PLMAF
        ideal_plmaf <- theory_cois$plmaf[which.max(diff[[1]])]

        # Get max value and determine the theory WSMAF
        theory_wsmaf <- theory_cois[which.max(diff[[1]]), i]
      } else {
        # Determine ideal PLMAF and WSMAF
        ideal_plmaf <- theory_cois$plmaf[length(theory_cois$plmaf)]
        theory_wsmaf <- theory_cois[length(theory_cois$plmaf), i]
      }

      # Using ideal PLMAF, determine which cut it would be part of
      # and then figure out m_variant value at this cut
      m_var <- processed_data$m_variant[
        Hmisc::cut2(ideal_plmaf, cuts = processed$cuts, minmax = F)
      ]

      # Find distance between theoretical and simulated curves
      dist[i] <- abs(theory_wsmaf - m_var)
    }
    # Name the distance vector so can extract COI information
    names(dist) <- colnames(theory_cois)[1:bound_coi]

    # Find coi by looking at minimum distance
    coi <- unlist(stringr::str_split(names(which.min(dist)), "_"))[2]
  } else if (comparison == "overall") {
    ## Method 3: Find distance between curves
    # Utilize helper function to compute overall distance between two curves
    overall_res <- distance_curves(processed_data, theory_cois, distance)

    # Extract information from the helper function
    coi <- overall_res$coi
    dist <- overall_res$dist
  }

  # Distance to probability
  dist <- as.numeric(dist)
  dist <- 1 / (dist + 1e-5)
  dist <- dist / sum(dist, na.rm = T)
  dist[is.nan(dist)] <- 0

  # List to return
  list(coi = as.numeric(coi), probability = dist)
}

#------------------------------------------------
#' Compute distance between two curves
#'
#' Compute the distance between two curves using several methods.
#' \describe{
#'   \item{\code{abs_sum}}{Absolute value of sum of difference.}
#'   \item{\code{sum_abs}}{Sum of absolute difference.}
#'   \item{\code{squared}}{Sum of squared difference.}
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
distance_curves <- function(processed_data, theory_cois, distance = "squared") {
  # Check inputs
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))

  # Find bound of COIs. Subtract 1 because theory_cois includes the PLMAF at
  # the end
  bound_coi <- ncol(theory_cois) - 1

  # Remove COI that indicates the PLMAF
  match_theory_cois <- theory_cois[, 1:bound_coi]

  # Find the difference between the theoretical and simulated curves
  gap <- match_theory_cois - processed_data$m_variant

  # Weigh the buckets by the number of points in each bucket
  gap <- gap * processed_data$bucket_size

  if (distance == "abs_sum") {
    # Find sum of differences
    gap <- abs(colSums(gap))
  } else if (distance == "sum_abs") {
    # Find absolute value of differences
    gap <- colSums(abs(gap))
  } else if (distance == "squared") {
    # Squared distance
    gap <- colSums(gap^2)
  }

  # Find COI by looking at minimum distance
  coi <- unlist(stringr::str_split(names(which.min(gap)), "_"))[2]

  # Prepare list to return
  list(coi = coi, dist = gap)
}
