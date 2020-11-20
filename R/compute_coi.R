#------------------------------------------------
#' Predict the COI
#'
#' Predict the COI of the sample.
#'
#' \loadmathjax
#' Compare the within sample allele frequency (WSAF) and the population level
#' allele frequency (PLAF) of the sample to what a theoretical WSAF and PLAF
#' should look like. By examining the sample's WSAF and PLAF to the theoretical
#' WSAF and PLAF, an estimation can be made about what the COI of the sample is.
#' We refer to the sample's WSAF vs PLAF as the "sample curve" and refer to the
#' theoretical WSAF vs PLAF as the "theoretical curve." To determine the
#' predicted COI value, one of three different methods can be selected:
#' \describe{
#'   \item{`end`}{Determines the distance between the theoretical and
#'   sample curve at a PLAF of 0.5. The COI is whichever theoretical COI
#'   curve has the smallest distance to the simulated data.}
#'   \item{`ideal`}{Determines the distance between the theoretical and
#'   sample curve at the ideal PLAF. The ideal PLAF is calculated by looking
#'   at the change between the COI of \mjseqn{i} and the COI of \mjseqn{i-1} and
#'   finding the PLAF for which this distance is maximized. The COI is whichever
#'   theoretical COI curve has the smallest distance to the simulated data at
#'   the ideal PLAF.}
#'   \item{`overall`}{Determines the distance between the theoretical and
#'   simulated curve for all PLAFs. Computes the distance between the
#'   theoretical curves and the simulated curve. The COI is whichever
#'   theoretical curve has the smallest distance to the simulated curve.
#'   There is an option to choose one of several distance metrics:
#'   * `abs_sum`: Absolute value of sum of difference.
#'   * `sum_abs`: Sum of absolute difference.
#'   * `squared`: Sum of squared difference.
#'   }}
#'
#' @inheritParams optimize_coi
#' @param cut A vector indicating how often the data is summarized.
#' @param comparison The method to be employed. One of `"end"`, `"ideal"`,
#' `"overall"`.
#' @param distance The distance method used to determine the distance between
#' the theoretical and simulated curves for the `"overall"` method. One of
#' `"abs_sum"`, `"sum_abs"`, `"squared"`.
#' @param weighted An indicator indicating whether to compute the weighted
#' distance.
#' @inheritParams theoretical_coi
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
                        seq_error = 0.01,
                        cut = seq(0, 0.5, 0.01),
                        comparison = "overall",
                        distance = "squared",
                        weighted = TRUE,
                        coi_method = "1") {
  ##Check inputs
  assert_in(data_type, c("sim", "real"))
  assert_single_string(data_type)
  assert_single_pos_int(max_coi)
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_single_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_logical(weighted)
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Warnings
  if (comparison != "overall") {
    message <- glue::glue("Please use the recommended method:",
                          '\n\u2139 The recommended method is "overall".',
                          '\n\u2716 User specified the "{comparison}" method.')
    warning(message, call. = FALSE)
  }
  if (distance != "squared") {
    message <- glue::glue("Please use the recommended distance metric:",
                          '\n\u2139 The recommended distance metric is "squared".',
                          '\n\u2716 User specified the "{distance}" metric.')
    warning(message, call. = FALSE)
  }

  # Process data
  if (data_type == "sim") {
    processed_data <- process_sim(data, seq_error, cut, coi_method)
  } else if (data_type == "real") {
    processed_data <- process_real(data$wsaf, data$plaf,
                                   seq_error,
                                   cut,
                                   coi_method)
  }

  # Calculate theoretical COI curves for the interval specified. Since we want
  # the theoretical curves and the simulated curves to have the PLAF values, we
  # compute the theoretical COI curves at processed_data$midpoints
  theory_cois <- theoretical_coi(1:max_coi,
                                 processed_data$midpoints,
                                 coi_method)

  # To check that PLAFs are the same
  assert_eq(theory_cois$plaf, processed_data$midpoints)

  ## Special cases for Method 2 where COI = 1
  # If there is no heterozygous data, it means that the COI = 1. Otherwise, we
  # can compare the expected number of loci and the number of loci our
  # simulation gives us.
  if (coi_method == "2") {
    # No heterozygous data present
    if (dim(processed_data)[1] == 0) {
      ret <- list(coi = 1, probability = c(1, rep(0, max_coi - 1)))
      return(ret)
    }

    # Size is the number of loci per bucket
    if (data_type == "sim") {
      size_plaf <- data$data$plaf
      size_wsaf <- data$data$wsaf
    } else if (data_type == "real") {
      size_plaf <- data$plaf
      size_wsaf <- data$wsaf
    }
    size <- data.frame(plaf_cut = cut(size_plaf, cut, include.lowest = TRUE),
                       variant = size_wsaf) %>%
      dplyr::group_by(.data$plaf_cut, .drop = FALSE) %>%
      dplyr::summarise(bucket_size = dplyr::n()) %>%
      stats::na.omit()
    size$midpoints <- cut[-length(cut)] + diff(cut) / 2

    breaks = size$midpoints
    nloci  = size$bucket_size

    # 95% CI for how many heterozygous loci we expect per bucket
    CI <- Hmisc::binconf((2 * breaks * (1 - breaks)) * nloci, nloci) * nloci
    expectation <- tibble::tibble(cbind(size, CI))

    # If the number of loci in our simulated data is less than the expected
    # value, we predict that our COI will be 1
    combined <- dplyr::left_join(expectation, processed_data,
                                 by = c("plaf_cut", "midpoints"),
                                 suffix = c("_expect", "_data")) %>%
      tidyr::replace_na(list(bucket_size_data = 0))

    if (sum(combined$Lower - combined$bucket_size_data) >= 0) {
      ret <- list(coi = 1, probability = c(1, rep(0, max_coi - 1)))
      return(ret)
    }
  }

  # Minus 1 because theory_cois now includes the PLAF at the end
  bound_coi = ncol(theory_cois) - 1

  if (comparison == "end") {
    ## Method 1: Compare last value
    # Get last row of theoretical COI curves and simulated data (PLAF of 0.5)
    # Last column is removed because it contains the PLAF
    last_theory <- theory_cois[nrow(theory_cois), 1:bound_coi]
    last_sim    <- processed_data$m_variant[nrow(processed_data)]

    # Find COI by looking at minimum distance
    dist <- abs(last_theory - last_sim)
    coi  <- unlist(stringr::str_split(names(which.min(dist)), "_"))[2]

  } else if (comparison == "ideal") {
    ## Method 2: Compute ideal PLAF
    # For each COI, find best PLAF and get theoretical and simulated values at
    # that PLAF
    dist <- list()
    for (i in 1:bound_coi) {
      if (i != 1) {
        # Find difference between i and i-1 curve
        diff = theory_cois[i] - theory_cois[i - 1]

        # Get max value and determine the PLAF
        ideal_plaf <- theory_cois$plaf[which.max(diff[[1]])]

        # Get max value and determine the theory WSAF
        theory_wsaf <- theory_cois[which.max(diff[[1]]), i]

      } else {
        # Determine ideal PLAF and WSAF
        ideal_plaf  <- theory_cois$plaf[length(theory_cois$plaf)]
        theory_wsaf <- theory_cois[length(theory_cois$plaf), i]
      }

      # Using ideal PLAF, determine which cut it would be part of
      # and then figure out m_variant value at this cut
      m_var <- processed_data$m_variant[cut(ideal_plaf, cut)]

      # Find distance between theoretical and simulated curves
      dist[i] <- abs(theory_wsaf - m_var)
    }
    # Name the distance vector so can extract COI information
    names(dist) <- colnames(theory_cois)[1:bound_coi]

    # Find coi by looking at minimum distance
    coi <- unlist(stringr::str_split(names(which.min(dist)), "_"))[2]

  } else if (comparison == "overall") {
    ## Method 3: Find distance between curves
    # Utilize helper function to compute overall distance between two curves
    overall_res <- distance_curves(processed_data, theory_cois,
                                   distance, weighted)

    # Extract information from the helper function
    coi  <- overall_res$coi
    dist <- overall_res$dist
  }

  # Distance to probability
  dist <- as.numeric(dist)
  dist <- 1 / (dist + 1e-5)
  dist <- dist / sum(dist, na.rm = T)
  dist[is.nan(dist)] <- 0

  # Prepare list to return
  ret <- list(as.numeric(coi), dist)
  names(ret) <- c("coi", "probability")
  return(ret)
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

distance_curves <- function(processed_data, theory_cois,
                            distance = "squared", weighted = TRUE) {
  # Check inputs
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_logical(weighted)

  # Find bound of COIs. Subtract 1 because theory_cois includes the PLAF at
  # the end
  bound_coi = ncol(theory_cois) - 1

  # Remove COI that indicates the PLAF
  match_theory_cois <- theory_cois[, 1:bound_coi]

  # First find difference between theoretical and simulate curve. Weigh
  # difference if wanted
  gap <- match_theory_cois - processed_data$m_variant
  if (weighted) {
    gap <- gap * processed_data$bucket_size
  }

  if (distance == "abs_sum") {
    # Find sum of differences
    gap <- abs(colSums(gap))

  } else if (distance == "sum_abs") {
    # Find absolute value of differences
    gap <- colSums(abs(gap))

  } else if (distance == "squared") {
    # Squared distance
    gap <- colSums(gap ^ 2)
  }

  # Find COI by looking at minimum distance
  coi <- unlist(stringr::str_split(names(which.min(gap)), "_"))[2]

  # Prepare list to return
  ret <- list(coi  = coi,
              dist = gap)
}
