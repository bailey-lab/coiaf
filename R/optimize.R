#------------------------------------------------
#' Likelihood of a COI
#'
#' A function to generate the likelihood of a specific COI value.
#'
#' The likelihood can be thought of the distance between two curves: the "real"
#' COI curve, generated from the inputted data, and the "simulated" COI curve,
#' which depends on the COI value specified. There are three different methods
#' implemented to compute the distance between two curves:
#' * `abs_sum`: Absolute value of sum of difference.
#' * `sum_abs`: Sum of absolute difference.
#' * `squared`: Sum of squared difference.
#'
#' @param coi The COI for which the likelihood will be generated.
#' @param processed_data The processed COI data. This is the output of
#' [process_sim()] or [process_real()].
#' @inheritParams compute_coi
#'
#' @return The likelihood for a specific COI value.
#' @family optimization functions
#' @export
likelihood <- function(coi,
                       processed_data,
                       distance = "squared",
                       coi_method = "variant") {
  # Check inputs
  assert_single_pos(coi)
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_string(coi_method)
  assert_in(coi_method, c("variant", "frequency"))

  # Compute theoretical curve
  if (coi_method == "variant") {
    theory_coi <- theoretical_coi(
      coi,
      processed_data$midpoints,
      coi_method = "variant"
    )
  } else {
    theory_coi <- theoretical_coi(
      coi,
      processed_data$midpoints,
      coi_method = "frequency"
    )
  }

  # Distance
  gap <- theory_coi - processed_data$m_variant

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

  # Gap is a named list with two entries: the coi and the PLMAF. We want to
  # return only the coi
  gap[1]
}



#------------------------------------------------
#' Optimize the COI
#'
#' A function to compute the COI of inputted data.
#'
#' The function utilizes [stats::optim()] In particular, the function utilizes
#' a quasi-Newton method to compute gradients and build a picture of the
#' surface to be optimized. The function uses a likelihood function as defined
#' by [likelihood()].
#'
#' @inheritParams compute_coi
#'
#' @return The predicted COI value.
#' @seealso [stats::optim()] for the complete documentation on the optimization
#' function.
#' @family optimization functions
#' @export
optimize_coi <- function(data,
                         data_type,
                         max_coi = 25,
                         seq_error = NULL,
                         bin_size = 20,
                         distance = "squared",
                         coi_method = "variant",
                         use_bins = FALSE) {

  # Check inputs
  assert_in(data_type, c("sim", "real"))
  assert_single_string(data_type)
  assert_single_pos_int(max_coi)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_string(coi_method)
  assert_in(coi_method, c("variant", "frequency"))
  assert_single_pos_int(bin_size)

  # removes NA from our data frame
  data <- remove_na_data(data, data_type)

  # Are we using bins or not
  if (!use_bins) {

    ret <- optimize_coi_regression(data,
                                  data_type,
                                  max_coi = max_coi,
                                  seq_error = seq_error,
                                  distance = distance,
                                  coi_method = coi_method,
                                  seq_error_bin_size = bin_size)
    return(ret)

  }

  # Warnings
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
    # If no coverage is provided, we assume that the coverage is uniform across
    # all loci
    if (!"coverage" %in% colnames(data)) {
      data$coverage <- rep(100, length(data$wsmaf))
    }

    processed <- process_real(
      data$wsmaf,
      data$plmaf,
      data$coverage,
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
      return(structure(
        1,
        notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method."
      ))
    }
  }

  # Compute COI
  # Details:
  #   par: The starting value of our parameter to be optimized
  #   fn: Our likelihood function. We pass in all variables for this function
  #   method: A modification of a quasi-Newton method that allows for bounds
  #   control:
  #     fnscale: Indicates that we want to minimize
  #     ndeps: The step sizes in the optimizer
  fit <- stats::optim(
    par = 2,
    fn = likelihood,
    processed_data = processed_data,
    distance = distance,
    coi_method = coi_method,
    method = "L-BFGS-B",
    lower = 1 + 1e-5,
    upper = max_coi,
    control = list(fnscale = 1, ndeps = 1e-5)
  )

  # Output warning if the model does not converge
  if (fit$convergence != 0) {
    if (fit$convergence == 1) {
      bullet <- "Iteration limit maxit has been reached."
    } else if (fit$convergence == 10) {
      bullet <- "Nelder-Mead simplex degeneracy."
    } else if (fit$convergence == 51) {
      bullet <- '"L-BFGS-B" method warning.'
    } else if (fit$convergence == 52) {
      bullet <- '"L-BFGS-B" method error'
    }
    message <- glue::glue(
      "The model did not converge:",
      "\n\u2716 {bullet}",
      "\n\u2716 Output of optim: {fit$message}."
    )
    warning(message, call. = FALSE)
  }

  # Return COI
  round(fit$par, 4)
}
