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
#' @inheritParams theoretical_coi
#' @inheritParams compute_coi
#'
#' @return The likelihood for a specific COI value.
#' @family optimization functions
#' @export

likelihood <- function(coi,
                       processed_data,
                       distance = "squared",
                       coi_method = "1") {
  # Check inputs
  assert_single_pos(coi)
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Compute theoretical curve
  if (coi_method == "1") {
    theory_coi <- theoretical_coi(coi,
                                  processed_data$midpoints,
                                  coi_method = "1")
  } else {
    theory_coi <- theoretical_coi(coi,
                                  processed_data$midpoints,
                                  coi_method = "2")
  }

  # Distance
  gap <- theory_coi - processed_data$m_variant

  if (distance == "abs_sum"){
    # Find sum of differences
    gap <- abs(colSums(gap))

  } else if (distance == "sum_abs"){
    # Find absolute value of differences
    gap <- colSums(abs(gap))

  } else if (distance == "squared"){
    # Squared distance
    gap <- colSums(gap ^ 2)
  }

  # Gap is a named list with two entries: the coi and the PLAF. We want to
  # return only the coi
  coi <- gap[1]
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
#' @param data The data for which the COI will be computed.
#' @param data_type The type of the data to be analyzed. One of
#' `"sim"` or `"real"`.
#' @inheritParams sensitivity
#' @inheritParams likelihood
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
                         coi_method = "1") {

  # Check inputs
  assert_in(data_type, c("sim", "real"))
  assert_single_string(data_type)
  assert_single_pos_int(max_coi)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Warnings
  if (distance != "squared") {
    message <- glue::glue("Please use the recommended distance metric:",
                          '\n\u2139 The recommended distance metric is "squared".',
                          '\n\u2716 User specified the "{distance}" metric.')
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
    processed <- process_real(data$wsaf, data$plaf, seq_error, bin_size, coi_method)
    processed_data <- processed$data
    seq_error <- processed$seq_error
    bin_size <- processed$bin_size
    cuts <- processed$cuts
  }

  # Special case where there are no heterozygous sites
  if (nrow(processed_data) == 0) return (coi <- 1)

  # Compute COI
  # Details:
  #   par: The starting value of our parameter to be optimized
  #   fn: Our likelihood function. We pass in all variables for this function
  #   method: A modification of a quasi-Newton method that allows for bounds
  #   control:
  #     fnscale: Indicates that we want to minimize
  #     ndeps: The step sizes in the optimizer
  fit <- stats::optim(par = 2,
                      fn = likelihood,
                      processed_data = processed_data,
                      distance = distance,
                      coi_method = coi_method,
                      method = "L-BFGS-B", lower = 1+1e-5, upper = max_coi,
                      control = list(fnscale = 1, ndeps = 1e-5))

  # Output warning if the model does not converge
  if (fit$convergence != 0){
    if (fit$convergence == 1) {
      bullet = "Iteration limit maxit has been reached."
    } else if (fit$convergence == 10) {
      bullet = "Nelder-Mead simplex degeneracy."
    } else if (fit$convergence == 51) {
      bullet = '"L-BFGS-B" method warning.'
    } else if (fit$convergence == 52) {
      bullet = '"L-BFGS-B" method error'
    }
    message <- glue::glue("The model did not converge:",
                          "\n\u2716 {bullet}",
                          "\n\u2716 Output of optim: {fit$message}.")
    warning(message, call. = FALSE)
  }

  # Return COI
  coi <- round(fit$par, 4)
}
