#------------------------------------------------
#' Likelihood of a COI
#'
#' A function to generate the likelihood of a specific COI value.
#' The likelihood can be thought of the distance between two curves: the "real"
#' COI curve, generated from the inputted data, and the "simulated" COI curve,
#' which depends on the COI value specified.
#'
#' @param coi The COI for which the likelihood will be generated.
#' @inheritParams compute_coi
#'
#' @return The likelihood for a specific COI value.
#'
#' @export

likelihood <- function(coi, processed_data,
                       dist_method = "squared", weighted = TRUE){

  # Check inputs
  assert_single_pos(coi)
  assert_single_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared"))
  assert_single_logical(weighted)

  # Compute theoretical curve
  theory_coi <- theoretical_coi(coi, processed_data$midpoints, method = "1")

  # Distance
  gap <- theory_coi - processed_data$m_variant
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
  }

  # Gap is a named list with two entries: the coi and the PLAF. We want to
  # return only the coi
  cio <- gap[1]
}

#------------------------------------------------
#' Optimize the COI
#'
#' A function to compute the COI of inputted data.
#' The function utilizes the [stats::optim()] function. In particular,
#' the function utilizes a quasi-Newton method to compute gradients and build a
#' picture of the surface to be optimized.
#'
#' @param data The data for which the COI will be computed.
#' @param data_type The type of the data to be analyzed. One of
#' `"sim"` or `"real"`.
#' @inheritParams coi_test
#'
#' @return The COI value.
#'
#' @seealso [stats::optim()]
#'
#' @export

optimize <- function(data,
                     data_type,
                     max_COI = 25,
                     seq_error = 0.01,
                     cut = seq(0, 0.5, 0.01),
                     dist_method = "squared",
                     weighted = TRUE) {

  # Check inputs
  assert_in(data_type, c("sim", "real"))
  assert_single_string(data_type)
  assert_single_pos_int(max_COI)
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_single_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared"))
  assert_logical(weighted)

  # Process data
  if (data_type == "sim"){
    processed_data <- process_sim(data, seq_error, cut)
  } else if (data_type == "real"){
    processed_data <- process_real(data$wsaf, data$plaf, seq_error, cut)
  }

  # Compute COI
  # Details:
  #   par: The starting value of our parameter to be optimized
  #   fn: Our likelihood function. We pass in all variables for this function
  #   method: A modification of a quasi-Newton method that allows for bounds
  #   control:
  #     fnscale: Indicates that we want to minimize
  #     ndeps: The step sizes in the optimizer
  fit <- stats::optim(par = 1,
                      fn = likelihood,
                      processed_data = processed_data,
                      dist_method = "squared",
                      weighted = TRUE,
                      method = "L-BFGS-B", lower = 1, upper = max_COI,
                      control = list(fnscale = 1, ndeps = 1e-5))

  # Output warning if the model does not converge
  if (fit$convergence != 0){
    message <- sprintf('The model did not converge. "%s"', fit$message) %>%
      stringr::str_squish() %>%
      stringr::str_wrap()
    warning(message, call. = FALSE)
  }

  # Return COI
  coi <- fit$par
}
