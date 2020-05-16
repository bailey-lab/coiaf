#------------------------------------------------
#' @title Run a single COI test
#'
#' @description Runs a single full COI test.
#'
#' @param COI_range A number indicating the range of COIs to compare the
#' simulated data to.
#' @inheritParams sim_biallelic
#' @inheritParams simulated_coi
#' @param cut A vector indicating how often the data is summarized.
#' @inheritParams compute_coi
#'
#' @return Predicted COI value.
#'
#' @keywords internal

run_coi_test <- function(COI = 3,
                         COI_range = 10,
                         PLAF = runif(1000, 0, 0.5),
                         coverage = 100,
                         alpha = 1,
                         overdispersion = 0,
                         epsilon = 0,
                         seq_error = 0.01,
                         cut = seq(0, 0.5, 0.01),
                         method = c("end", "ideal", "overall"),
                         dist_method = c("abs_sum", "sum_abs", "squared", "KL"),
                         weighted = FALSE){

  # Check inputs
  assert_single_pos_int(COI)
  assert_single_pos_int(COI_range)
  assert_vector(PLAF)
  assert_bounded(PLAF, left = 0, right = 0.5)
  assert_single_pos_int(coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion)
  assert_single_bounded(epsilon)
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_single_string(method)
  assert_in(method, c("end", "ideal", "overall"))
  assert_single_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared", "KL"))
  assert_single_logical(weighted)

  # Simulate data
  sim <- sim_biallelic(COI, PLAF, coverage, alpha, overdispersion, epsilon)

  # Simulated data results
  sim_results <- simulated_coi(sim, seq_error, cut)

  # Determine ideal COI range (what true COIs to compare simulaiton to)
  if (COI <= COI_range){
    theory_cois_interval <-seq(1, COI + COI_range)
  } else {
    theory_cois_interval <- seq(COI - COI_range, COI + COI_range)
  }

  # Compute COI
  calc_coi <- compute_coi(theory_cois_interval, sim_results, cut, method,
                         dist_method, weighted)

  return (calc_coi$coi)
}


#------------------------------------------------
#' @title Test COIs
#'
#' @description Runs several iterations of a full COI test with varying
#' parameters
#'
#' @param repetitions The number of times each sample will be run.
#' @param COI_range A number indicating the range of COIs to compare the
#' simulated data to.
#' @inheritParams sim_biallelic
#' @inheritParams simulated_coi
#' @param cut A vector indicating how often the data is summarized.
#' @inheritParams compute_coi
#'
#' @return A list of the following dataframes:
#' \describe{
#'   \item{\code{predicted_coi}}{A dataframe of the predicted COIs. COIs are
#'   predicted using \link{compute_coi}. Each column represents a separate set
#'   of parameters. Each row represents a predicted COI. Predictions are done
#'   many times, depending on the value of \code{repetitions}.}
#'   \item{\code{param_grid}}{The parameter grid. The parameter grid is all
#'   possible combinations of the parameters inputted. Each row represents a
#'   unique combination.}
#'   \item{\code{error_bias}}{A dataframe containing any parameter that was
#'   varied and the associated mean absolute error and bias (mean error). By
#'   showing only parameters that were varied, the output is easier to interpret
#'   and does not have information about parameters that were held constant.}
#' }
#'
#' @export

coi_test <- function(repetitions = 10,
                     COI = 3,
                     COI_range = 5,
                     PLAF = runif(1000, 0, 0.5),
                     coverage = 100,
                     alpha = 1,
                     overdispersion = 0,
                     epsilon = 0,
                     seq_error = 0.01,
                     cut = seq(0, 0.5, 0.01),
                     method = c("end", "ideal", "overall"),
                     dist_method = c("abs_sum", "sum_abs", "squared", "KL"),
                     weighted = FALSE){

  # Check inputs
  assert_pos_int(repetitions)
  assert_pos_int(COI)
  assert_pos_int(COI_range)
  assert_vector(PLAF)
  assert_bounded(PLAF, left = 0, right = 0.5)
  assert_pos_int(coverage)
  assert_pos(alpha, zero_allowed = FALSE)
  assert_pos(overdispersion)
  assert_bounded(epsilon)
  assert_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_string(method)
  assert_in(method, c("end", "ideal", "overall"))
  assert_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared", "KL"))
  assert_logical(weighted)

  # Create parameter grid
  param_grid <- expand.grid(COI = COI,
                            COI_range = COI_range,
                            coverage = coverage,
                            alpha = alpha,
                            overdispersion = overdispersion,
                            epsilon = epsilon,
                            seq_error = seq_error,
                            method = method,
                            dist_method = dist_method,
                            weighted = weighted,
                            stringsAsFactors = FALSE)

  # Function to determine if pbapply is installed. If it is installed, it will
  # display a progress bar
  list_apply <- function(x, fun, ...){
    if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(x, fun, ...)
    } else {
      lapply(x, fun, ...)
    }
  }

  # Run each row of param_grid
  coi_pred <- list_apply(seq_len(nrow(param_grid)), function(x) {

    # Run each sample repetitions times
    repeats <- vapply(seq_len(repetitions), function(y) {
      test_result <-
        run_coi_test(param_grid$COI[x],
                     param_grid$COI_range[x],
                     PLAF,
                     param_grid$coverage[x],
                     param_grid$alpha[x],
                     param_grid$overdispersion[x],
                     param_grid$epsilon[x],
                     param_grid$seq_error[x],
                     cut,
                     param_grid$method[x],
                     param_grid$dist_method[x],
                     param_grid$weighted[x])
      return (test_result)
    }, FUN.VALUE = numeric(1))

    return(repeats)
  })

  # Determine how many unique COIs there are
  num_cois = length(COI)

  # Calculate the number of times each COI is repeated in param_grid
  num_repeat_cois = length(param_grid$COI) / num_cois

  # Name the predicted cois
  names(coi_pred) <- paste("coi",
                           param_grid$COI,
                           rep(seq(num_repeat_cois), each = num_cois),
                           sep="_")

  # Calculate mean absolute error
  len <- nrow(param_grid)
  coi_error <- lapply(seq_len(len), function(x) {
    sum(abs(coi_pred[[x]] - param_grid$COI[x])) / length(coi_pred[[x]])
    })

  # Calculate bias (mean error)
  coi_bias <- lapply(seq_len(len), function(x) {
    sum(coi_pred[[x]] - param_grid$COI[x]) / length(coi_pred[[x]])
  })

  # Save error and bias with changing parameters
  error_bias <- param_grid %>%
    dplyr::select_if(function(x) dplyr::n_distinct(x) > 1)
  error_bias$error <- unlist(coi_error)
  error_bias$bias  <- unlist(coi_bias)

  # Return predicted COIs and param_grid
  ret <- list(predicted_coi = as.data.frame(coi_pred),
              param_grid    = param_grid,
              error_bias    = error_bias)
  return (ret)
}

