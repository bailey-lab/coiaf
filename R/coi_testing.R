#------------------------------------------------
#' @title Run a single COI test
#'
#' @description Runs a single full COI test.
#'
#' @param max_COI A number indicating the maximum COI to compare the
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
                         max_COI = 25,
                         PLAF = runif(1000, 0, 0.5),
                         coverage = 100,
                         alpha = 1,
                         overdispersion = 0,
                         epsilon = 0,
                         seq_error = 0.01,
                         cut = seq(0, 0.5, 0.01),
                         method = "overall",
                         dist_method ="squared",
                         weighted = TRUE){

  # Check inputs
  assert_single_pos_int(COI)
  assert_single_pos_int(max_COI)
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

  # Compute COI
  calc_coi <- compute_coi(theory_cois_interval = 1:max_COI,
                          sim_results, cut, method, dist_method, weighted)

  return (calc_coi$coi)
}


#------------------------------------------------
#' @title Test COIs
#'
#' @description Runs several iterations of a full COI test with varying
#' parameters
#'
#' @param repetitions The number of times each sample will be run.
#' @inheritParams run_coi_test
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
                     max_COI = 25,
                     PLAF = runif(1000, 0, 0.5),
                     coverage = 100,
                     alpha = 1,
                     overdispersion = 0,
                     epsilon = 0,
                     seq_error = 0.01,
                     cut = seq(0, 0.5, 0.01),
                     method = "overall",
                     dist_method = "squared",
                     weighted = TRUE){

  # Check inputs
  assert_pos_int(repetitions)
  assert_pos_int(COI)
  assert_single_pos_int(max_COI)
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
                            max_COI = max_COI,
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
                     param_grid$max_COI[x],
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

  # Name the predicted COIs
  names(coi_pred) <- paste("coi_", param_grid$COI, sep="")

  # Calculate mean absolute error
  len <- nrow(param_grid)
  boot_mae <- lapply(seq_len(len), function(x) {
    # Get the specific row (list of predicted COIs). Note that we convert to a
    # data frame because the boot package requires this.
    boot_data <- as.data.frame(coi_pred[x])

    # Create function to find mean absolute error for bootstrapping
    mae <- function(data, true_coi, indices){
      # Sample from the data
      sampled <- data[indices, ]

      # Compute mean absolute error for the sampled data and return
      sampled_mae <- sum(abs(sampled - true_coi)) / length(sampled)
      return(sampled_mae)
    }

    # Bootstrapping
    results <- boot::boot(data = boot_data, true_coi = param_grid$COI[x],
                          statistic = mae, R = 1000)

    # Get the normal confidence interval
    invisible(capture.output(CI <- boot::boot.ci(results, type = "norm")$norm))
    if (is.null(CI)){
      message(strwrap(sprintf('All values of the calculated statistic are the
      same. Could not calculate the bootstrapped confidence interval for a COI
                              of %s.', param_grid$COI[x])))
    }

    # Store the mean absolulte error and confidence interval bounds
    extract <- list(mae = results$t0, lower = CI[2], upper = CI[3])
    return(extract)
  })

  # Convert output of bootstrapping to dataframe structre
  boot_mae <- as.data.frame(do.call(rbind, boot_mae))

  # Calculate bias (mean error)
  coi_bias <- lapply(seq_len(len), function(x) {
    sum(coi_pred[[x]] - param_grid$COI[x]) / length(coi_pred[[x]])
  })

  # Save changing parameters with bootstrapping
  boot_error <- param_grid %>%
    dplyr::select_if(function(x) dplyr::n_distinct(x) > 1)
  boot_error <- dplyr::bind_cols(boot_error, boot_mae)
  boot_error$bias  <- unlist(coi_bias)

  # Return predicted COIs and param_grid
  ret <- list(predicted_coi = as.data.frame(coi_pred),
              param_grid    = param_grid,
              boot_error    = boot_error)
  return (ret)
}

