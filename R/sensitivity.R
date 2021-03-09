#------------------------------------------------
#' Single sensitivity analysis
#'
#' Runs a single full COI sensitivity analysis.
#'
#' @param max_coi A number indicating the maximum COI to compare the
#' simulated data to.
#' @inheritParams sim_biallelic
#' @inheritParams process_sim
#' @param cut A vector indicating how often the data is summarized.
#' @inheritParams compute_coi
#'
#' @return Predicted COI value.
#'
#' @keywords internal

single_sensitivity <- function(coi = 3,
                               max_coi = 25,
                               plaf = runif(1000, 0, 0.5),
                               coverage = 200,
                               alpha = 1,
                               overdispersion = 0,
                               relatedness = 0,
                               epsilon = 0,
                               seq_error = NULL,
                               bin_size = 20,
                               comparison = "overall",
                               distance ="squared",
                               weighted = TRUE,
                               coi_method = "1") {

  # Check inputs
  assert_single_pos_int(coi)
  assert_single_pos_int(max_coi)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_single_pos_int(coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion)
  assert_single_bounded(relatedness)
  assert_single_bounded(epsilon)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)
  assert_single_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_logical(weighted)
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Simulate data
  sim_data <- sim_biallelic(coi, plaf, coverage, alpha, overdispersion,
                            relatedness, epsilon)

  # Compute COI
  calc_coi <- compute_coi(sim_data, "sim", max_coi, seq_error, bin_size,
                          comparison, distance, coi_method)
}


#------------------------------------------------
#' Sensitivity analysis
#'
#' Runs several iterations of a full COI sensitivity analysis with varying
#' parameters.
#'
#' @param repetitions The number of times each sample will be run.
#' @inheritParams single_sensitivity
#'
#' @return A list of the following:
#' * `predicted_coi`: A dataframe of the predicted COIs. COIs are
#'   predicted using [compute_coi()]. Each column represents a separate set
#'   of parameters. Each row represents a predicted COI. Predictions are done
#'   many times, depending on the value of `repetitions`.
#' * `probability`:A list of matrices containing the probability
#'   that our model predicted each COI value. Each row contains the probability
#'   for a different run. The first row contains the average probabilities over
#'   all the runs.
#' * `param_grid`: The parameter grid. The parameter grid is all
#'   possible combinations of the parameters inputted. Each row represents a
#'   unique combination.
#' * `boot_error`: A dataframe containing information about the error
#'   of the algorithm. The first column indicates the COI that was fed into the
#'   simulation. The other columns indicate the mean absolute error (mae),
#'   the lower and upper bounds of the 95% confidence interval and the bias.
#'
#' @export

sensitivity <- function(repetitions = 10,
                        coi = 3,
                        max_coi = 25,
                        plaf = runif(1000, 0, 0.5),
                        coverage = 200,
                        alpha = 1,
                        overdispersion = 0,
                        relatedness = 0,
                        epsilon = 0,
                        seq_error = NULL,
                        bin_size = 20,
                        comparison = "overall",
                        distance = "squared",
                        weighted = TRUE,
                        coi_method = "1") {

  # Check inputs
  assert_pos_int(repetitions)
  assert_pos_int(coi)
  assert_single_pos_int(max_coi)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_pos_int(coverage)
  assert_pos(alpha, zero_allowed = FALSE)
  assert_pos(overdispersion)
  assert_bounded(relatedness)
  assert_bounded(epsilon)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)
  assert_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_logical(weighted)
  assert_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Create parameter grid
  param_grid <- expand.grid(coi = coi,
                            max_coi = max_coi,
                            coverage = coverage,
                            alpha = alpha,
                            overdispersion = overdispersion,
                            relatedness = relatedness,
                            epsilon = epsilon,
                            seq_error = seq_error,
                            bin_size = bin_size,
                            comparison = comparison,
                            distance = distance,
                            weighted = weighted,
                            coi_method = coi_method,
                            stringsAsFactors = FALSE)

  # Function to determine if pbapply is installed. If it is installed, it will
  # display a progress bar
  list_apply <- function(x, fun, ...) {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(x, fun, ...)
    } else {
      lapply(x, fun, ...)
    }
  }

  # Run each row of param_grid
  coi_pred <- list_apply(seq_len(nrow(param_grid)), function(x) {

    # Run each sample repetitions times
    repeats <- lapply(seq_len(repetitions), function(y) {
      test_result <-
        single_sensitivity(param_grid$coi[x],
                           param_grid$max_coi[x],
                           plaf,
                           param_grid$coverage[x],
                           param_grid$alpha[x],
                           param_grid$overdispersion[x],
                           param_grid$relatedness[x],
                           param_grid$epsilon[x],
                           param_grid$seq_error[x],
                           param_grid$bin_size[x],
                           param_grid$comparison[x],
                           param_grid$distance[x],
                           param_grid$weighted[x],
                           param_grid$coi_method[x])
      return (test_result)
    })

    return(repeats)
  })

  # Extract the COIs from the result list
  extracted_cois <- lapply(coi_pred, function(x) {
    unlist(lapply(x, function(i) { i$coi }))
  })

  # Extract the probabilities from the result list
  extracted_probs <- lapply(coi_pred, function(x) {
    # Get the probabilities
    matrix <- do.call(rbind, lapply(x, function(i) { i$probability }))

    # Name the matrix and return
    colnames(matrix) <- paste("coi", 1:max_coi, sep = "_")
    rownames(matrix) <- paste("rep", seq(repetitions), sep = "_")

    # Add a average row to the matrix
    matrix <- rbind(colMeans(matrix), matrix)
    rownames(matrix)[1] <- "average"

    return(matrix)
  })

  ## Naming
  # Determine how many unique COIs there are
  num_cois = length(coi)

  # Calculate the number of times each COI is repeated in param_grid
  num_repeat_cois = length(param_grid$coi) / num_cois

  # Name the predicted COIs
  names(extracted_cois) <- paste("coi",
                                 param_grid$coi,
                                 rep(seq(num_repeat_cois), each = num_cois),
                                 sep="_")
  names(extracted_probs) <- names(extracted_cois)

  ## Calculations
  # Calculate mean absolute error
  len <- nrow(param_grid)
  boot_mae <- lapply(seq_len(len), function(x) {
    # Get the specific row (list of predicted COIs). Note that we convert to a
    # data frame because the boot package requires this.
    boot_data <- as.data.frame(extracted_cois[x])

    # Create function to find mean absolute error for bootstrapping
    mae <- function(data, true_coi, indices) {
      # Sample from the data
      sampled <- data[indices, ]

      # Compute mean absolute error for the sampled data and return
      sampled_mae <- sum(abs(sampled - true_coi)) / length(sampled)
      return(sampled_mae)
    }

    # Bootstrapping
    results <- boot::boot(data = boot_data, true_coi = param_grid$coi[x],
                          statistic = mae, R = 1000)

    # Get the normal confidence interval
    invisible(utils::capture.output(CI <- boot::boot.ci(results,
                                                        type = "norm")$norm))

    # Store the mean absolute error and confidence interval bounds
    extract <- list(mae = results$t0, lower = CI[2], upper = CI[3])
    return(extract)
  })

  # Convert output of bootstrapping to data frame structure
  boot_mae <- as.data.frame(do.call(rbind, boot_mae))

  # Calculate bias (mean error)
  coi_bias <- lapply(seq_len(len), function(x) {
    sum(extracted_cois[[x]] - param_grid$coi[x]) / length(extracted_cois[[x]])
  })

  # Save changing parameters with bootstrapping
  boot_error <- param_grid %>%
    dplyr::select_if(function(x) dplyr::n_distinct(x) > 1)
  boot_error <- dplyr::bind_cols(boot_error, boot_mae)
  boot_error$bias  <- unlist(coi_bias)

  # Warnings for when could not compute CI. We want to show at most 5 cases
  # where there is an issue. Anymore, and the warning message becomes too
  # confusing.
  warn_tibble <- boot_error %>%
    tidyr::unchop(cols = tidyr::everything()) %>%
    dplyr::filter(is.na(.data$lower) | is.na(.data$upper))
  if (nrow(warn_tibble) >= 1) {
    bullet <- ""
    if (nrow(warn_tibble) < 5) {
      for (i in seq(nrow(warn_tibble))) {
        bullet <- glue::glue("{bullet}\n\u2716 COI of {warn_tibble$coi[i]} failed.")
      }
    } else if (nrow(warn_tibble) >= 5) {
      for (i in seq(5)) {
        bullet <- glue::glue("{bullet}\n\u2716 COI of {warn_tibble$coi[i]} failed.")
      }
      bullet <- glue::glue("{bullet}\n... and {nrow(warn_tibble) - 5} more failed")
    }

    message <- glue::glue("Can't calculate bootstrapped confidence interval:",
                          "{bullet}")
    warning(message, call. = FALSE)
  }

  # Return predicted COIs and param_grid
  ret <- list(predicted_coi = as.data.frame(extracted_cois),
              probability   = extracted_probs,
              param_grid    = param_grid,
              boot_error    = boot_error)
}

#------------------------------------------------
#' Continuous sensitivity analysis
#'
#' Runs several iterations of a full COI sensitivity analysis with varying
#' parameters.
#'
#' @param repetitions The number of times each sample will be run.
#' @inheritParams single_sensitivity
#'
#' @return A list of the following:
#' * `predicted_coi`: A dataframe of the predicted COIs. COIs are
#'   predicted using [compute_coi()]. Each column represents a separate set
#'   of parameters. Each row represents a predicted COI. Predictions are done
#'   many times, depending on the value of `repetitions`.
#' * `probability`:A list of matrices containing the probability
#'   that our model predicted each COI value. Each row contains the probability
#'   for a different run. The first row contains the average probabilities over
#'   all the runs.
#' * `param_grid`: The parameter grid. The parameter grid is all
#'   possible combinations of the parameters inputted. Each row represents a
#'   unique combination.
#' * `boot_error`: A dataframe containing information about the error
#'   of the algorithm. The first column indicates the COI that was fed into the
#'   simulation. The other columns indicate the mean absolute error (mae),
#'   the lower and upper bounds of the 95% confidence interval and the bias.
#'
#' @export
cont_sensitivity <- function(repetitions = 10,
                             coi = 3,
                             max_coi = 25,
                             plaf = runif(1000, 0, 0.5),
                             coverage = 200,
                             alpha = 1,
                             overdispersion = 0,
                             relatedness = 0,
                             epsilon = 0,
                             seq_error = NULL,
                             bin_size = 20,
                             comparison = "overall",
                             distance = "squared",
                             weighted = TRUE,
                             coi_method = "1") {

  # Check inputs
  assert_pos_int(repetitions)
  assert_pos_int(coi)
  assert_single_pos_int(max_coi)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_pos_int(coverage)
  assert_pos(alpha, zero_allowed = FALSE)
  assert_pos(overdispersion)
  assert_bounded(relatedness)
  assert_bounded(epsilon)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)
  assert_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_logical(weighted)
  assert_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Create parameter grid
  param_grid <- expand.grid(coi = coi,
                            max_coi = max_coi,
                            coverage = coverage,
                            alpha = alpha,
                            overdispersion = overdispersion,
                            relatedness = relatedness,
                            epsilon = epsilon,
                            seq_error = seq_error,
                            bin_size = bin_size,
                            comparison = comparison,
                            distance = distance,
                            weighted = weighted,
                            coi_method = coi_method,
                            stringsAsFactors = FALSE)

  # Function to determine if pbapply is installed. If it is installed, it will
  # display a progress bar
  list_apply <- function(x, fun, ...) {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(x, fun, ...)
    } else {
      lapply(x, fun, ...)
    }
  }

  # Run each row of param_grid
  coi_pred <- list_apply(seq_len(nrow(param_grid)), function(x) {

    # Run each sample repetitions times
    repeats <- lapply(seq_len(repetitions), function(y) {
      test_sim <- sim_biallelic(param_grid$coi[x],
                                plaf,
                                param_grid$coverage[x],
                                param_grid$alpha[x],
                                param_grid$overdispersion[x],
                                param_grid$relatedness[x],
                                param_grid$epsilon[x])

      test_result <- optimize_coi(test_sim,
                                  "sim",
                                  param_grid$max_coi[x],
                                  param_grid$seq_error[x],
                                  param_grid$bin_size[x],
                                  param_grid$distance[x],
                                  param_grid$weighted[x],
                                  param_grid$coi_method[x])
      return (test_result)
    })

    return(repeats)
  })

  # Extract the COIs from the result list
  extracted_cois <- lapply(coi_pred, function(x) {
    unlist(x)
  })

  ## Naming
  # Determine how many unique COIs there are
  num_cois = length(coi)

  # Calculate the number of times each COI is repeated in param_grid
  num_repeat_cois = length(param_grid$coi) / num_cois

  # Name the predicted COIs
  names(extracted_cois) <- paste("coi",
                                 param_grid$coi,
                                 rep(seq(num_repeat_cois), each = num_cois),
                                 sep="_")

  ## Calculations
  # Calculate mean absolute error
  len <- nrow(param_grid)
  boot_mae <- lapply(seq_len(len), function(x) {
    # Get the specific row (list of predicted COIs). Note that we convert to a
    # data frame because the boot package requires this.
    boot_data <- as.data.frame(extracted_cois[x])

    # Create function to find mean absolute error for bootstrapping
    mae <- function(data, true_coi, indices) {
      # Sample from the data
      sampled <- data[indices, ]

      # Compute mean absolute error for the sampled data and return
      sampled_mae <- sum(abs(sampled - true_coi)) / length(sampled)
      return(sampled_mae)
    }

    # Bootstrapping
    results <- boot::boot(data = boot_data, true_coi = param_grid$coi[x],
                          statistic = mae, R = 1000)

    # Get the normal confidence interval
    invisible(utils::capture.output(CI <- boot::boot.ci(results,
                                                        type = "norm")$norm))

    # Store the mean absolute error and confidence interval bounds
    extract <- list(mae = results$t0, lower = CI[2], upper = CI[3])
    return(extract)
  })

  # Convert output of bootstrapping to data frame structure
  boot_mae <- as.data.frame(do.call(rbind, boot_mae))

  # Calculate bias (mean error)
  coi_bias <- lapply(seq_len(len), function(x) {
    sum(extracted_cois[[x]] - param_grid$coi[x]) / length(extracted_cois[[x]])
  })

  # Save changing parameters with bootstrapping
  boot_error <- param_grid %>%
    dplyr::select_if(function(x) dplyr::n_distinct(x) > 1)
  boot_error <- dplyr::bind_cols(boot_error, boot_mae)
  boot_error$bias  <- unlist(coi_bias)

  # Warnings for when could not compute CI. We want to show at most 5 cases
  # where there is an issue. Anymore, and the warning message becomes too
  # confusing.
  warn_tibble <- boot_error %>%
    tidyr::unchop(cols = tidyr::everything()) %>%
    dplyr::filter(is.na(.data$lower) | is.na(.data$upper))
  if (nrow(warn_tibble) >= 1) {
    bullet <- ""
    if (nrow(warn_tibble) < 5) {
      for (i in seq(nrow(warn_tibble))) {
        bullet <- glue::glue("{bullet}\n\u2716 COI of {warn_tibble$coi[i]} failed.")
      }
    } else if (nrow(warn_tibble) >= 5) {
      for (i in seq(5)) {
        bullet <- glue::glue("{bullet}\n\u2716 COI of {warn_tibble$coi[i]} failed.")
      }
      bullet <- glue::glue("{bullet}\n... and {nrow(warn_tibble) - 5} more failed")
    }

    message <- glue::glue("Can't calculate bootstrapped confidence interval:",
                          "{bullet}")
    warning(message, call. = FALSE)
  }

  # Return predicted COIs and param_grid
  ret <- list(predicted_coi = as.data.frame(extracted_cois),
              param_grid    = param_grid,
              boot_error    = boot_error)
}
