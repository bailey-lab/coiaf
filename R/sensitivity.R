#------------------------------------------------
#' Single sensitivity analysis
#'
#' Runs a single full COI sensitivity analysis.
#'
#' @param disc_or_cont Whether to run a discrete or continuous COI estimation.
#' @inheritParams sim_biallelic
#' @inheritParams compute_coi
#'
#' @return Predicted COI value.
#'
#' @keywords internal
single_sensitivity <- function(disc_or_cont,
                               coi = 3,
                               max_coi = 25,
                               plmaf = runif(1000, 0, 0.5),
                               coverage = 200,
                               alpha = 1,
                               overdispersion = 0,
                               relatedness = 0,
                               epsilon = 0,
                               seq_error = 0.01,
                               bin_size = 20,
                               comparison = "overall",
                               distance = "squared",
                               coi_method = "variant",
                               use_bins = FALSE) {
  # Check inputs
  assert_single_pos_int(coi)
  assert_single_pos_int(max_coi)
  assert_vector(plmaf)
  assert_bounded(plmaf, left = 0, right = 0.5)
  assert_single_pos_int(coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion)
  assert_single_bounded(relatedness)
  assert_single_bounded(epsilon)
  if (!is.null(seq_error) & !is.na(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)
  assert_single_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_single_string(coi_method)
  assert_in(coi_method, c("variant", "frequency"))

  # Simulate data
  sim_data <- sim_biallelic(
    coi,
    plmaf,
    coverage,
    alpha,
    overdispersion,
    relatedness,
    epsilon
  )

  # Workaround for when seq_error = NULL. Have seq_error saved as NA so
  # our general sensitivity function can deal with it. But need to convert back
  # to NULL.
  if (is.na(seq_error)) seq_error <- NULL

  # Compute COI
  if (disc_or_cont == "disc") {
    compute_coi(
      sim_data,
      "sim",
      max_coi,
      seq_error,
      bin_size,
      comparison,
      distance,
      coi_method,
      use_bins
    )
  } else if (disc_or_cont == "cont") {
    optimize_coi(
      sim_data,
      "sim",
      max_coi,
      seq_error,
      bin_size,
      distance,
      coi_method,
      use_bins
    )
  }
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
                        plmaf = runif(1000, 0, 0.5),
                        coverage = 200,
                        alpha = 1,
                        overdispersion = 0,
                        relatedness = 0,
                        epsilon = 0,
                        seq_error = 0.01,
                        bin_size = 20,
                        comparison = "overall",
                        distance = "squared",
                        coi_method = "variant",
                        use_bins = FALSE) {
  # Check inputs
  assert_pos_int(repetitions)
  assert_pos_int(coi)
  assert_single_pos_int(max_coi)
  assert_vector(plmaf)
  assert_bounded(plmaf, left = 0, right = 0.5)
  assert_pos_int(coverage)
  assert_pos(alpha, zero_allowed = FALSE)
  assert_pos(overdispersion)
  assert_bounded(relatedness)
  assert_bounded(epsilon)
  if (!is.null(seq_error)) assert_bounded(seq_error)
  assert_pos_int(bin_size)
  assert_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_string(coi_method)
  assert_in(coi_method, c("variant", "frequency"))

  # Create parameter grid
  param_grid <- tidyr::expand_grid(
    # Simulation parameters
    coi = coi,
    coverage = coverage,
    alpha = alpha,
    overdispersion = overdispersion,
    relatedness = relatedness,
    epsilon = epsilon,

    # Estimation parameters
    max_coi = max_coi,
    seq_error = ifelse(is.null(seq_error), NA, seq_error),
    bin_size = bin_size,
    comparison = comparison,
    distance = distance,
    coi_method = coi_method,
    use_bins = use_bins,
  )

  # Run each row of param_grid
  coi_pred <- lapply(
    cli::cli_progress_along(seq_len(nrow(param_grid)), "Estimating the COI"),
    function(x) {
      # Run each sample repetitions times
      lapply(seq_len(repetitions), function(y) {
        single_sensitivity(
          disc_or_cont = "disc",
          coi = param_grid$coi[x],
          max_coi = param_grid$max_coi[x],
          plmaf = plmaf,
          coverage = param_grid$coverage[x],
          alpha = param_grid$alpha[x],
          overdispersion = param_grid$overdispersion[x],
          relatedness = param_grid$relatedness[x],
          epsilon = param_grid$epsilon[x],
          seq_error = param_grid$seq_error[x],
          bin_size = param_grid$bin_size[x],
          comparison = param_grid$comparison[x],
          distance = param_grid$distance[x],
          coi_method = param_grid$coi_method[x],
          use_bins = param_grid$use_bins[x]
        )
      })
    }
  )

  # Add names to lists
  names(coi_pred) <- paste0("param_set_", seq_len(nrow(param_grid)))
  coi_pred <- lapply(coi_pred, function(x) {
    names(x) <- paste0("rep_", seq_len(length(x)))
    x
  })

  # Extract the COIs from the result list
  extracted_cois <- lapply(coi_pred, function(x) {
    unlist(lapply(x, function(i) {
      i$coi
    }))
  })

  # Extract the probabilities from the result list
  extracted_probs <- lapply(coi_pred, function(x) {
    # Get the probabilities
    matrix <- do.call(rbind, lapply(x, function(i) {
      i$probability
    }))

    # Name the matrix and return
    colnames(matrix) <- paste("coi", 1:max_coi, sep = "_")
    rownames(matrix) <- paste("rep", seq(repetitions), sep = "_")

    # Add a average row to the matrix
    matrix <- rbind(colMeans(matrix), matrix)
    rownames(matrix)[1] <- "average"

    matrix
  })

  ## Naming
  # Determine how many unique COIs there are
  num_cois <- length(coi)

  # Calculate the number of times each COI is repeated in param_grid
  num_repeat_cois <- length(param_grid$coi) / num_cois

  # Name the predicted COIs
  names(extracted_cois) <- paste(
    "coi",
    param_grid$coi,
    rep(seq(num_repeat_cois), each = num_cois),
    sep = "_"
  )
  names(extracted_probs) <- names(extracted_cois)

  ## Calculations
  # Calculate mean absolute error
  boot_mae <- boot_mae(param_grid, extracted_cois)

  # Calculate bias (mean error)
  coi_bias <- bias(param_grid, extracted_cois)

  # Save changing parameters with bootstrapping
  boot_error <- dplyr::select(param_grid, where(~ dplyr::n_distinct(.x) > 1))
  boot_error <- dplyr::bind_cols(boot_error, boot_mae) %>%
    tidyr::unnest(cols = c(mae, lower, upper))
  boot_error$bias <- unlist(coi_bias)

  # Customize warnings for when could not compute CI.
  sa_warn(boot_error)

  # Return predicted COIs and param_grid
  list(
    predicted_coi = as.data.frame(extracted_cois),
    probability   = extracted_probs,
    param_grid    = param_grid,
    boot_error    = boot_error
  )
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
                             plmaf = runif(1000, 0, 0.5),
                             coverage = 200,
                             alpha = 1,
                             overdispersion = 0,
                             relatedness = 0,
                             epsilon = 0,
                             seq_error = 0.01,
                             bin_size = 20,
                             comparison = "overall",
                             distance = "squared",
                             coi_method = "variant",
                             use_bins = FALSE) {
  # Check inputs
  assert_pos_int(repetitions)
  assert_pos_int(coi)
  assert_single_pos_int(max_coi)
  assert_vector(plmaf)
  assert_bounded(plmaf, left = 0, right = 0.5)
  assert_pos_int(coverage)
  assert_pos(alpha, zero_allowed = FALSE)
  assert_pos(overdispersion)
  assert_bounded(relatedness)
  assert_bounded(epsilon)
  if (!is.null(seq_error)) assert_bounded(seq_error)
  assert_pos_int(bin_size)
  assert_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_string(coi_method)
  assert_in(coi_method, c("variant", "frequency"))

  # Create parameter grid
  param_grid <- tidyr::expand_grid(
    coi = coi,
    max_coi = max_coi,
    coverage = coverage,
    alpha = alpha,
    overdispersion = overdispersion,
    relatedness = relatedness,
    epsilon = epsilon,
    seq_error = ifelse(is.null(seq_error), NA, seq_error),
    bin_size = bin_size,
    comparison = comparison,
    distance = distance,
    coi_method = coi_method,
    use_bins = use_bins,
  )

  # Run each row of param_grid
  coi_pred <- lapply(
    cli::cli_progress_along(seq_len(nrow(param_grid)), "Estimating the COI"),
    function(x) {
      # Run each sample repetitions times
      lapply(seq_len(repetitions), function(y) {
        single_sensitivity(
          disc_or_cont = "cont",
          coi = param_grid$coi[x],
          max_coi = param_grid$max_coi[x],
          plmaf = plmaf,
          coverage = param_grid$coverage[x],
          alpha = param_grid$alpha[x],
          overdispersion = param_grid$overdispersion[x],
          relatedness = param_grid$relatedness[x],
          epsilon = param_grid$epsilon[x],
          seq_error = param_grid$seq_error[x],
          bin_size = param_grid$bin_size[x],
          comparison = param_grid$comparison[x],
          distance = param_grid$distance[x],
          coi_method = param_grid$coi_method[x],
          use_bins = param_grid$use_bins[x]
        )
      })
    }
  )

  # Add names to lists
  names(coi_pred) <- paste0("param_set_", seq_len(nrow(param_grid)))
  coi_pred <- lapply(coi_pred, function(x) {
    names(x) <- paste0("rep_", seq_len(length(x)))
    x
  })

  # Extract the COIs from the result list
  extracted_cois <- lapply(coi_pred, function(x) unlist(x))

  ## Naming
  # Determine how many unique COIs there are
  num_cois <- length(coi)

  # Calculate the number of times each COI is repeated in param_grid
  num_repeat_cois <- length(param_grid$coi) / num_cois

  # Name the predicted COIs
  names(extracted_cois) <- paste(
    "coi",
    param_grid$coi,
    rep(seq(num_repeat_cois), each = num_cois),
    sep = "_"
  )

  ## Calculations
  # Calculate mean absolute error
  boot_mae <- boot_mae(param_grid, extracted_cois)

  # Calculate bias (mean error)
  coi_bias <- bias(param_grid, extracted_cois)

  # Save changing parameters with bootstrapping
  boot_error <- dplyr::select(param_grid, where(~ dplyr::n_distinct(.x) > 1))
  boot_error <- dplyr::bind_cols(boot_error, boot_mae)
  boot_error$bias <- unlist(coi_bias)

  # Customize warnings for when could not compute CI.
  sa_warn(boot_error)

  # Return predicted COIs and param_grid
  list(
    predicted_coi = as.data.frame(extracted_cois),
    param_grid    = param_grid,
    boot_error    = boot_error
  )
}

# Calculate boostrapped mean absolute error
boot_mae <- function(param_grid, extracted_cois) {
  mae_list <- lapply(
    cli::cli_progress_along(seq_len(nrow(param_grid)), "Computing statistics"),
    function(x) {
      # Get the specific row (list of predicted COIs). Note that we convert to a
      # data frame because the boot package requires this.
      boot_data <- as.data.frame(extracted_cois[x])

      # If data has NaN in it, return all NaN
      if (any(is.nan(unlist(boot_data)))) {
        return(list(mae = NaN, lower = NaN, upper = NaN))
      }

      # Create function to find mean absolute error for bootstrapping
      mae <- function(data, true_coi, indices) {
        # Sample from the data
        sampled <- data[indices, ]

        # Compute mean absolute error for the sampled data and return
        sampled_mae <- sum(abs(sampled - true_coi)) / length(sampled)
        return(sampled_mae)
      }

      # Bootstrapping
      results <- boot::boot(
        data = boot_data,
        true_coi = param_grid$coi[x],
        statistic = mae,
        R = 1000
      )

      # Get the normal confidence interval
      invisible(utils::capture.output(
        CI <- boot::boot.ci(results, type = "norm")$normal
      ))

      # Store the mean absolute error and confidence interval bounds
      list(mae = results$t0, lower = CI[2], upper = CI[3])
    }
  )

  # Convert output of bootstrapping to data frame structure
  as.data.frame(do.call(rbind, mae_list))
}

# Calculate bias (mean error)
bias <- function(param_grid, extracted_cois) {
  lapply(seq_len(nrow(param_grid)), function(x) {
    sum(extracted_cois[[x]] - param_grid$coi[x]) / length(extracted_cois[[x]])
  })
}

# Customize warnings for when could not compute CI.
sa_warn <- function(boot_error) {
  warn_tibble <- boot_error %>%
    tidyr::unchop(cols = tidyr::everything()) %>%
    dplyr::filter(is.na(.data$lower) | is.na(.data$upper))
  warn_nan_tibble <- dplyr::filter(warn_tibble, is.nan(.data$mae))
  cli::cli_warn(c(
    "Unable to calculate bootstrapped confidence interval.",
    "x" = "{nrow(warn_tibble) - nrow(warn_nan_tibble)} computation{?s} failed.",
    "i" = "{nrow(warn_nan_tibble)} computation{?s} {?was/were} not applicable."
  ))
}
