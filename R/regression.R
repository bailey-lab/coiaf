#' Compute COI based on residuals of all loci against theoretical curves
#'
#' @inheritParams compute_coi
#' @param seq_error_bin_size Number of loci in smallest bin for estimating
#'   sequence error
#' @return A list of the following:
#' * `coi`: The predicted COI of the sample.
#' * `probability`: A probability density function representing the probability
#'  of each COI.
#'
#' @export
compute_coi_regression <- function(data,
                                   data_type,
                                   max_coi = 25,
                                   seq_error = 0.01,
                                   distance = "squared",
                                   coi_method = "variant",
                                   seq_error_bin_size = 20) {
  # Warnings for deprecated arguments
  if (distance != "squared") {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "compute_coi_regression(distance)",
      details = 'The distance method will be fixed to "squared" in the next release.'
    )
  }

  processed_data <- process_data_for_regression(
    data, data_type, max_coi, seq_error,
    distance, coi_method, seq_error_bin_size
  )

  # Special case for the Frequency Method where there is no data
  if (coi_method == "frequency" & nrow(processed_data) == 0) {
    return(list(
      coi = NaN,
      probability = c(1, rep(0, max_coi - 1)),
      notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method.",
      estimated_coi = 1,
      num_variant_loci = 0
    ))
  }

  # Calculate theoretical COI curves for the interval specified. Since we want
  # the theoretical curves and the simulated curves to have the PLMAF values, we
  # compute the theoretical COI curves at processed_data$midpoints
  theory_cois <- theoretical_coi(
    1:max_coi,
    processed_data$plmaf,
    coi_method
  )

  # create our distance metric
  overall_res <- distance_curves(processed_data, theory_cois, distance)

  # Extract information from the helper function
  coi <- overall_res$coi
  dist <- overall_res$dist

  # Distance to probability
  dist <- as.numeric(dist)
  dist <- 1 / (dist + 1e-5)
  dist <- dist / sum(dist, na.rm = T)
  dist[is.nan(dist)] <- 0

  # Special case for the Frequency Method
  if (coi_method == "frequency") {
    check <- switch(data_type,
      "sim" = check_freq_method(data$data$wsmaf, data$data$plmaf, seq_error),
      "real" = check_freq_method(data$wsmaf, data$plmaf, seq_error)
    )

    # If the actual number of variant sites is less than the lower bound of the
    # CI, the COI may be 1
    if (check[["variant"]] < check[["lower_ci"]]) {
      ret <- list(
        coi = NaN,
        probability = c(1, rep(0, max_coi - 1)),
        notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method.",
        estimated_coi = as.numeric(coi),
        num_variant_loci = check[["variant"]],
        expected_num_loci = check[["expected"]]
      )
      return(ret)
    }
  }

  # List to return
  return(list(coi = as.numeric(coi), probability = dist))
}

#' Compute COI based on all points fitted to best fitting curve for COI
#'
#' @inheritParams optimize_coi
#' @param seq_error_bin_size Number of loci in smallest bin for estimating
#'   sequence error
#' @return The predicted COI value.
#' @seealso [stats::optim()] for the complete documentation on the optimization
#' function.
#' @family optimization functions
#' @export
optimize_coi_regression <- function(data,
                                    data_type,
                                    max_coi = 25,
                                    seq_error = 0.01,
                                    distance = "squared",
                                    coi_method = "variant",
                                    seq_error_bin_size = 20) {
  # Warnings for deprecated arguments
  if (distance != "squared") {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "optimize_coi_regression(distance)",
      details = 'The distance method will be fixed to "squared" in the next release.'
    )
  }

  processed_data <- process_data_for_regression(
    data, data_type, max_coi, seq_error,
    distance, coi_method, seq_error_bin_size
  )

  # Special case for the Frequency Method where there is no data
  if (coi_method == "frequency" & nrow(processed_data) == 0) {
    return(structure(
      NaN,
      notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method.",
      estimated_coi = 1,
      num_variant_loci = 0
    ))
  }

  # correct names for likelihood function
  names(processed_data) <- c("midpoints", "m_variant", "coverage", "bucket_size")

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
    cli_warn(c(
      "The model did not converge:",
      "x" = "{bullet}",
      "x" = "Output of optim: {fit$message}."
    ))
  }

  # Estimated COI
  estimated_coi <- round(fit$par, 4)

  # Special case for the Frequency Method
  if (coi_method == "frequency") {
    check <- switch(data_type,
      "sim" = check_freq_method(data$data$wsmaf, data$data$plmaf, seq_error),
      "real" = check_freq_method(data$wsmaf, data$plmaf, seq_error)
    )

    # If the actual number of variant sites is less than the lower bound of the
    # CI, the COI may be 1
    if (check[["variant"]] < check[["lower_ci"]]) {
      return(structure(
        NaN,
        notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method.",
        estimated_coi = estimated_coi,
        num_variant_loci = check[["variant"]],
        expected_num_loci = check[["expected"]]
      ))
    }
  }

  estimated_coi
}

#' @noRd
process_data_for_regression <- function(data,
                                        data_type,
                                        max_coi,
                                        seq_error,
                                        distance,
                                        coi_method,
                                        seq_error_bin_size) {
  # Process data to get the wsmaf and plmaf
  if (data_type == "sim") {
    wsmaf <- data$data$wsmaf
    plmaf <- data$data$plmaf
    coverage <- data$data$coverage
  } else if (data_type == "real") {
    data <- check_input_data(data, "real")
    wsmaf <- data$wsmaf
    plmaf <- data$plmaf
    coverage <- data$coverage
  }

  # Infer value of seq_error if NULL
  if (is.null(seq_error)) {
    seq_error <- estimate_seq_error(wsmaf, plmaf, seq_error_bin_size)
  }

  # Process data to get in right format
  if (coi_method == "variant") {
    # Isolate PLMAF and whether a site is a variant, accounting for sequence
    # error
    df <- data.frame(
      plmaf = plmaf,
      m_variant = ifelse(wsmaf <= seq_error | wsmaf >= (1 - seq_error), 0, 1),
      coverage = coverage,
      bucket_size = 1
    )
  } else if (coi_method == "frequency") {
    # Subset to heterozygous sites
    df <- data.frame(
      plmaf = plmaf,
      m_variant = wsmaf,
      coverage = coverage,
      bucket_size = 1
    ) %>%
      dplyr::filter(
        .data$m_variant > seq_error & .data$m_variant < (1 - seq_error)
      )
  }

  df
}
