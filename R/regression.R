#' Compute COI based on residauals of all loci against theoretical curves
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
                                   seq_error = NULL,
                                   distance = "squared",
                                   coi_method = "variant",
                                   seq_error_bin_size = 20) {

  processed_data <- process_data_for_regression(
    data, data_type, max_coi, seq_error,
    distance, coi_method, seq_error_bin_size
  )

  # was this deemed to be COI = 1
  if ("coi" %in% names(processed_data)) {
    return(processed_data)
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

  # List to return
  return(list(coi = as.numeric(coi), probability = dist))

}

#' Compute COI based on all points fitted to best fitting curve for COI
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
                                    seq_error = NULL,
                                    distance = "squared",
                                    coi_method = "variant",
                                    seq_error_bin_size = 20) {


  processed_data <- process_data_for_regression(
    data, data_type, max_coi, seq_error,
    distance, coi_method, seq_error_bin_size
  )
  names(processed_data) <- c("midpoints", "m_variant", "coverage", "bucket_size")

  # was this deemed to be COI = 1
  if ("coi" %in% names(processed_data)) {
    return(processed_data)
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
  return(round(fit$par, 4))

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

    minor <- check_real_data(data$wsmaf, data$plmaf)
    wsmaf <- minor$wsmaf
    plmaf <- minor$plmaf
    coverage <- data$coverage

  }

  # Infer value of seq_error if NULL
  if (is.null(seq_error)) {
    seq_error <- estimate_seq_error(wsmaf, plmaf, seq_error_bin_size)
  }

  # Process data to get in right format
  if (coi_method == "variant") {

    # Isolate PLMAF and whether a site is a variant,
    # accounting for sequence error
    df <- data.frame(
      plmaf = plmaf,
      m_variant = ifelse(wsmaf <= seq_error | wsmaf >= (1 - seq_error), 0, 1),
      coverage = coverage
    )

  } else if (coi_method == "frequency") {

    # do a variant number check here
    check <- check_freq_method(wsmaf, plmaf, seq_error)

    # if checks out then return 1 here
    if (!check) {
      ret <- list(
        coi = 1,
        probability = c(1, rep(0, max_coi - 1)),
        notes = "Too few variant loci suggesting that the COI is 1 based on the Variant Method."
      )
      return(ret)
    }

    # Subset to heterozygous sites
    df <- data.frame(plmaf = plmaf, m_variant = wsmaf, coverage = coverage) %>%
      dplyr::filter(wsmaf > seq_error & wsmaf < (1 - seq_error))

  }

  df$bucket_size <- 1

  return(df)

}
