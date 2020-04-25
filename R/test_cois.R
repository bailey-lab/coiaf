#------------------------------------------------
#' @title Run a single COI test
#'
#' @description Runs a single full COI test.
#'
#' @param COI_range The range of COIs to test
#' @inheritParams sim_biallelic
#' @inheritParams simulated_coi
#' @param cut How often the data is summarized
#' @inheritParams compute_coi
#'
#' @return Computed COI values
#'
#' @keywords internal

run_coi_test <- function(COI = 3,
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

  return (calc_coi)
}


#------------------------------------------------
#' @title Test COIs
#'
#' @description Runs several iterations of a full COI test with varying
#' parameters
#'
#' @param repetitions The number of times each sample will be run
#' @param COI_range The range of COIs to test
#' @inheritParams sim_biallelic
#' @inheritParams simulated_coi
#' @param cut How often the data is summarized
#' @inheritParams compute_coi
#'
#' @return Computed COI values
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

  # Run each row of param_grid
  coi_pred <- lapply(seq_len(nrow(param_grid)), function(x) {

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
                     weighted)
      return (test_result)
    }, FUN.VALUE = numeric(1))

    return(repeats)
  })

  # Name the predicted COIs
  names(coi_pred) <- paste("coi_", param_grid$COI, sep="")

  # Return predicted COIs and param_grid
  ret <- list(predicted_coi = as.data.frame(coi_pred),
              param_grid    = param_grid)
  return (ret)
}
