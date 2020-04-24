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



