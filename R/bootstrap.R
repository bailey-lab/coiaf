#------------------------------------------------
#' Generate bootstrapped CI
#'
#' Generate bootstrapped confidence interval for COI estimates.
#'
#' @inheritParams compute_coi
#' @param solution_method Whether to estimate discrete or continuous COIs.
#' @param replicates The number of bootstrap replicates.
#' @param parallel Whether to parallelize the confidence interval calculation.
#' @param ncpus The number of processes to be used in parallel operation.
#'
#' @export
bootstrap_ci <- function(data,
                         max_coi = 25,
                         seq_error = 0.01,
                         coi_method = c("variant", "frequency"),
                         solution_method = c("discrete", "continuous"),
                         use_bins = FALSE,
                         bin_size = 20,
                         replicates = 100,
                         parallel = TRUE,
                         ncpus = 8) {
  UseMethod("bootstrap_ci")
}

#' @rdname bootstrap_ci
#' @export
bootstrap_ci.default <- function(data,
                                 max_coi = 25,
                                 seq_error = 0.01,
                                 coi_method = c("variant", "frequency"),
                                 solution_method = c("discrete", "continuous"),
                                 use_bins = FALSE,
                                 bin_size = 20,
                                 replicates = 100,
                                 parallel = TRUE,
                                 ncpus = 8) {
  # Argument matching
  coi_method <- rlang::arg_match(coi_method)
  solution_method <- rlang::arg_match(solution_method)

  # Define statistic function
  bootstrap_statistic <- function(data, indices, solution_method, ...) {
    data <- data[indices, ]

    # Determine function call based on solution_method
    switch(solution_method,
      discrete = compute_coi(data, "real", ...)$coi,
      continuous = optimize_coi(data, "real", ...)
    )
  }

  # Bootstrap
  boot_out <- boot::boot(
    data = data,
    statistic = bootstrap_statistic,
    R = replicates,
    parallel = if (parallel) "multicore" else "no",
    ncpus = ncpus,
    solution_method = solution_method,
    max_coi = 25,
    seq_error = 0.01,
    coi_method = coi_method,
    use_bins = FALSE,
    bin_size = 20
  )

  tryCatch(
    broom::tidy(boot_out, conf.int = TRUE) %>%
      dplyr::rename(coi = .data$statistic),
    error = function(e) {
      broom::tidy(boot_out, conf.int = FALSE) %>%
        tibble::add_column(conf.low = NaN, conf.high = NaN) %>%
        dplyr::rename(coi = .data$statistic)
    }
  )
}
