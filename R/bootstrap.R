#------------------------------------------------
#' Generate bootstrapped CI
#'
#' Generate bootstrapped confidence interval for COI estimates.
#'
#' @inheritParams compute_coi
#' @param solution_method Whether to estimate discrete or continuous COIs.
#' @param replicates The number of bootstrap replicates.
#' @param parallel Whether to parallelize the confidence interval calculation.
#'   Note that parallelization only works on non-Windows machines.
#' @param ncpus The number of processes to be used in parallel operation.
#'
#' @return
#' A [`tibble()`][tibble::tibble-package] with columns:
#' \describe{
#'   \item{coi}{The mean COI.}
#'   \item{bias}{Bias of the statistic.}
#'   \item{std.error}{The standard error of the statistic.}
#'   \item{conf.low}{The lower 95% confidence interval.}
#'   \item{conf.high}{The upper 95% confidence interval.}
#' }
#'
#' @seealso [boot::boot()], [boot::boot.ci()], [broom::tidy.boot()]
#' @export
#' @examples
#' sim_data <- sim_biallelic(coi = 5, plmaf = runif(100, 0, 0.5))
#' bootstrap_ci(sim_data, solution_method = "continuous")
bootstrap_ci <- function(data,
                         max_coi = 25,
                         seq_error = 0.01,
                         coi_method = c("variant", "frequency"),
                         solution_method = c("discrete", "continuous"),
                         use_bins = FALSE,
                         bin_size = 20,
                         replicates = 100,
                         parallel = FALSE,
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
                                 parallel = FALSE,
                                 ncpus = 8) {
  # Argument matching
  coi_method <- rlang::arg_match(coi_method)
  solution_method <- rlang::arg_match(solution_method)

  # Warning for deprecated arguments
  if (use_bins) {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "bootstrap_ci(use_bins)",
      details = "The ability to use bins to estimate the COI will be dropped in the next release."
    )
  }
  if (bin_size != 20) {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "bootstrap_ci(bin_size)",
      details = "The ability to use bins to estimate the COI will be dropped in the next release."
    )
  }

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
      dplyr::rename(coi = "statistic"),
    error = function(e) {
      broom::tidy(boot_out, conf.int = FALSE) %>%
        tibble::add_column(conf.low = NaN, conf.high = NaN) %>%
        dplyr::rename(coi = "statistic")
    }
  )
}

#' @rdname bootstrap_ci
#' @export
bootstrap_ci.sim <- function(data,
                             max_coi = 25,
                             seq_error = 0.01,
                             coi_method = c("variant", "frequency"),
                             solution_method = c("discrete", "continuous"),
                             use_bins = FALSE,
                             bin_size = 20,
                             replicates = 100,
                             parallel = FALSE,
                             ncpus = 8) {
  # Argument matching
  coi_method <- rlang::arg_match(coi_method)
  solution_method <- rlang::arg_match(solution_method)

  # Warning for deprecated arguments
  if (use_bins) {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "bootstrap_ci(use_bins)",
      details = "The ability to use bins to estimate the COI will be dropped in the next release."
    )
  }
  if (bin_size != 20) {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "bootstrap_ci(bin_size)",
      details = "The ability to use bins to estimate the COI will be dropped in the next release."
    )
  }

  # Define statistic function
  bootstrap_statistic <- function(data,
                                  indices,
                                  sim_data,
                                  solution_method,
                                  ...) {
    data <- sim_data
    data$data <- data$data[indices, ]

    # Determine function call based on solution_method
    switch(solution_method,
      discrete = compute_coi(data, "sim", ...)$coi,
      continuous = optimize_coi(data, "sim", ...)
    )
  }

  # Bootstrap
  boot_out <- boot::boot(
    data = seq_len(nrow(data$data)),
    statistic = bootstrap_statistic,
    R = replicates,
    parallel = if (parallel) "multicore" else "no",
    ncpus = ncpus,
    sim_data = data,
    solution_method = solution_method,
    max_coi = 25,
    seq_error = 0.01,
    coi_method = coi_method,
    use_bins = FALSE,
    bin_size = 20
  )

  tidy_boot_out <- tryCatch(
    broom::tidy(boot_out, conf.int = TRUE),
    error = function(e) {
      broom::tidy(boot_out, conf.int = FALSE) %>%
        tibble::add_column(conf.low = NaN, conf.high = NaN)
    }
  )

  tidy_boot_out %>%
    dplyr::rename(coi = "statistic") %>%
    tibble::add_column(estimates = list(boot_out$t), .after = "coi")
}
