check_freq_method <- function(wsmaf, plmaf, seq_error) {
  # Compute number of variant sites
  variant <- ifelse(wsmaf <= seq_error | wsmaf >= (1 - seq_error), 0, 1)
  n_variant <- sum(variant)

  # Compute expected number of variant sites using Hardy-Weinberg
  hardy_weinberg <- 2 * plmaf * (1 - plmaf)
  n_loci <- length(plmaf)

  # Find the 99% confidence interval
  bin_ci <- Hmisc::binconf(sum(hardy_weinberg), n_loci, alpha = 0.01) * n_loci

  # If the actual number of variant sites is less than the lower bound of the CI
  # return FALSE, otherwise return FALSE
  if (n_variant < bin_ci[1, "Lower"]) FALSE else TRUE
}
