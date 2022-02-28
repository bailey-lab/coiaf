check_freq_method <- function(wsmaf, plmaf, seq_error) {
  # Compute number of variant sites
  variant <- ifelse(wsmaf <= seq_error | wsmaf >= (1 - seq_error), 0, 1)
  n_variant <- sum(variant)

  # Compute expected number of variant sites using Hardy-Weinberg
  hardy_weinberg <- 2 * plmaf * (1 - plmaf)
  n_loci <- length(plmaf)

  # Find the 99% confidence interval
  bin_ci <- Hmisc::binconf(sum(hardy_weinberg), n_loci, alpha = 0.01) * n_loci
  names(bin_ci) <- c("expected", "lower_ci", "uper_ci")

  # Return number of variant sites and CI info
  c(variant = n_variant, bin_ci)
}
