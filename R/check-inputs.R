check_input_data <- function(data, data_type) {
  # Coverage value to be used when coverage column is missing
  mising_coverage <- 100

  if (data_type == "sim") {
    # Remove NA from our data frame
    data$data <- tidyr::drop_na(data$data)

    # Add coverage if it is missing
    if (!rlang::has_name(data$data, "coverage")) {
      data$data$coverage <- rep(mising_coverage, nrow(data$data))
    }
  } else if (data_type == "real") {
    # Remove NA from our data frame
    data <- tidyr::drop_na(data)

    # Add coverage if it is missing
    if (!rlang::has_name(data, "coverage")) {
      data$coverage <- rep(mising_coverage, nrow(data))
    }

    # Ensure we are calling the minor allele
    data <- check_minor_allele(data)
  }

  data
}

# In some cases we are fed in the major allele so we ensure we only examine
# the minor allele. If the PLAF is > 0.5, we know it is the major allele so
# we look at 1 - WSAF and 1 - PLAF.
check_minor_allele <- function(data) {
  dplyr::mutate(
    data,
    wsmaf = ifelse(.data$plmaf > 0.5, 1 - .data$wsmaf, .data$wsmaf),
    plmaf = ifelse(.data$plmaf > 0.5, 1 - .data$plmaf, .data$plmaf)
  )
}
